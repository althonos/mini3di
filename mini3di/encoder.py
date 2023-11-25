from __future__ import annotations

import abc
import enum
import functools
import struct
import typing

import numpy
import numpy.ma

from . import _unkerasify
from .layers import Layer, CentroidLayer, Model
from .utils import normalize

try:
    from importlib.resources import files as resource_files
except ImportError:
    from importlib_resources import files as resource_files

T = typing.TypeVar("T")
if typing.TYPE_CHECKING:
    from Bio.PDB import Chain
    from .utils import ArrayN, ArrayNx2, ArrayNx3, ArrayNx10, ArrayNxM

DISTANCE_ALPHA_BETA = 1.5336
ALPHABET = numpy.array(list("ACDEFGHIKLMNPQRSTVWYX"))


class _BaseEncoder(abc.ABC, typing.Generic[T]):
    @abc.abstractmethod
    def encode_atoms(
        self, 
        ca: ArrayNx3[numpy.floating], 
        cb: ArrayNx3[numpy.floating], 
        n: ArrayNx3[numpy.floating], 
        c: ArrayNx3[numpy.floating],
    ) -> T:
        raise NotImplementedError
    
    def encode_chain(
        self, 
        chain: Chain, 
        ca_residue: bool = True,
    ) -> T:
        # extract residues
        if ca_residue:
            residues = [residue for residue in chain.get_residues() if "CA" in residue]
        else:
            residues = chain.get_residues()
        # extract atom coordinates
        r = len(residues)
        ca = numpy.array(numpy.nan, dtype=numpy.float32).repeat(3 * r).reshape(r, 3)
        cb = ca.copy()
        n = ca.copy()
        c = ca.copy()
        for i, residue in enumerate(residues):
            ca[i, :] = residue["CA"].coord
            n[i, :] = residue["N"].coord
            c[i, :] = residue["C"].coord
            if "CB" in residue:
                cb[i, :] = residue["CB"].coord
        # encoder coordiantes
        return self.encode_atoms(ca, cb, n, c)


class VirtualCenterEncoder(_BaseEncoder["ArrayNx3[numpy.float32]"]):
   
    def __init__(
        self, 
        *, 
        alpha: float = 270.0, 
        beta: float = 0.0, 
        d: float = 2.0,
        distance_alpha_beta = DISTANCE_ALPHA_BETA,
    ) -> None:
        self._alpha = numpy.deg2rad(alpha)
        self._beta = numpy.deg2rad(beta)
        self._d = d
        self.distance_alpha_beta = distance_alpha_beta
    
    def _compute_virtual_center(
        self,
        ca: ArrayNx3[numpy.floating], 
        cb: ArrayNx3[numpy.floating], 
        n: ArrayNx3[numpy.floating], 
    ) -> ArrayNx3[numpy.floating]:
        assert ca.shape == n.shape
        assert ca.shape == cb.shape
        v = cb - ca
        a = cb - ca
        b = n - ca
        # normal angle
        k = normalize(numpy.cross(a, b, axis=-1))
        v = (
            v * numpy.cos(self._alpha)
            + numpy.cross(k, v) * numpy.sin(self._alpha)
            + k * (k * v).sum(axis=-1).reshape(-1, 1) * (1 - numpy.cos(self._alpha))
        )
        # dihedral angle
        k = normalize(n - ca)
        v = (
            v * numpy.cos(self._beta)
            + numpy.cross(k, v) * numpy.sin(self._beta)
            + k * (k * v).sum(axis=-1).reshape(-1, 1) * (1 - numpy.cos(self._beta))
        )
        return ca + v * self._d

    def _approximate_cb_position(
        self,
        ca: ArrayNx3[numpy.floating], 
        n: ArrayNx3[numpy.floating], 
        c: ArrayNx3[numpy.floating], 
    ) -> ArrayNx3[numpy.floating]:
        """Approximate the position of the Cβ from the backbone atoms."""
        assert ca.shape == n.shape
        assert ca.shape == c.shape
        v1 = normalize(c - ca)
        v2 = normalize(n - ca)
        v3 = v1 / 3.0
        b1 = v2 + v3
        b2 = numpy.cross(v1, b1, axis=-1)
        u1 = normalize(b1)
        u2 = normalize(b2)
        v4 = (numpy.sqrt(8) / 3.0) * ((-u1 / 2.0) - (u2 * numpy.sqrt(3) / 2.0)) - v3
        return ca + v4 * self.distance_alpha_beta

    def _create_nan_mask(
        self,
        *arrays: ArrayNxM[numpy.floating],
    ) -> ArrayN[numpy.bool_]:
        maximum = functools.partial(numpy.max, axis=-1)
        return ~functools.reduce( # type: ignore
            numpy.bitwise_or, map(maximum, map( numpy.isnan, arrays))
        ) 

    def encode_atoms(
        self, 
        ca: ArrayNx3[numpy.floating], 
        cb: ArrayNx3[numpy.floating], 
        n: ArrayNx3[numpy.floating], 
        c: ArrayNx3[numpy.floating],
    ) -> ArrayNx3[numpy.float32]:
        ca = numpy.asarray(ca)
        cb = numpy.asarray(cb)
        n = numpy.asarray(n)
        c = numpy.asarray(c)

        assert ca.shape == cb.shape
        assert ca.shape == c.shape
        assert ca.shape == n.shape

        # fix CB positions if needed
        nan_indices = numpy.isnan(cb)
        if numpy.any(nan_indices):
            cb_approx = self._approximate_cb_position(ca, n, c)
            # avoid writing to CB directly since it should be callee-save
            cb_approx[~nan_indices] = cb[~nan_indices]
            cb = cb_approx
        # compute virtual center
        vc = self._compute_virtual_center(ca, cb, n)
        # mask residues without coordinates
        mask = self._create_nan_mask(ca, n, c)
        return numpy.ma.masked_array(  # type: ignore
            vc, 
            mask=~mask.repeat(vc.shape[1]).reshape(vc.shape), 
            fill_value=numpy.nan,
        )  


class FeatureEncoder(_BaseEncoder["ArrayN[numpy.float32]"]):

    def __init__(self) -> None:
        self.vc_encoder = VirtualCenterEncoder()

    def _calc_conformation_descriptors(
        self,
        ca: ArrayNx3[numpy.floating], 
        mask: ArrayN[numpy.bool_], 
        partner_index: ArrayN[numpy.int64],
        dtype: typing.Type[numpy.floating] = numpy.float32,
    ) -> ArrayNx10[numpy.floating]:
        # build arrays of indices to use for vectorized angles
        n = ca.shape[0]
        I = numpy.arange(1, ca.shape[-2] - 1)
        J = partner_index[I]
        # compute conformational descriptors
        u1 = normalize(ca[..., I, :] - ca[..., I - 1, :])
        u2 = normalize(ca[..., I + 1, :] - ca[..., I, :])
        u3 = normalize(ca[..., J, :] - ca[..., J - 1, :])
        u4 = normalize(ca[..., J + 1, :] - ca[..., J, :])
        u5 = normalize(ca[..., J, :] - ca[..., I, :])
        desc = numpy.zeros((*ca.shape[:-1], 10), dtype=dtype)
        desc[I, 0] = numpy.sum(u1 * u2, axis=-1)
        desc[I, 1] = numpy.sum(u3 * u4, axis=-1)
        desc[I, 2] = numpy.sum(u1 * u5, axis=-1)
        desc[I, 3] = numpy.sum(u3 * u5, axis=-1)
        desc[I, 4] = numpy.sum(u1 * u4, axis=-1)
        desc[I, 5] = numpy.sum(u2 * u3, axis=-1)
        desc[I, 6] = numpy.sum(u1 * u3, axis=-1)
        desc[I, 7] = numpy.linalg.norm(ca[I] - ca[J], axis=-1)
        desc[I, 8] = numpy.clip(J - I, -4, 4)
        desc[I, 9] = numpy.copysign(numpy.log(numpy.abs(J - I) + 1), J - I)
        # update mask around invalid positions
        mask[1:-1] &= (
            mask[I - 1] & mask[I] & mask[I + 1] & mask[J - 1] & mask[J] & mask[J + 1]
        )
        mask[0] = mask[n - 1] = False
        return desc

    def _find_residue_partners(
        self,
        x: ArrayNx3[numpy.floating], 
        mask: ArrayN[numpy.bool_]
    ) -> ArrayN[numpy.int64]:
        axes = (*range(len(x.shape) - 2), -2, -1)
        # compute pairwise squared distance matrix
        r = numpy.sum(x * x, axis=-1).reshape(*x.shape[:-1], 1)
        r[0] = r[-1] = numpy.nan
        D = r - 2 * x @ x.T + r.T
        # avoid selecting residue itself as the best
        D[numpy.diag_indices_from(D)] = numpy.inf
        # avoid selecting a masked residue as the best
        D[~mask, :] = D[:, ~mask] = numpy.inf
        # find closest other atom for each residue
        return numpy.nan_to_num(D, copy=False, nan=numpy.inf).argmin(axis=1)

    def encode_atoms(
        self, 
        ca: ArrayNx3[numpy.floating], 
        cb: ArrayNx3[numpy.floating], 
        n: ArrayNx3[numpy.floating], 
        c: ArrayNx3[numpy.floating],
    ) -> ArrayN[numpy.uint8]:
        vc = self.vc_encoder.encode_atoms(ca, cb, n, c)
        mask = ~vc.mask[:, 0]
        # find closest neighbor for each residue
        partner_index = self._find_residue_partners(vc.data, mask)
        # build position features from residue angles
        descriptors = self._calc_conformation_descriptors(ca, mask, partner_index)
        return numpy.ma.masked_array(  # type: ignore
            descriptors, 
            mask=~mask.repeat(descriptors.shape[1]).reshape(descriptors.shape), 
            fill_value=numpy.nan,
        )  


class Encoder(_BaseEncoder["ArrayN[numpy.uint8]"]):

    _INVALID_STATE = 2
    _CENTROIDS: ArrayNx2[numpy.float32] = numpy.array(
        [
            [-1.0729, -0.3600],
            [-0.1356, -1.8914],
            [0.4948, -0.4205],
            [-0.9874, 0.8128],
            [-1.6621, -0.4259],
            [2.1394, 0.0486],
            [1.5558, -0.1503],
            [2.9179, 1.1437],
            [-2.8814, 0.9956],
            [-1.1400, -2.0068],
            [3.2025, 1.7356],
            [1.7769, -1.3037],
            [0.6901, -1.2554],
            [-1.1061, -1.3397],
            [2.1495, -0.8030],
            [2.3060, -1.4988],
            [2.5522, 0.6046],
            [0.7786, -2.1660],
            [-2.3030, 0.3813],
            [1.0290, 0.8772],
        ]
    )

    def __init__(self) -> None:
        self.feature_encoder = FeatureEncoder()
        with resource_files(__package__).joinpath("encoder_weights_3di.kerasify").open("rb") as f:
            layers = _unkerasify.load(f)
            layers.append(CentroidLayer(self._CENTROIDS))
        self.vae_encoder = Model(layers)

    def encode_atoms(
        self, 
        ca: ArrayNx3[numpy.floating], 
        cb: ArrayNx3[numpy.floating], 
        n: ArrayNx3[numpy.floating], 
        c: ArrayNx3[numpy.floating],
    ) -> ArrayN[numpy.uint8]:
        descriptors = self.feature_encoder.encode_atoms(ca, cb, n, c)
        states = self.vae_encoder(descriptors.data)
        states[descriptors.mask[:, 0]] = self._INVALID_STATE
        return numpy.ma.masked_array(
            states, 
            mask=descriptors.mask[:, 0], 
            fill_value=self._INVALID_STATE,
        )

    def build_sequence(self, states: ArrayN[numpy.uint8]) -> str:
        return "".join( ALPHABET[states] )

