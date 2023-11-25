from __future__ import annotations

import enum
import functools
import struct
import typing

import numpy

from .utils import normalize, ArrayN, ArrayNx2, ArrayNx3, ArrayNx10
from ._unkerasify import KerasifyParser, Layer

if typing.TYPE_CHECKING:
    from Bio.PDB import Chain


DISTANCE_ALPHA_BETA = 1.5336
ALPHABET = numpy.array(list("ACDEFGHIKLMNPQRSTVWYX"))


def approximate_cb_position(
    ca: ArrayNx3[numpy.floating], 
    n: ArrayNx3[numpy.floating], 
    c: ArrayNx3[numpy.floating], 
    distance: float = DISTANCE_ALPHA_BETA
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
    return ca + v4 * distance


def create_residue_mask(
    ca: ArrayNx3[numpy.floating],
    n: ArrayNx3[numpy.floating],
    c: ArrayNx3[numpy.floating], 
) -> ArrayN[numpy.bool_]:
    assert ca.shape == n.shape
    assert ca.shape == c.shape
    mask = (
        numpy.isnan(ca).max(axis=-1)
        | numpy.isnan(n).max(axis=-1)
        | numpy.isnan(c).max(axis=-1)
    )
    return ~mask


def find_residue_partners(
    cb: ArrayNx3[numpy.floating], 
    mask: ArrayN[numpy.bool_]
) -> ArrayN[numpy.int64]:
    axes = (*range(len(cb.shape) - 2), -2, -1)
    # compute pairwise squared distance matrix
    r = numpy.sum(cb * cb, axis=-1).reshape(*cb.shape[:-1], 1)
    r[0] = r[-1] = numpy.nan
    D = r - 2 * cb @ cb.T + r.T
    # avoid selecting residue itself as the best
    D[numpy.diag_indices_from(D)] = numpy.inf
    # avoid selecting a masked residue as the best
    D[~mask, :] = D[:, ~mask] = numpy.inf
    # find closest other atom for each residue
    return numpy.nan_to_num(D, copy=False, nan=numpy.inf).argmin(axis=1)


def calc_conformation_descriptors(
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


class FeatureEncoder:
    @classmethod
    def load(cls):
        with open("foldseek/data/encoder_weights_3di.kerasify", "rb") as f:
            parser = KerasifyParser(f)
            layers = list(parser)
        return cls(layers)

    def __init__(self, layers: typing.Iterable[Layer]):
        self.layers = list(layers)

    def __call__(self, X: ArrayNx10[numpy.floating]) -> ArrayNx2[numpy.floating]:
        return functools.reduce(lambda x, f: f(x), self.layers, X)


class VirtualCenterCalculator:
    def __init__(self, alpha: float = 270.0, beta: float = 0.0, d: float = 2.0):
        self._alpha = numpy.deg2rad(alpha)
        self._beta = numpy.deg2rad(beta)
        self._d = d

    def __call__(
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


class CentroidEncoder:
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

    @classmethod
    def load(cls, invalid_state: int = 2):
        return cls(cls._CENTROIDS.copy(), invalid_state=invalid_state)

    def __init__(
        self, 
        centroids: ArrayNx2[numpy.float32], 
        invalid_state: int = 2
    ) -> None:
        self.invalid_state = invalid_state
        self.centroids = numpy.asarray(centroids)
        self.r2 = numpy.sum(self.centroids**2, 1).reshape(-1, 1).T

    def __call__(
        self, 
        embeddings: ArrayNx2[numpy.floating], 
        mask: ArrayN[numpy.bool_]
    ) -> ArrayN[numpy.uint8]:
        # compute pairwise squared distance matrix
        r1 = numpy.sum(embeddings * embeddings, 1).reshape(-1, 1)
        D = r1 - 2 * embeddings @ self.centroids.T + self.r2
        # find closest centroid
        states = numpy.empty(D.shape[0], dtype=numpy.uint8)
        D.argmin(axis=1, out=states)
        # use invalid state for masked residues
        states[~mask] = self.invalid_state
        return states


class Encoder:
    def __init__(
        self, 
        *, 
        alpha: float = 270.0, 
        beta: float = 0.0, 
        d: float = 2.0
    ) -> None:
        self.vc_calculator = VirtualCenterCalculator(alpha=alpha, beta=beta, d=d)
        self.feature_encoder = FeatureEncoder.load()
        self.centroid_encoder = CentroidEncoder.load()

    def encode_chain(
        self, 
        chain: Chain, 
        ca_residue: bool = True,
    ) -> ArrayN[numpy.uint8]:
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

    def encode_atoms(
        self, 
        ca: ArrayNx3[numpy.floating], 
        cb: ArrayNx3[numpy.floating], 
        n: ArrayNx3[numpy.floating], 
        c: ArrayNx3[numpy.floating],
    ) -> ArrayN[numpy.uint8]:
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
            cb_approx = approximate_cb_position(ca, n, c)
            cb[nan_indices] = cb_approx[nan_indices]

        # compute virtual center
        vc = self.vc_calculator(ca, cb, n)
        # mask residues without coordinates
        mask = create_residue_mask(ca, n, c)
        # find closest neighbor for each residue
        partner_index = find_residue_partners(vc, mask)
        # build position features from residue angles
        descriptors = calc_conformation_descriptors(ca, mask, partner_index)
        # compute embeddings and decode states
        embeddings = self.feature_encoder(descriptors)
        return self.centroid_encoder(embeddings, mask)


if __name__ == "__main__":
    pass
