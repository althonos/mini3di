import enum
import functools
import struct
import importlib.resources

import numpy

from ._unkerasify import KerasifyParser
from .utils import normalize


DISTANCE_ALPHA_BETA = 1.5336
ALPHABET = numpy.array(list("ACDEFGHIKLMNPQRSTVWYX"))


def approximate_cb_position(ca, n, c, distance=DISTANCE_ALPHA_BETA):
    """Approximate the position of the Cβ from the backbone atoms.
    """
    assert ca.shape == n.shape
    assert ca.shape == c.shape
    v1 = normalize(c - ca)
    v2 = normalize(n - ca)
    v3 = v1 / 3.0
    b1 = v2 + v3
    b2 = numpy.cross(v1, b1, axis=-1)
    u1 = normalize(b1)
    u2 = normalize(b2)
    v4 = (numpy.sqrt(8)/3.0)*((-u1/2.0) - (u2*numpy.sqrt(3)/2.0)) - v3
    return ca + v4 * distance


def create_residue_mask(ca, n, c):
    """Create a mask for residues missing backbone coordinates.
    """
    assert ca.shape == n.shape
    assert ca.shape == c.shape
    mask = numpy.isnan(ca).max(axis=-1) | numpy.isnan(n).max(axis=-1) | numpy.isnan(c).max(axis=-1)
    return ~mask


def find_residue_partners(cb, mask):
    """Find the closest other residue for each residue.
    """
    assert cb.shape[0] == mask.shape[0]
    assert cb.shape[1] == 3
    # compute pairwise squared distance matrix
    r = numpy.sum(cb*cb, axis=-1).reshape(-1, 1)
    r[0] = r[-1] = numpy.nan
    D = r - 2*cb@cb.T + r.T
    # avoid selecting residue itself as the best
    D[numpy.diag_indices_from(D)] = numpy.inf
    # avoid selecting a masked residue as the best
    D[~mask, :] = D[:, ~mask] = numpy.inf
    # find closest other atom for each residue
    return numpy.nan_to_num(D, copy=False, nan=numpy.inf).argmin(axis=1)


def calc_conformation_descriptors(ca, mask, partner_index):
    """Compute conformation descriptors for each residues.
    """
    # build arrays of indices to use for vectorized angles
    n = ca.shape[0]
    I = numpy.arange(1, ca.shape[-2] - 1)
    J = partner_index[I]
    # compute features 
    u1 = normalize(ca[..., I, :] - ca[..., I-1, :])
    u2 = normalize(ca[..., I+1, :] - ca[..., I, :])
    u3 = normalize(ca[..., J, :] - ca[..., J-1, :])
    u4 = normalize(ca[..., J+1, :] - ca[..., J, :])
    u5 = normalize(ca[..., J, :] - ca[..., I, :])
    features = numpy.zeros((*ca.shape[:-1], 10))
    features[I, 0] = numpy.sum(u1 * u2, axis=-1)
    features[I, 1] = numpy.sum(u3 * u4, axis=-1)
    features[I, 2] = numpy.sum(u1 * u5, axis=-1)
    features[I, 3] = numpy.sum(u3 * u5, axis=-1)
    features[I, 4] = numpy.sum(u1 * u4, axis=-1)
    features[I, 5] = numpy.sum(u2 * u3, axis=-1)
    features[I, 6] = numpy.sum(u1 * u3, axis=-1)
    features[I, 7] = numpy.linalg.norm(ca[I] - ca[J], axis=-1)
    features[I, 8] = numpy.clip(J - I, -4, 4)
    features[I, 9] = numpy.copysign(numpy.log(numpy.abs(J - I) + 1), J - I)
    # update mask around invalid positions
    mask[1:-1] &= mask[I-1] & mask[I] & mask[I+1] & mask[J-1] & mask[J] & mask[J+1]
    mask[0] = mask[n-1] = False
    return features


class FeatureEncoder:

    @classmethod
    def load(cls):
        with importlib.resources.files(__package__).joinpath("encoder_weights_3di.kerasify").open("rb") as f:
            parser = KerasifyParser(f)
            layers = list(parser)
        return cls(layers)

    def __init__(self, layers):
        self.layers = layers

    def __call__(self, X):
        return functools.reduce(lambda x, f: f(x), self.layers, X)


class VirtualCenterCalculator:

    def __init__(self, alpha=270, beta=0, d=2):
        self._alpha = numpy.deg2rad(alpha)
        self._beta = numpy.deg2rad(beta)
        self._d = d

    def __call__(self, ca, cb, n):
        assert ca.shape == n.shape
        assert ca.shape == cb.shape
        v = cb - ca
        a = cb - ca
        b = n - ca
        # normal angle
        k = normalize(numpy.cross(a, b, axis=-1))
        v = v * numpy.cos(self._alpha) + numpy.cross(k, v) * numpy.sin(self._alpha) + k * (k * v).sum(axis=-1).reshape(-1, 1) * (1 - numpy.cos(self._alpha))
        # dihedral angle
        k = normalize(n - ca)
        v = v * numpy.cos(self._beta) + numpy.cross(k, v) * numpy.sin(self._beta) + k * (k * v).sum(axis=-1).reshape(-1, 1) * (1 - numpy.cos(self._beta))
        return ca + v * self._d


class CentroidEncoder:

    _CENTROIDS = numpy.array([
        [ -1.0729,  -0.3600],
        [ -0.1356,  -1.8914],
        [  0.4948,  -0.4205],
        [ -0.9874,   0.8128],
        [ -1.6621,  -0.4259],
        [  2.1394,   0.0486],
        [  1.5558,  -0.1503],
        [  2.9179,   1.1437],
        [ -2.8814,   0.9956],
        [ -1.1400,  -2.0068],
        [  3.2025,   1.7356],
        [  1.7769,  -1.3037],
        [  0.6901,  -1.2554],
        [ -1.1061,  -1.3397],
        [  2.1495,  -0.8030],
        [  2.3060,  -1.4988],
        [  2.5522,   0.6046],
        [  0.7786,  -2.1660],
        [ -2.3030,   0.3813],
        [  1.0290,   0.8772],
    ])

    def __init__(self, centroids=None, invalid_state=2):
        self.invalid_state = invalid_state
        self.centroids = self._CENTROIDS.copy() if centroids is None else numpy.asarray(centroids)
        self.r2 = numpy.sum(self.centroids**2, 1).reshape(-1, 1).T

    def __call__(self, embeddings, mask):
        # compute pairwise squared distance matrix
        r1 = numpy.sum(embeddings*embeddings, 1).reshape(-1, 1)
        D = r1 - 2*embeddings@self.centroids.T + self.r2
        # find closest centroid
        states = D.argmin(axis=1)
        # use invalid state for masked residues
        states[~mask] = self.invalid_state
        return states






if __name__ == "__main__":
    pass