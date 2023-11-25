from __future__ import annotations

import abc
import enum
import functools
import itertools
import struct
import typing
from typing import BinaryIO, Iterable

import numpy

from .utils import relu

if typing.TYPE_CHECKING:
    from .utils import ArrayNxM


class LayerType(enum.IntEnum):
    DENSE = 1
    CONVOLUTION2D = 2
    FLATTEN = 3
    ELU = 4
    ACTIVATION = 5
    MAXPOOLING2D = 6
    LSTM = 7
    EMBEDDING = 8


class ActivationType(enum.IntEnum):
    LINEAR = 1
    RELU = 2
    SOFTPLUS = 3
    SIGMOID = 4
    TANH = 5
    HARD_SIGMOID = 6


class Layer(abc.ABC):
    @abc.abstractmethod
    def __call__(self, X: ArrayNxM[numpy.floating]) -> ArrayNxM[numpy.floating]:
        raise NotImplementedError


class DenseLayer(Layer):
    def __init__(self, weights, biases=None, activation=ActivationType.RELU):
        self.activation = activation
        self.weights = numpy.asarray(weights)
        if biases is None:
            self.biases = numpy.zeros(self.weights.shape[1])
        else:
            self.biases = numpy.asarray(biases)

    def __call__(self, X: ArrayNxM[numpy.floating]) -> ArrayNxM[numpy.floating]:
        _X = numpy.asarray(X)
        out = _X @ self.weights
        out += self.biases

        if self.activation == ActivationType.RELU:
            return relu(out, out=out)
        else:
            return out


class Model:

    @classmethod
    def load(cls, f: BinaryIO) -> Model:
        parser = KerasifyParser(f)
        return cls(parser)

    def __init__(self, layers: Iterable[Layer] = ()):
        self.layers = list(layers)
    
    def __call__(self, X: ArrayNxM[numpy.floating]) -> ArrayNxM[numpy.floating]:
        return functools.reduce(lambda x, f: f(x), self.layers, X)


class KerasifyParser:
    """An incomplete parser for model files serialized with `kerasify`.

    Note:
        Only dense layers are supported, since the ``foldseek`` VQ-VAE model
        is only using 3 dense layers.

    """

    def __init__(self, file: typing.BinaryIO) -> None:
        self.file = file
        self.buffer = bytearray(1024)
        (self.n_layers,) = self._get("I")

    def __iter__(self) -> KerasifyParser:
        return self

    def __next__(self) -> Layer:
        layer = self.read()
        if layer is None:
            raise StopIteration
        return layer

    def _read(self, format: str) -> memoryview:
        n = struct.calcsize(format)
        if len(self.buffer) < n:
            self.buffer.extend(
                itertools.islice(itertools.repeat(0), n - len(self.buffer))
            )
        v = memoryview(self.buffer)[:n]
        self.file.readinto(v)  # type: ignore
        return v

    def _get(self, format: str) -> typing.Tuple[typing.Any, ...]:
        v = self._read(format)
        return struct.unpack(format, v)

    def read(self) -> typing.Optional[Layer]:
        if self.n_layers == 0:
            return None

        self.n_layers -= 1
        layer_type = LayerType(self._get("I")[0])
        if layer_type == LayerType.DENSE:
            (w0,) = self._get("I")
            (w1,) = self._get("I")
            (b0,) = self._get("I")
            weights = (
                numpy.frombuffer(self._read(f"={w0*w1}f"), dtype="f4")
                .reshape(w0, w1)
                .copy()
            )
            biases = numpy.frombuffer(self._read(f"={b0}f"), dtype="f4").copy()
            activation = ActivationType(self._get("I")[0])
            return DenseLayer(weights, biases, activation)
        else:
            raise NotImplementedError(f"Unsupported layer type: {layer_type!r}")
