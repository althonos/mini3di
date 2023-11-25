from typing import Annotated, Literal, TypeVar

import numpy
import numpy.typing

DType = TypeVar("DType", bound=numpy.generic)
ArrayN = Annotated[numpy.typing.NDArray[DType], Literal["N"]]
ArrayNx2 = Annotated[numpy.typing.NDArray[DType], Literal["N", 2]]
ArrayNx3 = Annotated[numpy.typing.NDArray[DType], Literal["N", 3]]
ArrayNx10 = Annotated[numpy.typing.NDArray[DType], Literal["N", 10]]
ArrayNxM = Annotated[numpy.typing.NDArray[DType], Literal["N", "M"]]

def normalize(x, *, inplace=False):
    norm = numpy.linalg.norm(x, axis=-1).reshape(*x.shape[:-1], 1)
    return numpy.divide(x, norm, out=x if inplace else None)

def relu(
    x,
    /,
    out=None,
    *,
    where=True,
    casting="same_kind",
    order="K",
    dtype=None,
    subok=True,
):
    return numpy.maximum(
        0.0,
        x,
        out=out,
        where=where,
        casting=casting,
        order=order,
        dtype=dtype,
        subok=subok,
    )
