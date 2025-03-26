"""
Gene regulatory network inference from single-cell data.

Executable gene regulatory network (GRN) inference method adapted to
time-stamped scRNA-seq datasets. The algorithm consists in calibrating
the parameters of a mechanistic model of gene expression, which can then
be simulated to reproduce the dataset used for inference.
"""
from importlib.metadata import version as _version
from cardamom.model import NetworkModel

__all__ = [
    'NetworkModel',
]

try:
    __version__ = _version('cardamom')
except Exception:
    __version__ = 'unknown version'
