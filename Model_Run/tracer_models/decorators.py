from functools import wraps
from logging import getLogger

logger = getLogger(__name__)
def njit(function):
    try:
        from numba import njit as jit
        return jit(function)
    except ImportError:
        return function
