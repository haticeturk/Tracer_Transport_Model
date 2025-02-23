"""This module contains utility functions for working with Pastas models."""

import logging
from datetime import datetime, timedelta
from logging import handlers

import numpy as np
from pandas import Series, Timedelta, Timestamp, date_range, to_datetime
from pandas.tseries.frequencies import to_offset
from scipy import interpolate

logger = logging.getLogger(__name__)
def check_numba():
    try:
        from numba import njit
    except ImportError:
        logger.warning("Numba is not installed. Installing Numba is "
                       "recommended for significant speed-ups.")
