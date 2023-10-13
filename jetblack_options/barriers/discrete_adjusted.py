"""Barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, cast

from ..distributions import CND, CBND
from ..european.black_scholes.plain_vanilla import GBlackScholes, EGBlackScholes, EGBlackScholes_OutPutFlag


# Discrete barrier monitoring adjustment
def DiscreteAdjustedBarrier(S: float, H: float, v: float, dt: float) -> float:
    if H > S:
        return H * exp(0.5826 * v * sqrt(dt))
    elif H < S:
        return H * exp(-0.5826 * v * sqrt(dt))
    else:
        raise ValueError('no solution')
