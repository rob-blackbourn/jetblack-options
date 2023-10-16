"""Barriers"""

from math import exp, sqrt

# Discrete barrier monitoring adjustment
def DiscreteAdjustedBarrier(S: float, H: float, v: float, dt: float) -> float:
    if H > S:
        return H * exp(0.5826 * v * sqrt(dt))
    elif H < S:
        return H * exp(-0.5826 * v * sqrt(dt))
    else:
        raise ValueError('no solution')
