"""Garman and Kolhagen (1983) Currency options"""

from math import exp, log, pi, sqrt
from typing import Callable

from ..distributions import CND


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        rf: float,
        v: float,
        cdf: Callable[[float], float] = CND
) -> float:
    # Garman and Kolhagen (1983) Currency options
                
    d1 = (log(S / K) + (r - rf + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return S * exp(-rf * T) * cdf(d1) - K * exp(-r * T) * cdf(d2)
    else:
        return K * exp(-r * T) * cdf(-d2) - S * exp(-rf * T) * cdf(-d1)
