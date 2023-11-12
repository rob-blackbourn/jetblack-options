"""Black Scholes variance analytic solutions"""

from math import exp, log, sqrt
from statistics import NormalDist
from typing import Callable

from ..distributions import CDF

norm = NormalDist()
cdf = norm.cdf


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
) -> float:
    # The generalized Black and Scholes formula on variance form

    d1 = (log(S / K) + (b + v / 2) * T) / sqrt(v * T)
    d2 = d1 - sqrt(v * T)

    if is_call:
        return S * exp((b - r) * T) * cdf(d1) - K * exp(-r * T) * cdf(d2)
    else:
        return K * exp(-r * T) * cdf(-d2) - S * exp((b - r) * T) * cdf(-d1)
