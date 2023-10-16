"""Black 76"""

from math import exp, log, sqrt
from typing import Callable

from ..distributions import CDF


def price(
        is_call: bool,
        F: float,
        K: float,
        T: float,
        r: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    """Black (1976) Options on futures/forwards.

    Args:
        is_call (bool): True for a call, false for a put.
        F (float): The price of the future.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.
        cdf (Callable[[float], float], optional): The cumulative probability
            distribution function. Defaults to CDF.

    Returns:
        float: The option price.
    """

    d1 = (log(F / K) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return exp(-r * T) * (F * cdf(d1) - K * cdf(d2))
    else:
        return exp(-r * T) * (K * cdf(-d2) - F * cdf(-d1))
