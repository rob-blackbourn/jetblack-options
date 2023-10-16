"""Garman and Kolhagen (1983) Currency options"""

from math import exp, log, sqrt
from typing import Callable

from ..distributions import CDF


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        rf: float,
        v: float,
        cdf: Callable[[float], float] = CDF
) -> float:
    """Garman and Kolhagen (1983) Currency options.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate of the base currency.
        rf (float): The risk free rate of the quote currency.
        v (float): The asset volatility.
        cdf (Callable[[float], float], optional): The cumulative probability
            distribution function. Defaults to CDF.

    Returns:
        float: _description_
    """
    # Garman and Kolhagen (1983) Currency options
                
    d1 = (log(S / K) + (r - rf + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return S * exp(-rf * T) * cdf(d1) - K * exp(-r * T) * cdf(d2)
    else:
        return K * exp(-r * T) * cdf(-d2) - S * exp(-rf * T) * cdf(-d1)
