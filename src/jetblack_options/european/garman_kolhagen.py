"""Garman and Kolhagen (1983) Currency options"""

from math import exp, log, sqrt
from statistics import NormalDist

norm = NormalDist()
cdf = norm.cdf
pdf = norm.pdf
inv_cdf = norm.inv_cdf


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        rf: float,
        v: float,
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
