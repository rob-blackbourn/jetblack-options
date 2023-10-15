from math import exp, log, sqrt

from ...distributions import CND


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float
) -> float:
    """Black-Scholes on a non-dividend paying stock.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.

    Returns:
        float: The price of the option.
    """
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return S * CND(d1) - K * exp(-r*T)* CND(d2)
    else:
        return K * exp(-r * T) * CND(-d2) - S * CND(-d1)