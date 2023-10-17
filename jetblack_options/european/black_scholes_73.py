"""Black-Scholes 1973"""

from math import exp, log, sqrt
from typing import Callable

from ..distributions import CDF
from ..implied_volatility import solve_ivol


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    """Black-Scholes for a non-dividend paying stock.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        v (float): The asset volatility.
        cdf (Callable[[float], float], optional): The cumulative probability
            distribution function. Defaults to CDF.

    Returns:
        float: The price of the option.
    """
    d1 = (log(S / K) + (r + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return S * cdf(d1) - K * exp(-r*T) * cdf(d2)
    else:
        return K * exp(-r * T) * cdf(-d2) - S * cdf(-d1)


def ivol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        p: float,
        *,
        cdf: Callable[[float], float] = CDF,
        max_iterations: int = 35,
        epsilon=1e-8
) -> float:
    """Calculate the volatility of a Black-Scholes 73 option that is implied by
    the price.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        p (float): The option price.
        cdf (Callable[[float], float], optional): The cumulative probability
            distribution function. Defaults to CDF.
        max_iterations (int, Optional): The maximum number of iterations before
            a price is returned. Defaults to 35.
        epsilon (float, Optional): The largest acceptable error. Defaults to 1e-8.

    Returns:
        float: The implied volatility.
    """
    return solve_ivol(
        p,
        lambda v: price(is_call, S, K, T, r, v, cdf=cdf),
        max_iterations=max_iterations,
        epsilon=epsilon
    )
