"""implied volatility for Black-Scholes-Merton"""

from typing import Callable

from ..distributions import CDF

from .black_scholes_merton import price

def ivol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        p: float,
        *,
        cdf: Callable[[float], float] = CDF,
        max_iterations: int = 20,
        epsilon = 1e-8
) -> float:
    """Calculate the volatility of an option that is implied by the price.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        p (float): The option price.
        cdf (Callable[[float], float], optional): The cumulative probability
            distribution function. Defaults to CDF.
        max_iterations (int, Optional): The maximum number of iterations before
            a price is returned. Defaults to 20.
        epsilon (float, Optional): The largest acceptable error. Defaults to 1e-8.

    Returns:
        float: The implied volatility.
    """

    v_lo = 0.005
    v_hi = 4
    p_lo = price(is_call, S, K, T, r, b, v_lo, cdf=cdf)
    p_hi = price(is_call, S, K, T, r, b, v_hi, cdf=cdf)

    n = 0
    v = v_lo + (p - p_lo) * (v_hi - v_lo) / (p_hi - p_lo)
    p1 = price(is_call, S, K, T, r, b, v, cdf=cdf)
    while abs(p - p1) > epsilon and n < max_iterations:
        n += 1
        
        if p1 < p:
            v_lo = v
        else:
            v_hi = v

        p_lo = price(is_call, S, K, T, r, b, v_lo, cdf=cdf)
        p_hi = price(is_call, S, K, T, r, b, v_hi, cdf=cdf)
        v = v_lo + (p - p_lo) * (v_hi - v_lo) / (p_hi - p_lo)
        p1 = price(is_call, S, K, T, r, b, v, cdf=cdf)

    return v
