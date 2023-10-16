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
        cdf: Callable[[float], float] = CDF
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

    Returns:
        float: The implied volatility.
    """

    vLow = 0.005
    vHigh = 4
    epsilon = 0.00000001
    cLow = price(is_call, S, K, T, r, b, vLow, cdf=cdf)
    cHigh = price(is_call, S, K, T, r, b, vHigh, cdf=cdf)
    N = 0
    vi = vLow + (p - cLow) * (vHigh - vLow) / (cHigh - cLow)
    while abs(p - price(is_call, S, K, T, r, b, vi, cdf=cdf)) > epsilon:
        N = N + 1
        if N > 20:
            break
        
        if price(is_call, S, K, T, r, b, vi, cdf=cdf) < p:
            vLow = vi
        else:
            vHigh = vi

        cLow = price(is_call, S, K, T, r, b, vLow, cdf=cdf)
        cHigh = price(is_call, S, K, T, r, b, vHigh, cdf=cdf)
        vi = vLow + (p - cLow) * (vHigh - vLow) / (cHigh - cLow)

    return vi
