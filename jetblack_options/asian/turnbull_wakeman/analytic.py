"""Turnbull-Wakeman analytic"""

from math import exp, log, sqrt
from typing import Callable

from ...distributions import CND
from ...european.black_scholes_merton import price as bs_price


def price(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CND
) -> float:
    """Arithmetic average rate option based on Turnbull-Wakeman.

    Args:
        is_call (bool): Call or put.
        S (float): Asset price.
        SA (float): Realized average so far.
        X (float): Strike price
        T (float): Time to maturity in years.
        T2 (float): Original time in average period in years, constant over life
            of option.
        r (float): Risk free rate.
        b (float): Cost of carry underlying asset can be positive and negative.
        v (float): Annualized volatility of asset price.
        cnd (Callable[[float], float], optional): The cumulative normal
            distribution function. Defaults to CND.

    Raises:
        ValueError: _description_

    Returns:
        float: The price of the option.
    """

    # Time to start of average period in years
    t1 = max(0, T - T2)
    # Remaining time of average periods
    tau = T2 - T
   
    if b == 0:    
        M1 = 1
    else:
        M1 = (exp(b * T) - exp(b * t1)) / (b * (T - t1))

    # Take into account when option wil be exercised
    if tau > 0:
    
        if T2 / T * K - tau / T * SA < 0:
    
            if is_call:
                SA = SA * (T2 - T) / T2 + S * M1 * T / T2 # Expected average at maturity
                return max(0, SA - K) * exp(-r * T)
            else:
                return 0

    if b == 0: # Extended to hold for options on futures 16 May 1999 Espen Haug
       M2 = (
           2 * exp(v * v * T) / (v ** 4 * (T - t1) ** 2) -
           2 * exp(v * v * t1) * (1 + v * v * (T - t1)) / (v ** 4 * (T - t1) ** 2)
        )
    else:
        M2 = (
            2 * exp((2 * b + v * v) * T) / ((b + v * v) * (2 * b + v * v) * (T - t1) ** 2) +
            2 * exp((2 * b + v * v) * t1) / (b * (T - t1) ** 2) * (1 / (2 * b + v * v) - exp(b * (T - t1)) / (b + v * v))
        )

    bA = log(M1) / T
    vA = sqrt(log(M2) / T - 2 * bA)
    if tau > 0:
        K = T2 / T * K - tau / T * SA
        return bs_price(is_call, S, K, T, r, bA, vA, cdf=cdf) * T / T2
    else:
        return bs_price(is_call, S, K, T, r, bA, vA, cdf=cdf)
