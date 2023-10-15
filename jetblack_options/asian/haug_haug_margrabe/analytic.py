"""Asian options"""

from math import exp, log, sqrt
from typing import Callable

from ...distributions import CND
from ...european.black_scholes.analytic import price as bs_price


def price(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CND
) -> float:

    # This is a modified version of the Levy formula, this is the formula published in "Asian Pyramid Power" By
    #  Haug, Haug and Margrabe

    h = (T - t1) / (n - 1)

    if b == 0:
        EA = S
    else:
        EA = S / n * exp(b * t1) * (1 - exp(b * h * n)) / (1 - exp(b * h))
   
    if m > 0:
        if SA > n / m * X: #  Exercise is certain for call, put must be out-of-the-money
        
            if not is_call:
                return 0
            else:
                SA = SA * m / n + EA * (n - m) / n
                return (SA - X) * exp(-r * T)

    if m == n - 1: # Only one fix left use Black-Scholes weighted with time
   
        X = n * X - (n - 1) * SA
        return bs_price(is_call, S, X, T, r, b, v, cdf=cdf) * 1 / n

    if b == 0:
        EA2 = (
            S * S * exp(v * v * t1) / (n * n) * 
            (
                (1 - exp(v * v * h * n)) / (1 - exp(v * v * h)) +
                2 / (1 - exp(v * v * h)) * (
                    n - (1 - exp(v * v * h * n)) / (1 - exp(v * v * h))
                )
            )
        )
    else:
        EA2 = (
            S * S * exp((2 * b + v * v) * t1) / (n * n) *
            (
                (1 - exp((2 * b + v * v) * h * n)) / (1 - exp((2 * b + v * v) * h)) +
                2 / (1 - exp((b + v * v) * h)) *
                (
                    (1 - exp(b * h * n)) / (1 - exp(b * h)) -
                    (1 - exp((2 * b + v * v) * h * n)) / (1 - exp((2 * b + v * v) * h))
                )
            )
        )

    vA = sqrt((log(EA2) - 2 * log(EA)) / T)

    if m > 0:
        X = n / (n - m) * X - m / (n - m) * SA
    
    d1 = (log(EA / X) + vA ** 2 / 2 * T) / (vA * sqrt(T))
    d2 = d1 - vA * sqrt(T)

    if is_call:
        value = exp(-r * T) * (EA * cdf(d1) - X * cdf(d2))
    else:
        value = exp(-r * T) * (X * cdf(-d2) - EA * cdf(-d1))

    return value * (n - m) / n
