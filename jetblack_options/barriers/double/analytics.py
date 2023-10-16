"""Double barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, Callable

from ...distributions import CDF
from ...european.black_scholes_merton import price as bs_price

# Double barrier options
def price(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    F = U * exp(delta1 * T)
    E = L * exp(delta2 * T)
    Sum1 = 0
    Sum2 = 0
    
    if is_call:
        for n in range(-5, 5+1):
            d1 = (log(S * U ** (2 * n) / (K * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d2 = (log(S * U ** (2 * n) / (F * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d3 = (log(L ** (2 * n + 2) / (K * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d4 = (log(L ** (2 * n + 2) / (F * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / v ** 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / v ** 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / v ** 2 + 1
            Sum1 = Sum1 + (U ** n / L ** n) ** mu1 * (L / S) ** mu2 * (cdf(d1) - cdf(d2)) - (L ** (n + 1) / (U ** n * S)) ** mu3 * (cdf(d3) - cdf(d4))
            Sum2 = Sum2 + (U ** n / L ** n) ** (mu1 - 2) * (L / S) ** mu2 * (cdf(d1 - v * sqrt(T)) - cdf(d2 - v * sqrt(T))) - (L ** (n + 1) / (U ** n * S)) ** (mu3 - 2) * (cdf(d3 - v * sqrt(T)) - cdf(d4 - v * sqrt(T)))
        out_value = S * exp((b - r) * T) * Sum1 - K * exp(-r * T) * Sum2
    else:
        for n in range(-5, 5+1):
            d1 = (log(S * U ** (2 * n) / (E * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d2 = (log(S * U ** (2 * n) / (K * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d3 = (log(L ** (2 * n + 2) / (E * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d4 = (log(L ** (2 * n + 2) / (K * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / v ** 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / v ** 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / v ** 2 + 1
            Sum1 = Sum1 + (U ** n / L ** n) ** mu1 * (L / S) ** mu2 * (cdf(d1) - cdf(d2)) - (L ** (n + 1) / (U ** n * S)) ** mu3 * (cdf(d3) - cdf(d4))
            Sum2 = Sum2 + (U ** n / L ** n) ** (mu1 - 2) * (L / S) ** mu2 * (cdf(d1 - v * sqrt(T)) - cdf(d2 - v * sqrt(T))) - (L ** (n + 1) / (U ** n * S)) ** (mu3 - 2) * (cdf(d3 - v * sqrt(T)) - cdf(d4 - v * sqrt(T)))

        out_value = K * exp(-r * T) * Sum2 - S * exp((b - r) * T) * Sum1

    if not is_in:
        return out_value
    else:
        return bs_price(is_call, S, K, T, r, b, v) - out_value
