"""Double barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, Callable

from ...distributions import CND
from ...european.black_scholes.analytic import price as bs_price

# Double barrier options
def price(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    
    F = U * exp(delta1 * T)
    E = L * exp(delta2 * T)
    Sum1 = 0
    Sum2 = 0
    
    if is_call:
        for n in range(-5, 5+1):
            d1 = (log(S * U ** (2 * n) / (X * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d2 = (log(S * U ** (2 * n) / (F * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d3 = (log(L ** (2 * n + 2) / (X * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d4 = (log(L ** (2 * n + 2) / (F * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / v ** 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / v ** 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / v ** 2 + 1
            Sum1 = Sum1 + (U ** n / L ** n) ** mu1 * (L / S) ** mu2 * (cnd(d1) - cnd(d2)) - (L ** (n + 1) / (U ** n * S)) ** mu3 * (cnd(d3) - cnd(d4))
            Sum2 = Sum2 + (U ** n / L ** n) ** (mu1 - 2) * (L / S) ** mu2 * (cnd(d1 - v * sqrt(T)) - cnd(d2 - v * sqrt(T))) - (L ** (n + 1) / (U ** n * S)) ** (mu3 - 2) * (cnd(d3 - v * sqrt(T)) - cnd(d4 - v * sqrt(T)))
        out_value = S * exp((b - r) * T) * Sum1 - X * exp(-r * T) * Sum2
    else:
        for n in range(-5, 5+1):
            d1 = (log(S * U ** (2 * n) / (E * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d2 = (log(S * U ** (2 * n) / (X * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d3 = (log(L ** (2 * n + 2) / (E * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d4 = (log(L ** (2 * n + 2) / (X * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / v ** 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / v ** 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / v ** 2 + 1
            Sum1 = Sum1 + (U ** n / L ** n) ** mu1 * (L / S) ** mu2 * (cnd(d1) - cnd(d2)) - (L ** (n + 1) / (U ** n * S)) ** mu3 * (cnd(d3) - cnd(d4))
            Sum2 = Sum2 + (U ** n / L ** n) ** (mu1 - 2) * (L / S) ** mu2 * (cnd(d1 - v * sqrt(T)) - cnd(d2 - v * sqrt(T))) - (L ** (n + 1) / (U ** n * S)) ** (mu3 - 2) * (cnd(d3 - v * sqrt(T)) - cnd(d4 - v * sqrt(T)))

        out_value = X * exp(-r * T) * Sum2 - S * exp((b - r) * T) * Sum1

    if not is_in:
        return out_value
    else:
        return bs_price(is_call, S, X, T, r, b, v) - out_value
