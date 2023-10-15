"""Soft Barrier"""

from math import exp, log, sqrt
from typing import Callable

from ...distributions import CND
from ...european.black_scholes.analytic import price as bs_price

# Soft barrier options
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
        *,
        cnd: Callable[[float], float] = CND
) -> float:

    if is_call:
        eta = 1
    else:
        eta = -1
    
    mu = (b + v ** 2 / 2) / v ** 2
    lambda1 = exp(-1 / 2 * v ** 2 * T * (mu + 0.5) * (mu - 0.5))
    lambda2 = exp(-1 / 2 * v ** 2 * T * (mu - 0.5) * (mu - 1.5))
    d1 = log(U ** 2 / (S * X)) / (v * sqrt(T)) + mu * v * sqrt(T)
    d2 = d1 - (mu + 0.5) * v * sqrt(T)
    d3 = log(U ** 2 / (S * X)) / (v * sqrt(T)) + (mu - 1) * v * sqrt(T)
    d4 = d3 - (mu - 0.5) * v * sqrt(T)
    e1 = log(L ** 2 / (S * X)) / (v * sqrt(T)) + mu * v * sqrt(T)
    e2 = e1 - (mu + 0.5) * v * sqrt(T)
    e3 = log(L ** 2 / (S * X)) / (v * sqrt(T)) + (mu - 1) * v * sqrt(T)
    e4 = e3 - (mu - 0.5) * v * sqrt(T)
    
    value = (
        eta * 1 / (U - L) * (
            S * exp((b - r) * T)
            * S ** (-2 * mu) * (S * X) ** (mu + 0.5)
            / (2 * (mu + 0.5)) * (
                (U ** 2 / (S * X)) ** (mu + 0.5) * cnd(eta * d1)
                - lambda1 * cnd(eta * d2)
                - (L ** 2 / (S * X)) ** (mu + 0.5) * cnd(eta * e1)
                + lambda1 * cnd(eta * e2)
            )
            - X * exp(-r * T) * S ** (-2 * (mu - 1))
            * (S * X) ** (mu - 0.5) / (2 * (mu - 0.5))
            * (
                (U ** 2 / (S * X)) ** (mu - 0.5) * cnd(eta * d3)
                - lambda2 * cnd(eta * d4)
                - (L ** 2 / (S * X)) ** (mu - 0.5) * cnd(eta * e3)
                + lambda2 * cnd(eta * e4)
            )
        )
    )
    
    if is_in:
        return value
    else:
        return bs_price(is_call, S, X, T, r, b, v, cdf=cnd) - value
