"""Look backs"""

from math import exp, log, sqrt
from typing import Callable

from ...distributions import CDF, CBND


def price(
        is_call: bool,
        S: float,
        K: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        cdf: Callable[[float], float] = CDF,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:
    # Partial-time fixed strike look-back options

    d1 = (log(S / K) + (b + v ** 2 / 2) * T2) / (v * sqrt(T2))
    d2 = d1 - v * sqrt(T2)
    e1 = ((b + v ** 2 / 2) * (T2 - t1)) / (v * sqrt(T2 - t1))
    e2 = e1 - v * sqrt(T2 - t1)
    f1 = (log(S / K) + (b + v ** 2 / 2) * t1) / (v * sqrt(t1))
    f2 = f1 - v * sqrt(t1)
    if is_call:
        return (
            S * exp((b - r) * T2) * cdf(d1)
            - exp(-r * T2) * K * cdf(d2)
            + S * exp(-r * T2) * v ** 2 / (2 * b) * (
                -(S / K) ** (-2 * b / v ** 2) * cbnd(
                    d1 - 2 * b * sqrt(T2) / v,
                    -f1 + 2 * b * sqrt(t1) / v,
                    -sqrt(t1 / T2)
                )
                + exp(b * T2) * cbnd(e1, d1, sqrt(1 - t1 / T2))
            )
            - S * exp((b - r) * T2) * cbnd(-e1, d1, -sqrt(1 - t1 / T2))
            - K * exp(-r * T2) * cbnd(f2, -d2, -sqrt(t1 / T2))
            + exp(-b * (T2 - t1)) * (1 - v ** 2 / (2 * b)) * S
            * exp((b - r) * T2) * cdf(f1) * cdf(-e2))
    else:
        return (
            K * exp(-r * T2) * cdf(-d2)
            - S * exp((b - r) * T2) * cdf(-d1)
            + S * exp(-r * T2) * v ** 2 / (2 * b) * (
                (S / K) ** (-2 * b / v ** 2) * cbnd(
                    -d1 + 2 * b * sqrt(T2) / v,
                    f1 - 2 * b * sqrt(t1) / v,
                    -sqrt(t1 / T2)
                )
                - exp(b * T2) * cbnd(-e1, -d1, sqrt(1 - t1 / T2))
            )
            + S * exp((b - r) * T2) * cbnd(e1, -d1, -sqrt(1 - t1 / T2))
            + K * exp(-r * T2) * cbnd(-f2, d2, -sqrt(t1 / T2))
            - exp(-b * (T2 - t1)) * (1 - v ** 2 / (2 * b)) * S
            * exp((b - r) * T2) * cdf(-f1) * cdf(e2))
