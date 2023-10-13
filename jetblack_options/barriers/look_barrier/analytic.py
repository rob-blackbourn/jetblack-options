"""Barriers"""

from math import exp, log, sqrt
from typing import Callable

from ...distributions import CND, CBND
from .partial_fixed_lookback import price as lookback_price

def price(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:
    # Look-barrier options

    hh = log(H / S)
    k = log(X / S)
    mu1 = b - v ** 2 / 2
    mu2 = b + v ** 2 / 2
    rho = sqrt(t1 / T2)
    
    if is_call:
        eta = 1
        m = min(hh, k)
    else:
        eta = -1
        m = max(hh, k)
    
    g1 = (
        cnd(eta * (hh - mu2 * t1) / (v * sqrt(t1))) -
        exp(2 * mu2 * hh / v ** 2) * cnd(eta * (-hh - mu2 * t1) / (v * sqrt(t1)))
    ) - (
        cnd(eta * (m - mu2 * t1) / (v * sqrt(t1))) -
        exp(2 * mu2 * hh / v ** 2) * cnd(eta * (m - 2 * hh - mu2 * t1) / (v * sqrt(t1)))
    )
    g2 = (
        cnd(eta * (hh - mu1 * t1) / (v * sqrt(t1))) -
        exp(2 * mu1 * hh / v ** 2) * cnd(eta * (-hh - mu1 * t1) / (v * sqrt(t1)))
    ) - (
        cnd(eta * (m - mu1 * t1) / (v * sqrt(t1))) -
        exp(2 * mu1 * hh / v ** 2) * cnd(eta * (m - 2 * hh - mu1 * t1) / (v * sqrt(t1)))
    )

    part1 = (
        S *
        exp((b - r) * T2) *
        (1 + v ** 2 / (2 * b)) *
        (
            cbnd(
                eta * (m - mu2 * t1) / (v * sqrt(t1)),
                eta * (-k + mu2 * T2) / (v * sqrt(T2)),
                -rho
            ) -
            exp(2 * mu2 * hh / v ** 2) *
            cbnd(
                eta * (m - 2 * hh - mu2 * t1) / (v * sqrt(t1)),
                eta * (2 * hh - k + mu2 * T2) / (v * sqrt(T2)),
                -rho
            )
        )
    )
    part2 = (
        -exp(-r * T2) * X * (
            cbnd(
                eta * (m - mu1 * t1) / (v * sqrt(t1)),
                eta * (-k + mu1 * T2) / (v * sqrt(T2)),
                -rho
            ) -
            exp(2 * mu1 * hh / v ** 2) *
            cbnd(
                eta * (m - 2 * hh - mu1 * t1) / (v * sqrt(t1)),
                eta * (2 * hh - k + mu1 * T2) / (v * sqrt(T2)),
                -rho
            )
        )
    )
    part3 = (
        -exp(-r * T2) * v ** 2 / (2 * b) * (
            S * (S / X) ** (-2 * b / v ** 2) *
            cbnd(
                eta * (m + mu1 * t1) / (v * sqrt(t1)),
                eta * (-k - mu1 * T2) / (v * sqrt(T2)),
                -rho
            ) -
            H * (H / X) ** (-2 * b / v ** 2) *
            cbnd(
                eta * (m - 2 * hh + mu1 * t1) / (v * sqrt(t1)),
                eta * (2 * hh - k - mu1 * T2) / (v * sqrt(T2)),
                -rho
            )
        )
    )
    part4 = (
        S * exp((b - r) * T2) * (
            (1 + v ** 2 / (2 * b)) *
            cnd(eta * mu2 * (T2 - t1) / (v * sqrt(T2 - t1))) +
            exp(-b * (T2 - t1)) * (1 - v ** 2 / (2 * b)) *
            cnd(eta * (-mu1 * (T2 - t1)) / (v * sqrt(T2 - t1)))
        ) * g1 - exp(-r * T2) * X * g2)
    out_value = eta * (part1 + part2 + part3 + part4)

    if not is_in:
        return out_value
    else:
        return lookback_price(is_call, S, X, t1, T2, r, b, v) - out_value
    