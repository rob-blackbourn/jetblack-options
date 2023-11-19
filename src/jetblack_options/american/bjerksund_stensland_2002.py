"""The Bjerksund and Stensland (2002)"""

from math import exp, log, sqrt
from statistics import NormalDist
from typing import Callable

from ..distributions import CBND as cbnd
from ..european.generalised_black_scholes import price as bs_price

norm = NormalDist()
cdf = norm.cdf


def _phi(
        S: float,
        T: float,
        gamma_: float,
        h: float,
        i: float,
        r: float,
        b: float,
        v: float,
) -> float:
    lambda_ = (-r + gamma_ * b + 0.5 * gamma_ * (gamma_ - 1) * v ** 2) * T
    d = -(log(S / h) + (b + (gamma_ - 0.5) * v ** 2) * T) / (v * sqrt(T))
    kappa = 2 * b / v ** 2 + 2 * gamma_ - 1
    return (
        exp(lambda_) * S ** gamma_ * (
            cdf(d)
            - (i / S) ** kappa * cdf(d - 2 * log(i / S) / (v * sqrt(T))))
    )


def _ksi(
        S: float,
        T2: float,
        gamma_: float,
        h: float,
        I2: float,
        I1: float,
        t1: float,
        r: float,
        b: float,
        v: float,
) -> float:
    e1 = (log(S / I1) + (b + (gamma_ - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e2 = (log(I2 ** 2 / (S * I1)) +
          (b + (gamma_ - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e3 = (log(S / I1) - (b + (gamma_ - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e4 = (log(I2 ** 2 / (S * I1)) -
          (b + (gamma_ - 0.5) * v ** 2) * t1) / (v * sqrt(t1))

    f1 = (log(S / h) + (b + (gamma_ - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f2 = (log(I2 ** 2 / (S * h)) +
          (b + (gamma_ - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f3 = (log(I1 ** 2 / (S * h)) +
          (b + (gamma_ - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f4 = (log(S * I1 ** 2 / (h * I2 ** 2)) +
          (b + (gamma_ - 0.5) * v ** 2) * T2) / (v * sqrt(T2))

    rho = sqrt(t1 / T2)
    lambda_ = -r + gamma_ * b + 0.5 * gamma_ * (gamma_ - 1) * v ** 2
    kappa = 2 * b / (v ** 2) + (2 * gamma_ - 1)

    return (
        exp(lambda_ * T2) *
        S ** gamma_ *
        (
            cbnd(-e1, -f1, rho) -
            (I2 / S) ** kappa * cbnd(-e2, -f2, rho) -
            (I1 / S) ** kappa * cbnd(-e3, -f3, -rho) +
            (I1 / I2) ** kappa * cbnd(-e4, -f4, -rho)
        )
    )


def _call_price(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
) -> float:

    t1 = 1 / 2 * (sqrt(5) - 1) * T

    if b >= r:
        # Use Black-Scholes as it is never optimal to exercise before maturity.
        return bs_price(True, S, K, T, r, b, v)

    beta = (
        (1 / 2 - b / v ** 2)
        + sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
    )
    b_infinity = beta / (beta - 1) * K
    b0 = max(K, r / (r - b) * K)

    ht1 = -(b * t1 + 2 * v * sqrt(t1)) * K ** 2 / ((b_infinity - b0) * b0)
    ht2 = -(b * T + 2 * v * sqrt(T)) * K ** 2 / ((b_infinity - b0) * b0)
    I1 = b0 + (b_infinity - b0) * (1 - exp(ht1))
    I2 = b0 + (b_infinity - b0) * (1 - exp(ht2))
    alfa1 = (I1 - K) * I1 ** (-beta)
    alfa2 = (I2 - K) * I2 ** (-beta)

    if S >= I2:
        return S - K
    else:
        return (
            alfa2 * S ** beta
            - alfa2 * _phi(S, t1, beta, I2, I2, r, b, v)
            + _phi(S, t1, 1, I2, I2, r, b, v)
            - _phi(S, t1, 1, I1, I2, r, b, v)
            - K * _phi(S, t1, 0, I2, I2, r, b, v)
            + K * _phi(S, t1, 0, I1, I2, r, b, v)
            + alfa1 * _phi(S, t1, beta, I1, I2, r, b, v)
            - alfa1 * _ksi(S, T, beta, I1, I2, I1, t1, r, b, v)
            + _ksi(S, T, 1, I1, I2, I1, t1, r, b, v)
            - _ksi(S, T, 1, K, I2, I1, t1, r, b, v)
            - K * _ksi(S, T, 0, I1, I2, I1, t1, r, b, v)
            + K * _ksi(S, T, 0, K, I2, I1, t1, r, b, v)
        )


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
) -> float:
    """The Bjerksund and Stensland (2002) American approximation.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The current asset price.
        K (float): The option strike price
        T (float): The time to maturity of the option in years.
        r (float): The risk free rate.
        b (float): The cost of carry of the asset.
        v (float): The volatility of the asset.

    Returns:
        float: The price of the option.
    """

    if is_call:
        return _call_price(S, K, T, r, b, v)
    else:
        # Use the Bjerksund and Stensland put-call transformation
        return _call_price(K, S, T, r - b, -b, v)
