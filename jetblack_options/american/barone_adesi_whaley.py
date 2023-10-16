"""Barone-Adesi and Whaley"""

from math import exp, log, sqrt
from typing import Callable

from ..distributions import CDF, PDF
from ..european.black_scholes_merton import price as bs_price


def _kc(
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF,
        cdf: Callable[[float], float] = CDF
) -> float:
    """Newton Raphson algorithm to solve for the critical commodity price for a
    Call.

    Args:
        K (float): The strike.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        b (float): The asset growth.
        v (float): The volatility.
        pdf (Callable[[float], float], optional): A function returning the normal
            distribution. Defaults to PDF.
        cdf (Callable[[float], float], optional): A function returning the
            cumulative normal distribution. Defaults to CDF.

    Returns:
        float: The price.
    """

    # Calculate the seed value Si
    N = 2 * b / v ** 2
    m = 2 * r / v ** 2
    q2u = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * m)) / 2
    su = K / (1 - 1 / q2u)
    h2 = -(b * T + 2 * v * sqrt(T)) * K / (su - K)
    Si = K + (su - K) * (1 - exp(h2))

    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Si / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q2 = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * k)) / 2
    LHS = Si - K
    RHS = (
        bs_price(True, Si, K, T, r, b, v, cdf=cdf) +
        (1 - exp((b - r) * T) * cdf(d1)) * Si / Q2
    )
    bi = (
        exp((b - r) * T) * cdf(d1) * (1 - 1 / Q2) +
        (1 - exp((b - r) * T) * cdf(d1) / (v * sqrt(T))) / Q2
    )
    E = 0.000001
    # Newton Raphson algorithm for finding critical price Si
    while abs(LHS - RHS) / K > E:
        Si = (K + RHS - bi * Si) / (1 - bi)
        d1 = (log(Si / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        LHS = Si - K
        RHS = (
            bs_price(True, Si, K, T, r, b, v, cdf=cdf) +
            (1 - exp((b - r) * T) * cdf(d1)) * Si / Q2
        )
        bi = (
            exp((b - r) * T) * cdf(d1) * (1 - 1 / Q2) +
            (1 - exp((b - r) * T) * pdf(d1) / (v * sqrt(T))) / Q2
        )

    return Si


def _call_price(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF,
        cdf: Callable[[float], float] = CDF
) -> float:

    if b >= r:
        return bs_price(True, S, K, T, r, b, v, cdf=cdf)
    else:
        Sk = _kc(K, T, r, b, v, pdf=pdf, cdf=cdf)
        N = 2 * b / v ** 2
        k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
        d1 = (log(Sk / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        Q2 = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * k)) / 2
        a2 = (Sk / Q2) * (1 - exp((b - r) * T) * cdf(d1))
        if S < Sk:
            return (
                bs_price(True, S, K, T, r, b, v, cdf=cdf)
                + a2 * (S / Sk) ** Q2
            )
        else:
            return S - K


def _kp(
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF,
        cdf: Callable[[float], float] = CDF
) -> float:
    # Newton Raphson algorithm to solve for the critical commodity price for a Put

    # Calculation of seed value, Si
    N = 2 * b / v ** 2
    m = 2 * r / v ** 2
    q1u = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * m)) / 2
    su = K / (1 - 1 / q1u)
    h1 = (b * T - 2 * v * sqrt(T)) * K / (K - su)
    Si = su + (K - su) * exp(h1)

    k = 2 * r / (v * 2 * (1 - exp(-r * T)))
    d1 = (log(Si / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q1 = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * k)) / 2
    LHS = K - Si
    RHS = (
        bs_price(False, Si, K, T, r, b, v, cdf=cdf)
        - (1 - exp((b - r) * T) * cdf(-d1)) * Si / Q1
    )
    bi = (
        -exp((b - r) * T) * cdf(-d1) * (1 - 1 / Q1)
        - (1 + exp((b - r) * T) * pdf(-d1) / (v * sqrt(T))) / Q1
    )
    E = 0.000001
    # Newton Raphson algorithm for finding critical price Si
    while abs(LHS - RHS) / K > E:
        Si = (K - RHS + bi * Si) / (1 + bi)
        d1 = (log(Si / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        LHS = K - Si
        RHS = (
            bs_price(False, Si, K, T, r, b, v, cdf=cdf)
            - (1 - exp((b - r) * T) * cdf(-d1)) * Si / Q1
        )
        bi = (
            -exp((b - r) * T) * cdf(-d1) * (1 - 1 / Q1)
            - (1 + exp((b - r) * T) * cdf(-d1) / (v * sqrt(T))) / Q1
        )

    return Si


def _put_price(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF,
        cdf: Callable[[float], float] = CDF
) -> float:

    Sk = _kp(K, T, r, b, v, pdf=pdf, cdf=cdf)
    N = 2 * b / v ** 2
    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Sk / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q1 = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * k)) / 2
    a1 = -(Sk / Q1) * (1 - exp((b - r) * T) * cdf(-d1))

    if S > Sk:
        return bs_price(False, S, K, T, r, b, v, cdf=cdf) + a1 * (S / Sk) ** Q1
    else:
        return K - S


def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        pdf: Callable[[float], float] = PDF,
        cdf: Callable[[float], float] = CDF
) -> float:
    """The Barone-Adesi and Whaley (1987) American approximation.

    Args:
        is_call (bool): True for a call, false for a put.
        S (float): The asset price.
        K (float): The strike price.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        b (float): The cost of carry.
        v (float): The asset volatility.
        cdf (Callable[[float], float], optional): The cumulative density function. Defaults to CDF.
        pdf (Callable[[float], float], optional): The probability density function. Defaults to PDF.

    Returns:
        float: The price of the option.
    """
    # The Barone-Adesi and Whaley (1987) American approximation
    if is_call:
        return _call_price(S, K, T, r, b, v, pdf=pdf, cdf=cdf)
    else:
        return _put_price(S, K, T, r, b, v, pdf=pdf, cdf=cdf)
