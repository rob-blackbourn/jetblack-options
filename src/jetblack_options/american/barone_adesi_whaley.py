"""Barone-Adesi and Whaley"""

from math import exp, log, sqrt
from statistics import NormalDist

from ..european.generalised_black_scholes import price as bs_price

norm = NormalDist()
cdf = norm.cdf
pdf = norm.pdf
inv_cdf = norm.inv_cdf


def _kc(
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
) -> float:
    """Newton Raphson algorithm to solve for the critical commodity price for a
    Call.

    Args:
        K (float): The strike.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        b (float): The asset growth.
        v (float): The volatility.

    Returns:
        float: The price.
    """

    # Calculate the seed value Si
    n = 2 * b / v ** 2
    m = 2 * r / v ** 2
    q2u = (-(n - 1) + sqrt((n - 1) ** 2 + 4 * m)) / 2
    su = K / (1 - 1 / q2u)
    h2 = -(b * T + 2 * v * sqrt(T)) * K / (su - K)
    Si = K + (su - K) * (1 - exp(h2))

    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Si / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    q2 = (-(n - 1) + sqrt((n - 1) ** 2 + 4 * k)) / 2
    lhs = Si - K
    rhs = (
        bs_price(True, Si, K, T, r, b, v) +
        (1 - exp((b - r) * T) * cdf(d1)) * Si / q2
    )
    bi = (
        exp((b - r) * T) * cdf(d1) * (1 - 1 / q2) +
        (1 - exp((b - r) * T) * cdf(d1) / (v * sqrt(T))) / q2
    )
    epsilon = 0.000001
    # Using the Newton Raphson algorithm solve for Si
    while abs(lhs - rhs) / K > epsilon:
        Si = (K + rhs - bi * Si) / (1 - bi)
        d1 = (log(Si / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        lhs = Si - K
        rhs = (
            bs_price(True, Si, K, T, r, b, v) +
            (1 - exp((b - r) * T) * cdf(d1)) * Si / q2
        )
        bi = (
            exp((b - r) * T) * cdf(d1) * (1 - 1 / q2) +
            (1 - exp((b - r) * T) * pdf(d1) / (v * sqrt(T))) / q2
        )

    return Si


def _call_price(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
) -> float:

    if b >= r:
        return bs_price(True, S, K, T, r, b, v)

    Sk = _kc(K, T, r, b, v)
    n = 2 * b / v ** 2
    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Sk / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    q2 = (-(n - 1) + sqrt((n - 1) ** 2 + 4 * k)) / 2
    a2 = (Sk / q2) * (1 - exp((b - r) * T) * cdf(d1))
    if S < Sk:
        return (
            bs_price(True, S, K, T, r, b, v)
            + a2 * (S / Sk) ** q2
        )
    else:
        return S - K


def _kp(
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
) -> float:
    # Newton Raphson algorithm to solve for the critical commodity price for a Put

    # Calculation of seed value, Si
    n = 2 * b / v ** 2
    m = 2 * r / v ** 2
    q1u = (-(n - 1) - sqrt((n - 1) ** 2 + 4 * m)) / 2
    su = K / (1 - 1 / q1u)
    h1 = (b * T - 2 * v * sqrt(T)) * K / (K - su)
    Si = su + (K - su) * exp(h1)

    k = 2 * r / (v * 2 * (1 - exp(-r * T)))
    d1 = (log(Si / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    q1 = (-(n - 1) - sqrt((n - 1) ** 2 + 4 * k)) / 2
    lhs = K - Si
    rhs = (
        bs_price(False, Si, K, T, r, b, v)
        - (1 - exp((b - r) * T) * cdf(-d1)) * Si / q1
    )
    bi = (
        -exp((b - r) * T) * cdf(-d1) * (1 - 1 / q1)
        - (1 + exp((b - r) * T) * pdf(-d1) / (v * sqrt(T))) / q1
    )
    epsilon = 0.000001
    # Using the Newton Raphson algorithm, solve for Si.
    while abs(lhs - rhs) / K > epsilon:
        Si = (K - rhs + bi * Si) / (1 + bi)
        d1 = (log(Si / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        lhs = K - Si
        rhs = (
            bs_price(False, Si, K, T, r, b, v)
            - (1 - exp((b - r) * T) * cdf(-d1)) * Si / q1
        )
        bi = (
            -exp((b - r) * T) * cdf(-d1) * (1 - 1 / q1)
            - (1 + exp((b - r) * T) * cdf(-d1) / (v * sqrt(T))) / q1
        )

    return Si


def _put_price(
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
) -> float:

    Sk = _kp(K, T, r, b, v)
    n = 2 * b / v ** 2
    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Sk / K) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    q1 = (-(n - 1) - sqrt((n - 1) ** 2 + 4 * k)) / 2
    a1 = -(Sk / q1) * (1 - exp((b - r) * T) * cdf(-d1))

    if S > Sk:
        return bs_price(False, S, K, T, r, b, v) + a1 * (S / Sk) ** q1
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

    Returns:
        float: The price of the option.
    """

    if is_call:
        return _call_price(S, K, T, r, b, v)
    else:
        return _put_price(S, K, T, r, b, v)
