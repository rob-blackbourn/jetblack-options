"""The Bjerksund and Stensland (1993)"""

from math import exp, log, sqrt
from statistics import NormalDist

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
            cdf(d) - (i / S) ** kappa * cdf(d - 2 * log(i / S) / (v * sqrt(T)))
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

    if b >= r:
        # We can use Black-Scholes as it is never optimal to exercise before
        # maturity.
        return bs_price(True, S, K, T, r, b, v)

    beta = (
        (1 / 2 - b / v ** 2)
        + sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
    )
    b_infinity = beta / (beta - 1) * K
    b0 = max(K, r / (r - b) * K)
    ht = -(b * T + 2 * v * sqrt(T)) * b0 / (b_infinity - b0)
    i = b0 + (b_infinity - b0) * (1 - exp(ht))
    alpha = (i - K) * i ** (-beta)
    if S >= i:
        return S - K
    else:
        return (
            alpha * S ** beta
            - alpha * _phi(S, T, beta, i, i, r, b, v)
            + _phi(S, T, 1, i, i, r, b, v)
            - _phi(S, T, 1, K, i, r, b, v)
            - K * _phi(S, T, 0, i, i, r, b, v)
            + K * _phi(S, T, 0, K, i, r, b, v)
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
    """The Bjerksund and Stensland (1993) American approximation.

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
