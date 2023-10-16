"""American"""

from math import exp, log, sqrt
from typing import Callable, Literal, Optional

from ..distributions import CDF
from ..european.black_scholes_merton import price as bs_price

def _phi(
        S: float,
        T: float,
        gamma_: float,
        h: float,
        i: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF,
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
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if b >= r: # Never optimal to exercise before maturity
        return bs_price(True, S, K, T, r, b, v, cdf=cdf)
    else:
        beta = (
            (1 / 2 - b / v ** 2)
            + sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
        )
        b_infinity = beta / (beta - 1) * K
        B0 = max(K, r / (r - b) * K)
        ht = -(b * T + 2 * v * sqrt(T)) * B0 / (b_infinity - B0)
        i = B0 + (b_infinity - B0) * (1 - exp(ht))
        alpha = (i - K) * i ** (-beta)
        if S >= i:
            return S - K
        else:
            return (
                alpha * S ** beta
                - alpha * _phi(S, T, beta, i, i, r, b, v, cdf=cdf)
                + _phi(S, T, 1, i, i, r, b, v, cdf=cdf)
                - _phi(S, T, 1, K, i, r, b, v, cdf=cdf)
                - K * _phi(S, T, 0, i, i, r, b, v, cdf=cdf)
                + K * _phi(S, T, 0, K, i, r, b, v, cdf=cdf)
            )

# The Bjerksund and Stensland (1993) American approximation
def price(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CDF
) -> float:
    if is_call:
        return _call_price(S, K, T, r, b, v, cdf=cdf)
    else:
        # Use the Bjerksund and Stensland put-call transformation
        return _call_price(K, S, T, r - b, -b, v, cdf=cdf)
    