"""American"""

from math import exp, log, sqrt
from typing import Callable, Literal, Optional

from ...distributions import CND
from ...european.black_scholes.analytic import price as bs_price

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
        cnd: Callable[[float], float] = CND,
) -> float:
    lambda_ = (-r + gamma_ * b + 0.5 * gamma_ * (gamma_ - 1) * v ** 2) * T
    d = -(log(S / h) + (b + (gamma_ - 0.5) * v ** 2) * T) / (v * sqrt(T))
    kappa = 2 * b / v ** 2 + 2 * gamma_ - 1
    return (
        exp(lambda_) * S ** gamma_ * (
            cnd(d) - (i / S) ** kappa * cnd(d - 2 * log(i / S) / (v * sqrt(T)))
        )
    )
    

def _call_price(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    
    if b >= r: # Never optimal to exercise before maturity
        return bs_price(True, S, X, T, r, b, v, cnd=cnd)
    else:
        beta = (
            (1 / 2 - b / v ** 2)
            + sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
        )
        b_infinity = beta / (beta - 1) * X
        B0 = max(X, r / (r - b) * X)
        ht = -(b * T + 2 * v * sqrt(T)) * B0 / (b_infinity - B0)
        i = B0 + (b_infinity - B0) * (1 - exp(ht))
        alpha = (i - X) * i ** (-beta)
        if S >= i:
            return S - X
        else:
            return (
                alpha * S ** beta
                - alpha * _phi(S, T, beta, i, i, r, b, v, cnd=cnd)
                + _phi(S, T, 1, i, i, r, b, v, cnd=cnd)
                - _phi(S, T, 1, X, i, r, b, v, cnd=cnd)
                - X * _phi(S, T, 0, i, i, r, b, v, cnd=cnd)
                + X * _phi(S, T, 0, X, i, r, b, v, cnd=cnd)
            )

# The Bjerksund and Stensland (1993) American approximation
def price(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:
    if is_call:
        return _call_price(S, X, T, r, b, v, cnd=cnd)
    else:
        # Use the Bjerksund and Stensland put-call transformation
        return _call_price(X, S, T, r - b, -b, v, cnd=cnd)
    