"""American"""

from math import exp, log, sqrt
from typing import Callable, Literal, Optional

from ...distributions import CND, CBND, ND
from ...european.black_scholes.analytic import price as bs_price

def phi(
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
    return exp(lambda_) * S ** gamma_ * (cnd(d) - (i / S) ** kappa * cnd(d - 2 * log(i / S) / (v * sqrt(T))))

def ksi(
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
        *,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    e1 = (log(S / I1) + (b + (gamma_ - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e2 = (log(I2 ** 2 / (S * I1)) + (b + (gamma_ - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e3 = (log(S / I1) - (b + (gamma_ - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e4 = (log(I2 ** 2 / (S * I1)) - (b + (gamma_ - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    
    f1 = (log(S / h) + (b + (gamma_ - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f2 = (log(I2 ** 2 / (S * h)) + (b + (gamma_ - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f3 = (log(I1 ** 2 / (S * h)) + (b + (gamma_ - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f4 = (log(S * I1 ** 2 / (h * I2 ** 2)) + (b + (gamma_ - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    
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

    

def BSAmericanCallApprox(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
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
                - alpha * phi(S, T, beta, i, i, r, b, v, cnd=cnd)
                + phi(S, T, 1, i, i, r, b, v, cnd=cnd)
                - phi(S, T, 1, X, i, r, b, v, cnd=cnd)
                - X * phi(S, T, 0, i, i, r, b, v, cnd=cnd)
                + X * phi(S, T, 0, X, i, r, b, v, cnd=cnd)
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
        return BSAmericanCallApprox(S, X, T, r, b, v)
    else:
        # Use the Bjerksund and Stensland put-call transformation
        return BSAmericanCallApprox(X, S, T, r - b, -b, v)
    