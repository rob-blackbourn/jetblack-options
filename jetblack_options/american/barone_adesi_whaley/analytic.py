"""American"""

from math import exp, log, sqrt
from typing import Callable

from ...distributions import CND, ND
from ...european.black_scholes.analytic import price as bs_price


def _kc(
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    """Newton Raphson algorithm to solve for the critical commodity price for a
    Call.

    Args:
        X (float): The strike.
        T (float): The time to expiry in years.
        r (float): The risk free rate.
        b (float): The asset growth.
        v (float): The volatility.
        nd (Callable[[float], float], optional): A function returning the normal
            distribution. Defaults to ND.
        cnd (Callable[[float], float], optional): A function returning the
            cumulative normal distribution. Defaults to CND.

    Returns:
        float: The price.
    """

    # Calculate the seed value Si
    N = 2 * b / v ** 2
    m = 2 * r / v ** 2
    q2u = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * m)) / 2
    su = X / (1 - 1 / q2u)
    h2 = -(b * T + 2 * v * sqrt(T)) * X / (su - X)
    Si = X + (su - X) * (1 - exp(h2))

    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Si / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q2 = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * k)) / 2
    LHS = Si - X
    RHS = (
        bs_price(True, Si, X, T, r, b, v, cdf=cnd) +
        (1 - exp((b - r) * T) * cnd(d1)) * Si / Q2
    )
    bi = (
        exp((b - r) * T) * cnd(d1) * (1 - 1 / Q2) +
        (1 - exp((b - r) * T) * cnd(d1) / (v * sqrt(T))) / Q2
    )
    E = 0.000001
    # Newton Raphson algorithm for finding critical price Si
    while abs(LHS - RHS) / X > E:
        Si = (X + RHS - bi * Si) / (1 - bi)
        d1 = (log(Si / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        LHS = Si - X
        RHS = (
            bs_price(True, Si, X, T, r, b, v, cdf=cnd) +
            (1 - exp((b - r) * T) * cnd(d1)) * Si / Q2
        )
        bi = (
            exp((b - r) * T) * cnd(d1) * (1 - 1 / Q2) +
            (1 - exp((b - r) * T) * nd(d1) / (v * sqrt(T))) / Q2
        )

    return Si


def _call_price(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:

    if b >= r:
        return bs_price(True, S, X, T, r, b, v, cdf=cnd)
    else:
        Sk = _kc(X, T, r, b, v, nd=nd, cnd=cnd)
        N = 2 * b / v ** 2
        k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
        d1 = (log(Sk / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        Q2 = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * k)) / 2
        a2 = (Sk / Q2) * (1 - exp((b - r) * T) * cnd(d1))
        if S < Sk:
            return (
                bs_price(True, S, X, T, r, b, v, cdf=cnd)
                + a2 * (S / Sk) ** Q2
            )
        else:
            return S - X


def _kp(
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    # Newton Raphson algorithm to solve for the critical commodity price for a Put

    # Calculation of seed value, Si
    N = 2 * b / v ** 2
    m = 2 * r / v ** 2
    q1u = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * m)) / 2
    su = X / (1 - 1 / q1u)
    h1 = (b * T - 2 * v * sqrt(T)) * X / (X - su)
    Si = su + (X - su) * exp(h1)

    k = 2 * r / (v * 2 * (1 - exp(-r * T)))
    d1 = (log(Si / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q1 = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * k)) / 2
    LHS = X - Si
    RHS = (
        bs_price(False, Si, X, T, r, b, v, cdf=cnd)
        - (1 - exp((b - r) * T) * cnd(-d1)) * Si / Q1
    )
    bi = (
        -exp((b - r) * T) * cnd(-d1) * (1 - 1 / Q1)
        - (1 + exp((b - r) * T) * nd(-d1) / (v * sqrt(T))) / Q1
    )
    E = 0.000001
    # Newton Raphson algorithm for finding critical price Si
    while abs(LHS - RHS) / X > E:
        Si = (X - RHS + bi * Si) / (1 + bi)
        d1 = (log(Si / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        LHS = X - Si
        RHS = (
            bs_price(False, Si, X, T, r, b, v, cdf=cnd)
            - (1 - exp((b - r) * T) * cnd(-d1)) * Si / Q1
        )
        bi = (
            -exp((b - r) * T) * cnd(-d1) * (1 - 1 / Q1)
            - (1 + exp((b - r) * T) * cnd(-d1) / (v * sqrt(T))) / Q1
        )

    return Si


def _put_price(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:

    Sk = _kp(X, T, r, b, v, nd=nd, cnd=cnd)
    N = 2 * b / v ** 2
    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Sk / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q1 = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * k)) / 2
    a1 = -(Sk / Q1) * (1 - exp((b - r) * T) * cnd(-d1))

    if S > Sk:
        return bs_price(False, S, X, T, r, b, v, cdf=cnd) + a1 * (S / Sk) ** Q1
    else:
        return X - S


def price(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    # The Barone-Adesi and Whaley (1987) American approximation
    if is_call:
        return _call_price(S, X, T, r, b, v, nd=nd, cnd=cnd)
    else:
        return _put_price(S, X, T, r, b, v, nd=nd, cnd=cnd)
