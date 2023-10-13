"""Turnbull-Wakeman analytic"""

from math import exp, log, sqrt
from typing import Literal, Optional, Callable

from ...distributions import CND

from .analytic import price

def delta(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, X, T, T2, r, b, v, cnd=cnd) -
        price(is_call, S - dS, SA, X, T, T2, r, b, v, cnd=cnd)
    ) / (2 * dS)

def elasticity(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, X, T, T2, r, b, v, cnd=cnd)
        - price(is_call, S - dS, SA, X, T, T2, r, b, v, cnd=cnd)
    ) / (2 * dS) * S / price(is_call, S, SA, X, T, T2, r, b, v, cnd=cnd)

def gamma(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, X, T, T2, r, b, v, cnd=cnd)
        - 2 * price(is_call, S, SA, X, T, T2, r, b, v, cnd=cnd)
        + price(is_call, S - dS, SA, X, T, T2, r, b, v, cnd=cnd)
    ) / dS ** 2

def dgamma_dvol(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, X, T, T2, r, b, v + dv, cnd=cnd)
        - 2 * price(is_call, S, SA, X, T, T2, r, b, v + dv, cnd=cnd)
        + price(is_call, S - dS, SA, X, T, T2, r, b, v + dv, cnd=cnd)
        - price(is_call, S + dS, SA, X, T, T2, r, b, v - dv, cnd=cnd)
        + 2 * price(is_call, S, SA, X, T, T2, r, b, v - dv, cnd=cnd)
        - price(is_call, S - dS, SA, X, T, T2, r, b, v - dv, cnd=cnd)
    ) / (2 * dv * dS ** 2) / 100

def gammap(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return S / 100 * (
        price(is_call, S + dS, SA, X, T, T2, r, b, v, cnd=cnd)
        - 2 * price(is_call, S, SA, X, T, T2, r, b, v, cnd=cnd)
        + price(is_call, S - dS, SA, X, T, T2, r, b, v, cnd=cnd)
    ) / dS ** 2

def ddelta_dvol(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, X, T, T2, r, b, v + dv, cnd=cnd)
        -  price(is_call, S + dS, SA, X, T, T2, r, b, v - dv, cnd=cnd)
        - price(is_call, S - dS, SA, X, T, T2, r, b, v + dv, cnd=cnd)
        + price(is_call, S - dS, SA, X, T, T2, r, b, v - dv, cnd=cnd)
    ) / (4 * dS * 0.01) / 100

def vega(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X, T, T2, r, b, v + dv, cnd=cnd)
        - price(is_call, S, SA, X, T, T2, r, b, v - dv, cnd=cnd)
    ) / 2

def vomma(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X, T, T2, r, b, v + dv, cnd=cnd)
        - 2 * price(is_call, S, SA, X, T, T2, r, b, v, cnd=cnd)
        + price(is_call, S, SA, X, T, T2, r, b, v - dv, cnd=cnd)
    ) / dv ** 2 / 10000

def vegap(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X, T, T2, r, b, v + dv)
        - price(is_call, S, SA, X, T, T2, r, b, v - dv)
    ) * v / 0.1 / 2

def dvega_dvol(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X, T, T2, r, b, v + dv, cnd=cnd)
        - 2 * price(is_call, S, SA, X, T, T2, r, b, v, cnd=cnd)
        + price(is_call, S, SA, X, T, T2, r, b, v - dv, cnd=cnd)
    )

def theta(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365,
        cnd: Callable[[float], float] = CND
) -> float:
    if T <= dT:
        return (
            price(is_call, S, SA, X, 0.00001, T2, r, b, v)
            - price(is_call, S, SA, X, T, T2, r, b, v)
        )
    else:
        return (
            price(is_call, S, SA, X, T - dT, T2, r, b, v)
            - price(is_call, S, SA, X, T, T2, r, b, v)
        )

def rho(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X, T, T2, r + dr, b + dr, v)
        - price(is_call, S, SA, X, T, T2, r - dr, b - dr, v)
    ) / 2

def futures_rho(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X, T, T2, r + dr, b, v, cnd=cnd)
        - price(is_call, S, SA, X, T, T2, r - dr, b, v, cnd=cnd)
    ) / 2

def rho2(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X, T, T2, r, b - db, v, cnd=cnd)
        - price(is_call, S, SA, X, T, T2, r, b + db, v, cnd=cnd)
    ) / 2

def carry(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X, T, T2, r, b + db, v, cnd=cnd)
        - price(is_call, S, SA, X, T, T2, r, b - db, v, cnd=cnd)
    ) / 2

def speed(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return 1 / dS ** 3 * (
        price(is_call, S + 2 * dS, SA, X, T, T2, r, b, v, cnd=cnd)
        - 3 * price(is_call, S + dS, SA, X, T, T2, r, b, v, cnd=cnd)
        + 3 * price(is_call, S, SA, X, T, T2, r, b, v, cnd=cnd)
        - price(is_call, S - dS, SA, X, T, T2, r, b, v, cnd=cnd)
    )

def strike_delta(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X + dS, T, T2, r, b, v, cnd=cnd)
        - price(is_call, S, SA, X - dS, T, T2, r, b, v, cnd=cnd)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, X + dS, T, T2, r, b, v, cnd=cnd)
        - 2 * price(is_call, S, SA, X, T, T2, r, b, v, cnd=cnd)
        + price(is_call, S, SA, X - dS, T, T2, r, b, v, cnd=cnd)
    ) / dS ** 2
