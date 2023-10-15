"""Barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, Callable

from ...distributions import CND, CBND

from .analytic import price

def delta(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S + dS, X, L, U, T, r, b, v, cnd=cnd)
        - price(is_call, is_in, S - dS, X, L, U, T, r, b, v, cnd=cnd)
    ) / (2 * dS)


def ddelta_dvol(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S + dS, X, L, U, T, r, b, v + dv, cnd=cnd)
        - price(is_call, is_in, S + dS, X, L, U, T, r, b, v - dv, cnd=cnd)
        - price(is_call, is_in, S - dS, X, L, U, T, r, b, v + dv, cnd=cnd)
        + price(is_call, is_in, S - dS, X, L, U, T, r, b, v - dv, cnd=cnd)
    ) / (4 * dS * dv) / 100

def gamma(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S + dS, X, L, U, T, r, b, v, cnd=cnd)
        - 2 * price(is_call, is_in, S, X, L, U, T, r, b, v, cnd=cnd)
        + price(is_call, is_in, S - dS, X, L, U, T, r, b, v, cnd=cnd)
    ) / (dS ** 2)

def gammap(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND
) -> float:
    return S / 100 * gamma(is_call, is_in, S + dS, X, L, U, T, r, b, v, cnd=cnd)

def dgamma_dvol(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        dv=0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S + dS, X, L, U, T, r, b, v + dv, cnd=cnd)
        - 2 * price(is_call, is_in, S, X, L, U, T, r, b, v + dv, cnd=cnd)
        + price(is_call, is_in, S - dS, X, L, U, T, r, b, v + dv, cnd=cnd)
        - price(is_call, is_in, S + dS, X, L, U, T, r, b, v - dv, cnd=cnd)
        + 2 * price(is_call, is_in, S, X, L, U, T, r, b, v - dv, cnd=cnd)
        - price(is_call, is_in, S - dS, X, L, U, T, r, b, v - dv, cnd=cnd)
    ) / (2 * dv * dS ** 2) / 100

def vega(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S, X, L, U, T, r, b, v + dv, cnd=cnd)
        - price(is_call, is_in, S, X, L, U, T, r, b, v - dv, cnd=cnd)
    ) / 2

def vomma(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S, X, L, U, T, r, b, v + dv)
        - 2 * price(is_call, is_in, S, X, L, U, T, r, b, v)
        + price(is_call, is_in, S, X, L, U, T, r, b, v - dv)
    ) / dv ** 2 / 10000

def vegap(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        vega(is_call, is_in, S + dS, X, L, U, T, r, b, v, dv=dv, cnd=cnd)
        * v / 0.1
    )

def rho(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S, X, L, U, T, r + dr, b + dr, v, cnd=cnd)
        - price(is_call, is_in, S, X, L, U, T, r - dr, b - dr, v, cnd=cnd)
    ) / 2

def rho2(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S, X, L, U, T, r, b - db, v, cnd=cnd)
        - price(is_call, is_in, S, X, L, U, T, r, b + db, v, cnd=cnd)
    ) / 2

def carry(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S, X, L, U, T, r, b + db, v, cnd=cnd)
        - price(is_call, is_in, S, X, L, U, T, r, b - db, v, cnd=cnd)
    ) / 2

def theta(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365,
        cnd: Callable[[float], float] = CND
) -> float:
    if T <= dT:
        return (
            price(is_call, is_in, S, X, L, U, 0.00001, r, b, v, cnd=cnd)
            - price(is_call, is_in, S, X, L, U, T, r, b, v, cnd=cnd)
        )
    else:
        return (
            price(is_call, is_in, S, X, L, U, T - dT, r, b, v, cnd=cnd)
            - price(is_call, is_in, S, X, L, U, T, r, b, v, cnd=cnd)
        )

def strike_delta(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S, X + dS, L, U, T, r, b, v, cnd=cnd)
        - price(is_call, is_in, S, X - dS, L, U, T, r, b, v, cnd=cnd)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, is_in, S, X + dS, L, U, T, r, b, v, cnd=cnd)
        - 2 * price(is_call, is_in, S, X, L, U, T, r, b, v, cnd=cnd)
        + price(is_call, is_in, S, X - dS, L, U, T, r, b, v, cnd=cnd)
    ) / (dS ** 2)