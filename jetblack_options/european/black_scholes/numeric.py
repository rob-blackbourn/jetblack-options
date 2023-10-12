"""Plain Vanilla"""

from math import exp, log, pi, sqrt
from typing import Callable, Literal, Optional

from ...distributions import CND, ND, CNDEV, CHIINV

from .analytic import price

def delta(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd)
    ) / (2 * dS)

def deltap(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S * (1 + dS), X, T, r, b, v, cnd=cnd) -
        price(is_call, S * (1 - dS), X, T, r, b, v, cnd=cnd)
    ) * 2 / S

def elasticity(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd)
    ) / (2 * dS) * S / price(is_call, S, X, T, r, b, v, cnd=cnd)

def gamma(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd)
    ) / dS ** 2

def dgamma_dvol(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v + dv, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v + dv, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v + dv, cnd=cnd) -
        price(is_call, S + dS, X, T, r, b, v - dv, cnd=cnd) +
        2 * price(is_call, S, X, T, r, b, v - dv, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v - dv, cnd=cnd)
    ) / (2 * dv * dS ** 2) / 100

def gammap(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return S / 100 * (
        price(is_call, S + dS, X, T, r, b, v, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd)
    ) / dS ** 2

def ddelta_dvol(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v + dv, cnd=cnd) -
        price(is_call, S + dS, X, T, r, b, v - dv, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v + dv, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v - dv, cnd=cnd)
    ) / (4 * dS * dv) / 100

def vega(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, cnd=cnd) -
        price(is_call, S, X, T, r, b, v - dv, cnd=cnd)
    ) / 2

def vegap(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, cnd=cnd) -
        price(is_call, S, X, T, r, b, v - dv, cnd=cnd)
    ) * v / 0.1 / 2

def dvega_dvol(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd) +
        price(is_call, S, X, T, r, b, v - dv, cnd=cnd)
    )

def theta(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dt: float = 1.0 / 365.0,
        cnd: Callable[[float], float] = CND
) -> float:
    if T <= dt:
        return (
            price(is_call, S, X, 0.00001, r, b, v, cnd=cnd) -
            price(is_call, S, X, T, r, b, v, cnd=cnd)
        )
    else:
        return (
            price(is_call, S, X, T - dt, r, b, v, cnd=cnd) -
            price(is_call, S, X, T, r, b, v, cnd=cnd)
        )

def rho(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r + dr, b + dr, v, cnd=cnd) -
        price(is_call, S, X, T, r - dr, b - dr, v, cnd=cnd)
    ) / 2

def futures_rho(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r + dr, 0, v, cnd=cnd) -
        price(is_call, S, X, T, r - dr, 0, v, cnd=cnd)
    ) / 2

def futures_rho2(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b - db, v, cnd=cnd) -
        price(is_call, S, X, T, r, b + db, v, cnd=cnd)
    ) / 2

def carry(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b + db, v, cnd=cnd) -
        price(is_call, S, X, T, r, b - db, v, cnd=cnd)
    ) / 2

def speed(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + 2 * dS, X, T, r, b, v, cnd=cnd) -
        3 * price(is_call, S + dS, X, T, r, b, v, cnd=cnd) +
        3 * price(is_call, S, X, T, r, b, v, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd)
    ) / dS ** 3

def strike_delta(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dX: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X + dX, T, r, b, v, cnd=cnd) -
        price(is_call, S, X - dX, T, r, b, v, cnd=cnd)
    ) / (2 * dX)

def strike_gamma(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dX: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X + dX, T, r, b, v, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd) +
        price(is_call, S, X - dX, T, r, b, v, cnd=cnd)
    ) / dX ** 2
