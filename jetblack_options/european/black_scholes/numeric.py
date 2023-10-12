"""Plain Vanilla"""

from math import exp, log, pi, sqrt
from typing import Literal, Optional

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
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v) -
        price(is_call, S - dS, X, T, r, b, v)
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
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S * (1 + dS), X, T, r, b, v) -
        price(is_call, S * (1 - dS), X, T, r, b, v)
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
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v) -
        price(is_call, S - dS, X, T, r, b, v)
    ) / (2 * dS) * S / price(is_call, S, X, T, r, b, v)

def gamma(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v) -
        2 * price(is_call, S, X, T, r, b, v) +
        price(is_call, S - dS, X, T, r, b, v)
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
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v + 0.01) -
        2 * price(is_call, S, X, T, r, b, v + 0.01) +
        price(is_call, S - dS, X, T, r, b, v + 0.01) -
        price(is_call, S + dS, X, T, r, b, v - 0.01) +
        2 * price(is_call, S, X, T, r, b, v - 0.01) -
        price(is_call, S - dS, X, T, r, b, v - 0.01)
    ) / (2 * 0.01 * dS ** 2) / 100

def gammap(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return S / 100 * (
        price(is_call, S + dS, X, T, r, b, v) -
        2 * price(is_call, S, X, T, r, b, v) +
        price(is_call, S - dS, X, T, r, b, v)
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
        dv: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v + dv) -
        price(is_call, S + dS, X, T, r, b, v - dv) -
        price(is_call, S - dS, X, T, r, b, v + dv) +
        price(is_call, S - dS, X, T, r, b, v - dv)
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
        dv: float = 0.01
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv) -
        price(is_call, S, X, T, r, b, v - dv)
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
        dv: float = 0.01
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv) -
        price(is_call, S, X, T, r, b, v - dv)
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
        dv: float = 0.01
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv) -
        2 * price(is_call, S, X, T, r, b, v) +
        price(is_call, S, X, T, r, b, v - dv)
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
        dt: float = 1.0 / 365.0
) -> float:
    if T <= dt:
        return (
            price(is_call, S, X, 0.00001, r, b, v) -
            price(is_call, S, X, T, r, b, v)
        )
    else:
        return (
            price(is_call, S, X, T - dt, r, b, v) -
            price(is_call, S, X, T, r, b, v)
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
        dr: float = 0.01
) -> float:
    return (
        price(is_call, S, X, T, r + dr, b + dr, v) -
        price(is_call, S, X, T, r - dr, b - dr, v)
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
        dr: float = 0.01
) -> float:
    return (
        price(is_call, S, X, T, r + dr, 0, v) -
        price(is_call, S, X, T, r - dr, 0, v)
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
        db: float = 0.01
) -> float:
    return (
        price(is_call, S, X, T, r, b - db, v) -
        price(is_call, S, X, T, r, b + db, v)
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
        db: float = 0.01
) -> float:
    return (
        price(is_call, S, X, T, r, b + db, v) -
        price(is_call, S, X, T, r, b - db, v)
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
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + 2 * dS, X, T, r, b, v) -
        3 * price(is_call, S + dS, X, T, r, b, v) +
        3 * price(is_call, S, X, T, r, b, v) -
        price(is_call, S - dS, X, T, r, b, v)
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
        dX: float = 0.01
) -> float:
    return (
        price(is_call, S, X + dX, T, r, b, v) -
        price(is_call, S, X - dX, T, r, b, v)
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
        dX: float = 0.01
) -> float:
    return (
        price(is_call, S, X + dX, T, r, b, v) -
        2 * price(is_call, S, X, T, r, b, v) +
        price(is_call, S, X - dX, T, r, b, v)
    ) / dX ** 2
