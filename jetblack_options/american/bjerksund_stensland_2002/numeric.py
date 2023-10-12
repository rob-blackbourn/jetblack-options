"""American"""

from typing import Callable

from ...distributions import CND, CBND
from ...european.black_scholes.analytic import price as bs_price

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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / (2 * dS)

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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / (2 * dS) * S / price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd)


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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd) +
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd)
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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v + 0.01, cnd=cnd, cbnd=cbnd) -
        2 * price(is_call, S, X, T, r, b, v + 0.01, cnd=cnd, cbnd=cbnd) +
        price(is_call, S - dS, X, T, r, b, v + 0.01, cnd=cnd, cbnd=cbnd) -
        price(is_call, S + dS, X, T, r, b, v - 0.01, cnd=cnd, cbnd=cbnd) +
        2 * price(is_call, S, X, T, r, b, v - 0.01, cnd=cnd, cbnd=cbnd) -
        price(is_call, S - dS, X, T, r, b, v - 0.01, cnd=cnd, cbnd=cbnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return S / 100 * (
        price(is_call, S + dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd) +
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / dS ** 2


def time_gamma(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X, T + dT, r, b, v, cnd=cnd, cbnd=cbnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd) +
        price(is_call, S, X, T - dT, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / dT ** 2


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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return 1 / (4 * dS * 0.01) * (
        price(is_call, S + dS, X, T, r, b, v + dv, cnd=cnd, cbnd=cbnd) -
        price(is_call, S + dS, X, T, r, b, v - dv, cnd=cnd, cbnd=cbnd) -
        price(is_call, S - dS, X, T, r, b, v + dv, cnd=cnd, cbnd=cbnd) +
        price(is_call, S - dS, X, T, r, b, v - dv, cnd=cnd, cbnd=cbnd)
    ) / 100


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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, cnd=cnd, cbnd=cbnd) -
        price(is_call, S, X, T, r, b, v - dv, cnd=cnd, cbnd=cbnd)
    ) / 2


def vomma(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, cnd=cnd, cbnd=cbnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd) +
        price(is_call, S, X, T, r, b, v - dv, cnd=cnd, cbnd=cbnd)
    ) / dv ** 2 / 10000


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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return v / 0.1 * (
        price(is_call, S, X, T, r, b, v + dv, cnd=cnd, cbnd=cbnd) -
        price(is_call, S, X, T, r, b, v - dv, cnd=cnd, cbnd=cbnd)
    ) / 2


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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, cnd=cnd, cbnd=cbnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd) +
        price(is_call, S, X, T, r, b, v - dv, cnd=cnd, cbnd=cbnd)
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
        dT: float = 1 / 365,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    if T <= dT:
        return (
            price(is_call, S, X, 0.00001, r, b, v, cnd=cnd, cbnd=cbnd) -
            price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd)
        )
    else:
        return (
            price(is_call, S, X, T - dT, r, b, v, cnd=cnd, cbnd=cbnd) -
            price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd)
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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X, T, r + dr, b + dr, v, cnd=cnd, cbnd=cbnd) -
        price(is_call, S, X, T, r - dr, b - dr, v, cnd=cnd, cbnd=cbnd)
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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X, T, r + dr, b, v, cnd=cnd, cbnd=cbnd) -
        price(is_call, S, X, T, r - dr, b, v, cnd=cnd, cbnd=cbnd)
    ) / 2

def rho2(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X, T, r, b - db, v, cnd=cnd, cbnd=cbnd) -
        price(is_call, S, X, T, r, b + db, v, cnd=cnd, cbnd=cbnd)
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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X, T, r, b + db, v, cnd=cnd, cbnd=cbnd) -
        price(is_call, S, X, T, r, b - db, v, cnd=cnd, cbnd=cbnd)
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
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return 1 / dS ** 3 * (
        price(is_call, S + 2 * dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        3 * price(is_call, S + dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd) +
        3 * price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        price(is_call, S - dS, X, T, r, b, v, cnd=cnd, cbnd=cbnd)
    )

def strike_delta(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X + dS, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        price(is_call, S, X - dS, T, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float:
    return (
        price(is_call, S, X + dS, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        2 * price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd) +
        price(is_call, S, X - dS, T, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / dS ** 2

def price_diff(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND,
) -> float: # Difference in value between BS Approx and Black-Scholes Merton value
    return (
        price(is_call, S, X, T, r, b, v, cnd=cnd, cbnd=cbnd) -
        bs_price(is_call, S, X, T, r, b, v, cnd=cnd)
    )

