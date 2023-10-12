"""American"""

from math import exp, log, sqrt
from typing import Callable, Literal, Optional

from ...distributions import CND, ND
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
            price(is_call, S + dS, X, T, r, b, v, nd=nd, cnd=cnd) -
            price(is_call, S - dS, X, T, r, b, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v, nd=nd, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v, nd=nd, cnd=cnd)
    ) / (2 * dS) * S / price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd)

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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v, nd=nd, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v + 0.01, nd=nd, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v + 0.01, nd=nd, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v + 0.01, nd=nd, cnd=cnd) -
        price(is_call, S + dS, X, T, r, b, v - 0.01, nd=nd, cnd=cnd) +
        2 * price(is_call, S, X, T, r, b, v - 0.01, nd=nd, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v - 0.01, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return S / 100 * (
        price(is_call, S + dS, X, T, r, b, v, nd=nd, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T + dT, r, b, v, nd=nd, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd) +
        price(is_call, S, X, T - dT, r, b, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return 1 / (4 * dS * 0.01) * (
        price(is_call, S + dS, X, T, r, b, v + 0.01, nd=nd, cnd=cnd) -
        price(is_call, S + dS, X, T, r, b, v - 0.01, nd=nd, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v + 0.01, nd=nd, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v - 0.01, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, nd=nd, cnd=cnd) -
        price(is_call, S, X, T, r, b, v - dv, nd=nd, cnd=cnd)
    ) / 2

def dvega_dvolXX(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float: # DvegaDvol/vomma
    return (
        price(is_call, S, X, T, r, b, v + dv) -
        2 * price(is_call, S, X, T, r, b, v) +
        price(is_call, S, X, T, r, b, v - dv)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, nd=nd, cnd=cnd) -
        price(is_call, S, X, T, r, b, v - dv, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + dv, nd=nd, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd) +
        price(is_call, S, X, T, r, b, v - dv, nd=nd, cnd=cnd)
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
        dT: float = 1 /365,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    if T <= dT:
        return (
            price(is_call, S, X, 0.00001, r, b, v, nd=nd, cnd=cnd)
            - price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd)
        )
    else:
        return (
            price(is_call, S, X, T - dT, r, b, v, nd=nd, cnd=cnd)
            - price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r + dr, b + dr, v, nd=nd, cnd=cnd)
        - price(is_call, S, X, T, r - dr, b - dr, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r + dr, b, v, nd=nd, cnd=cnd)
        - price(is_call, S, X, T, r - dr, b, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b - db, v, nd=nd, cnd=cnd)
        - price(is_call, S, X, T, r, b + db, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b + db, v, nd=nd, cnd=cnd)
        - price(is_call, S, X, T, r, b - db, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + 2 * dS, X, T, r, b, v, nd=nd, cnd=cnd)
        - 3 * price(is_call, S + dS, X, T, r, b, v, nd=nd, cnd=cnd)
        + 3 * price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd)
        - price(is_call, S - dS, X, T, r, b, v, nd=nd, cnd=cnd)
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
        dS: float = 0.01,
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X + dS, T, r, b, v, nd=nd, cnd=cnd)
        - price(is_call, S, X - dS, T, r, b, v, nd=nd, cnd=cnd)
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
        nd: Callable[[float], float] = ND,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X + dS, T, r, b, v, nd=nd, cnd=cnd)
        - 2 * price(is_call, S, X, T, r, b, v, nd=nd, cnd=cnd)
        + price(is_call, S, X - dS, T, r, b, v, nd=nd, cnd=cnd)
    ) / dS ** 2

    # elif OutPutFlag == "di": #Difference in value between BS Approx and Black-Scholes Merton value
    #     return (
    #         price(is_call, S, X, T, r, b, v) -
    #         price(is_call, S, X, T, r, b, v)
    #     )
