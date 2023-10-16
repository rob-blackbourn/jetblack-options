"""Asian options"""

from math import exp, log, sqrt
from typing import Callable, Optional

from ...distributions import CND
from ...european.black_scholes_merton import price as bs_price

from .analytic import price

def delta(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) / (2 * dS)


def gamma(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) / dS ** 2


def theta(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365,
        cdf: Callable[[float], float] = CND
) -> float:
    if t1 > dT and T > dT:
        return (
            price(is_call, S, SA, K, t1 - dT, T - dT, n, m, r, b, v, cdf=cdf)
            - price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        )
    else:
        raise ValueError("must have at least 1 day")


def vega(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) / 2


def rho(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r + dr, b + dr, v, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r - dr, b - dr, v, cdf=cdf)
    ) / 2

def carry(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b + db, v, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r, b - db, v, cdf=cdf)
    ) / 2


def elasticity(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) / (2 * dS) * S / price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)


def speed(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + 2 * dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - 3 * price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + 3 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) / dS ** 3


def gammap(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) *S / 100 / dS ** 2


def vegap(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return v / 0.1 * (
        price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) / 2


def ddelta_dvol(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        + price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) / (4 * dS * dv) / 100


def dgamma_dvol(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        + price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
        + 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) / (2 * dv * dS ** 2) / 100


def dvega_dvol(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    )


def vomma(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return dvega_dvol(is_call, S, SA, K, t1, T, n, m, r, b, v, dv=dv, cdf=cdf)  / dv ** 2 / 10000

def futures_rho(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r + dr, b, v, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r - dr, b, v, cdf=cdf)
    ) / 2


def rho2(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b - db, v, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r, b + db, v, cdf=cdf)
    ) / 2

def strike_delta(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K + dS, t1, T, n, m, r, b, v, cdf=cdf)
        - price(is_call, S, SA, K - dS, t1, T, n, m, r, b, v, cdf=cdf)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K + dS, t1, T, n, m, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + price(is_call, S, SA, K - dS, t1, T, n, m, r, b, v, cdf=cdf)
    ) / dS ** 2
