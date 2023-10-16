"""Haug, Haug & Margrabe - numeric"""

from math import exp, log, sqrt
from typing import Callable

from ...distributions import CDF

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) / (2 * dS)

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) / (2 * dS) * S / price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)

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
        cdf: Callable[[float], float] = CDF
) -> float:
        return (
            price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
            - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
            + price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        ) / dS ** 2

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        + price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
        + 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) / (2 * dv * dS ** 2) / 100

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) * S / 100 / dS ** 2

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        + price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) / (4 * dS * dv) / 100

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) / 2

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) / dv ** 2 / 10000

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    ) *v / 0.1 / 2

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b, v + dv, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + price(is_call, S, SA, K, t1, T, n, m, r, b, v - dv, cdf=cdf)
    )

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
        cdf: Callable[[float], float] = CDF
) -> float:
    if t1 > dT and T > dT:
        return (
            price(is_call, S, SA, K, t1 - dT, T - dT, n, m, r, b, v, cdf=cdf)
            - price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        )
    else:
        raise ValueError("insufficient time to expiry")

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r + dr, b + dr, v, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r - dr, b - dr, v, cdf=cdf)
    ) / 2

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
        cdf: Callable[[float], float] = CDF
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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b - db, v, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r, b + db, v, cdf=cdf)
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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S, SA, K, t1, T, n, m, r, b + db, v, cdf=cdf)
        - price(is_call, S, SA, K, t1, T, n, m, r, b - db, v, cdf=cdf)
    ) / 2

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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S + 2 * dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - 3 * price(is_call, S + dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + 3 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        - price(is_call, S - dS, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
    ) / dS ** 3

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
        cdf: Callable[[float], float] = CDF
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
        cdf: Callable[[float], float] = CDF
) -> float:
    return (
        price(is_call, S, SA, K + dS, t1, T, n, m, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, SA, K, t1, T, n, m, r, b, v, cdf=cdf)
        + price(is_call, S, SA, K - dS, t1, T, n, m, r, b, v, cdf=cdf)
    ) / dS ** 2
