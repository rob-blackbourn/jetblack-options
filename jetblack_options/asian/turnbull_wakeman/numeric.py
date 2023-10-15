"""Turnbull-Wakeman numeric"""

from typing import Callable

from ...distributions import CND

from .analytic import price

def delta(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, T, T2, r, b, v, cdf=cdf) -
        price(is_call, S - dS, SA, K, T, T2, r, b, v, cdf=cdf)
    ) / (2 * dS)


def gamma(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, T, T2, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, SA, K, T, T2, r, b, v, cdf=cdf)
        + price(is_call, S - dS, SA, K, T, T2, r, b, v, cdf=cdf)
    ) / dS ** 2


def theta(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365,
        cdf: Callable[[float], float] = CND
) -> float:
    if T <= dT:
        return (
            price(is_call, S, SA, K, 0.00001, T2, r, b, v, cdf=cdf)
            - price(is_call, S, SA, K, T, T2, r, b, v, cdf=cdf)
        )
    else:
        return (
            price(is_call, S, SA, K, T - dT, T2, r, b, v, cdf=cdf)
            - price(is_call, S, SA, K, T, T2, r, b, v, cdf=cdf)
        )


def vega(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, T, T2, r, b, v + dv, cdf=cdf)
        - price(is_call, S, SA, K, T, T2, r, b, v - dv, cdf=cdf)
    ) / 2


def rho(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, T, T2, r + dr, b + dr, v, cdf=cdf)
        - price(is_call, S, SA, K, T, T2, r - dr, b - dr, v, cdf=cdf)
    ) / 2


def carry(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, T, T2, r, b + db, v, cdf=cdf)
        - price(is_call, S, SA, K, T, T2, r, b - db, v, cdf=cdf)
    ) / 2


def elasticity(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, T, T2, r, b, v, cdf=cdf)
        - price(is_call, S - dS, SA, K, T, T2, r, b, v, cdf=cdf)
    ) / (2 * dS) * S / price(is_call, S, SA, K, T, T2, r, b, v, cdf=cdf)


def speed(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return 1 / dS ** 3 * (
        price(is_call, S + 2 * dS, SA, K, T, T2, r, b, v, cdf=cdf)
        - 3 * price(is_call, S + dS, SA, K, T, T2, r, b, v, cdf=cdf)
        + 3 * price(is_call, S, SA, K, T, T2, r, b, v, cdf=cdf)
        - price(is_call, S - dS, SA, K, T, T2, r, b, v, cdf=cdf)
    )


def gammap(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return S / 100 * (
        price(is_call, S + dS, SA, K, T, T2, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, SA, K, T, T2, r, b, v, cdf=cdf)
        + price(is_call, S - dS, SA, K, T, T2, r, b, v, cdf=cdf)
    ) / dS ** 2


def vegap(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, T, T2, r, b, v + dv, cdf=cdf)
        - price(is_call, S, SA, K, T, T2, r, b, v - dv, cdf=cdf)
    ) * v / 0.1 / 2


def ddelta_dvol(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, T, T2, r, b, v + dv, cdf=cdf)
        -  price(is_call, S + dS, SA, K, T, T2, r, b, v - dv, cdf=cdf)
        - price(is_call, S - dS, SA, K, T, T2, r, b, v + dv, cdf=cdf)
        + price(is_call, S - dS, SA, K, T, T2, r, b, v - dv, cdf=cdf)
    ) / (4 * dS * 0.01) / 100


def dgamma_dvol(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, SA, K, T, T2, r, b, v + dv, cdf=cdf)
        - 2 * price(is_call, S, SA, K, T, T2, r, b, v + dv, cdf=cdf)
        + price(is_call, S - dS, SA, K, T, T2, r, b, v + dv, cdf=cdf)
        - price(is_call, S + dS, SA, K, T, T2, r, b, v - dv, cdf=cdf)
        + 2 * price(is_call, S, SA, K, T, T2, r, b, v - dv, cdf=cdf)
        - price(is_call, S - dS, SA, K, T, T2, r, b, v - dv, cdf=cdf)
    ) / (2 * dv * dS ** 2) / 100


def dvega_dvol(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, T, T2, r, b, v + dv, cdf=cdf)
        - 2 * price(is_call, S, SA, K, T, T2, r, b, v, cdf=cdf)
        + price(is_call, S, SA, K, T, T2, r, b, v - dv, cdf=cdf)
    )


def vomma(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return dvega_dvol(is_call, S, SA, K, T, T2, r, b, v, dv=dv, cdf=cdf) / dv ** 2 / 10000

def futures_rho(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, T, T2, r + dr, b, v, cdf=cdf)
        - price(is_call, S, SA, K, T, T2, r - dr, b, v, cdf=cdf)
    ) / 2

def rho2(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K, T, T2, r, b - db, v, cdf=cdf)
        - price(is_call, S, SA, K, T, T2, r, b + db, v, cdf=cdf)
    ) / 2

def strike_delta(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K + dS, T, T2, r, b, v, cdf=cdf)
        - price(is_call, S, SA, K - dS, T, T2, r, b, v, cdf=cdf)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        S: float,
        SA: float,
        K: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, SA, K + dS, T, T2, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, SA, K, T, T2, r, b, v, cdf=cdf)
        + price(is_call, S, SA, K - dS, T, T2, r, b, v, cdf=cdf)
    ) / dS ** 2
