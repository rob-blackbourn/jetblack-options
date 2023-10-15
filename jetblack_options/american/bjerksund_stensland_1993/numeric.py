"""American"""

from math import exp
from typing import Callable

from ...distributions import CND
from ...european.black_scholes.analytic import price as bs_price
from ...european.black_scholes.implied_volatility import ivol as bs_ivol
from .analytic import price

def delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
        return (
            price(is_call, S + dS, K, T, r, b, v, cdf=cdf) -
            price(is_call, S - dS, K, T, r, b, v, cdf=cdf)
        ) / (2 * dS)

def gamma(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, K, T, r, b, v, cdf=cdf) -
        2 * price(is_call, S, K, T, r, b, v, cdf=cdf) +
        price(is_call, S - dS, K, T, r, b, v, cdf=cdf)
    ) / dS ** 2

def theta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365,
        cdf: Callable[[float], float] = CND
) -> float:
    if T <= dT:
        return (
            price(is_call, S, K, 0.00001, r, b, v, cdf=cdf)
            - price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
    else:
        return (
            price(is_call, S, K, T - dT, r, b, v, cdf=cdf)
            - price(is_call, S, K, T, r, b, v, cdf=cdf)
        )

def vega(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, K, T, r, b, v + dv, cdf=cdf)
        - price(is_call, S, K, T, r, b, v - dv, cdf=cdf)
    ) / 2

def rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, K, T, r + dr, b + dr, v, cdf=cdf)
        - price(is_call, S, K, T, r - dr, b - dr, v, cdf=cdf)
    ) / 2

def carry(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, K, T, r, b + db, v, cdf=cdf)
        - price(is_call, S, K, T, r, b - db, v, cdf=cdf)
    ) / 2

def elasticity(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, K, T, r, b, v, cdf=cdf) -
        price(is_call, S - dS, K, T, r, b, v, cdf=cdf)
    ) / (2 * dS) * S / price(is_call, S, K, T, r, b, v, cdf=cdf)

def speed(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return 1 / dS ** 3 * (
        price(is_call, S + 2 * dS, K, T, r, b, v, cdf=cdf)
        - 3 * price(is_call, S + dS, K, T, r, b, v, cdf=cdf)
        + 3 * price(is_call, S, K, T, r, b, v, cdf=cdf)
        - price(is_call, S - dS, K, T, r, b, v, cdf=cdf)
    )

def gammap(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
        return S / 100 * (
            price(is_call, S + dS, K, T, r, b, v, cdf=cdf)
            - 2 * price(is_call, S, K, T, r, b, v, cdf=cdf)
            + price(is_call, S - dS, K, T, r, b, v, cdf=cdf)
        ) / dS ** 2

def vegap(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
        return v / 0.1 * (
            price(is_call, S, K, T, r, b, v + dv, cdf=cdf)
            - price(is_call, S, K, T, r, b, v - dv, cdf=cdf)
        ) / 2

def ddelta_dvol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return 1 / (4 * dS * 0.01) * (
        price(is_call, S + dS, K, T, r, b, v + 0.01, cdf=cdf)
        - price(is_call, S + dS, K, T, r, b, v - 0.01, cdf=cdf)
        - price(is_call, S - dS, K, T, r, b, v + 0.01, cdf=cdf)
        + price(is_call, S - dS, K, T, r, b, v - 0.01, cdf=cdf)
    ) / 100

def dgamma_dvol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, K, T, r, b, v + 0.01, cdf=cdf) -
        2 * price(is_call, S, K, T, r, b, v + 0.01, cdf=cdf) +
        price(is_call, S - dS, K, T, r, b, v + 0.01, cdf=cdf) -
        price(is_call, S + dS, K, T, r, b, v - 0.01, cdf=cdf) +
        2 * price(is_call, S, K, T, r, b, v - 0.01, cdf=cdf) -
        price(is_call, S - dS, K, T, r, b, v - 0.01, cdf=cdf)
    ) / (2 * 0.01 * dS ** 2) / 100

def dvega_dvol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
        return (
            price(is_call, S, K, T, r, b, v + dv, cdf=cdf)
            - 2 * price(is_call, S, K, T, r, b, v, cdf=cdf)
            + price(is_call, S, K, T, r, b, v - dv, cdf=cdf)
        )

def vomma(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return dvega_dvol(is_call, S, K, T, r, b, v, dv=dv, cdf=cdf) / dv ** 2 / 10000

def time_gamma(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, K, T + dT, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, K, T, r, b, v, cdf=cdf)
        + price(is_call, S, K, T - dT, r, b, v, cdf=cdf)
    ) / dT ** 2

def futures_rho(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, K, T, r + dr, b, v, cdf=cdf)
        - price(is_call, S, K, T, r - dr, b, v, cdf=cdf)
    ) / 2

def rho2(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, K, T, r, b - db, v, cdf=cdf)
        - price(is_call, S, K, T, r, b + db, v, cdf=cdf)
    ) / 2

def strike_delta(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, K + dS, T, r, b, v, cdf=cdf)
        - price(is_call, S, K - dS, T, r, b, v, cdf=cdf)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01,
        cdf: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, K + dS, T, r, b, v, cdf=cdf)
        - 2 * price(is_call, S, K, T, r, b, v, cdf=cdf)
        + price(is_call, S, K - dS, T, r, b, v, cdf=cdf)
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
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v, cdf=cnd)
        - bs_price(is_call, S, X, T, r, b, v, cdf=cnd)
    )

def ivol(
        is_call: bool,
        S: float,
        K: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        cdf: Callable[[float], float] = CND
) -> float:
    is_call = S >= S * exp(b * T)
    p = price(is_call, S, K, T, r, b, v, cdf=cdf)
    return bs_ivol(is_call, S, K, T, r, b, p, cdf=cdf)

