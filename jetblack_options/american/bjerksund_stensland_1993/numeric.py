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
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S + dS, X, T, r, b, v + 0.01, cnd=cnd) -
        2 * price(is_call, S, X, T, r, b, v + 0.01, cnd=cnd) +
        price(is_call, S - dS, X, T, r, b, v + 0.01, cnd=cnd) -
        price(is_call, S + dS, X, T, r, b, v - 0.01, cnd=cnd) +
        2 * price(is_call, S, X, T, r, b, v - 0.01, cnd=cnd) -
        price(is_call, S - dS, X, T, r, b, v - 0.01, cnd=cnd)
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
        cnd: Callable[[float], float] = CND
) -> float:
        return S / 100 * (
            price(is_call, S + dS, X, T, r, b, v, cnd=cnd)
            - 2 * price(is_call, S, X, T, r, b, v, cnd=cnd)
            + price(is_call, S - dS, X, T, r, b, v, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T + 1 / 365, r, b, v, cnd=cnd)
        - 2 * price(is_call, S, X, T, r, b, v, cnd=cnd)
        + price(is_call, S, X, T - 1 / 365, r, b, v, cnd=cnd)
    ) / (1 / 365) ** 2

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
        cnd: Callable[[float], float] = CND
) -> float:
    return 1 / (4 * dS * 0.01) * (
        price(is_call, S + dS, X, T, r, b, v + 0.01, cnd=cnd)
        - price(is_call, S + dS, X, T, r, b, v - 0.01, cnd=cnd)
        - price(is_call, S - dS, X, T, r, b, v + 0.01, cnd=cnd)
        + price(is_call, S - dS, X, T, r, b, v - 0.01, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + 0.01, cnd=cnd)
        - price(is_call, S, X, T, r, b, v - 0.01, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b, v + 0.01, cnd=cnd)
        - 2 * price(is_call, S, X, T, r, b, v, cnd=cnd)
        + price(is_call, S, X, T, r, b, v - 0.01, cnd=cnd)
    ) / 0.01 ** 2 / 10000

def vegap(
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
        return v / 0.1 * (
            price(is_call, S, X, T, r, b, v + 0.01, cnd=cnd)
            - price(is_call, S, X, T, r, b, v - 0.01, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
        return (
            price(is_call, S, X, T, r, b, v + 0.01, cnd=cnd)
            - 2 * price(is_call, S, X, T, r, b, v, cnd=cnd)
            + price(is_call, S, X, T, r, b, v - 0.01, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    if T <= 1 / 365:
        return (
            price(is_call, S, X, 0.00001, r, b, v, cnd=cnd)
            - price(is_call, S, X, T, r, b, v, cnd=cnd)
        )
    else:
        return (
            price(is_call, S, X, T - 1 / 365, r, b, v, cnd=cnd)
            - price(is_call, S, X, T, r, b, v, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r + 0.01, b + 0.01, v, cnd=cnd)
        - price(is_call, S, X, T, r - 0.01, b - 0.01, v, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r + 0.01, b, v, cnd=cnd)
        - price(is_call, S, X, T, r - 0.01, b, v, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b - 0.01, v, cnd=cnd)
        - price(is_call, S, X, T, r, b + 0.01, v, cnd=cnd)
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
        dS: float = 0.01,
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X, T, r, b + 0.01, v, cnd=cnd)
        - price(is_call, S, X, T, r, b - 0.01, v, cnd=cnd)
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
    return 1 / dS ** 3 * (
        price(is_call, S + 2 * dS, X, T, r, b, v, cnd=cnd)
        - 3 * price(is_call, S + dS, X, T, r, b, v, cnd=cnd)
        + 3 * price(is_call, S, X, T, r, b, v, cnd=cnd)
        - price(is_call, S - dS, X, T, r, b, v, cnd=cnd)
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
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X + dS, T, r, b, v, cnd=cnd)
        - price(is_call, S, X - dS, T, r, b, v, cnd=cnd)
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
        cnd: Callable[[float], float] = CND
) -> float:
    return (
        price(is_call, S, X + dS, T, r, b, v, cnd=cnd)
        - 2 * price(is_call, S, X, T, r, b, v, cnd=cnd)
        + price(is_call, S, X - dS, T, r, b, v, cnd=cnd)
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
        price(is_call, S, X, T, r, b, v, cnd=cnd)
        - bs_price(is_call, S, X, T, r, b, v, cnd=cnd)
    )

def ivol(
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
    is_call = S >= S * exp(b * T)
    p = price(is_call, S, X, T, r, b, v, cnd=cnd)
    return bs_ivol(is_call, S, X, T, r, b, p, cnd=cnd)

