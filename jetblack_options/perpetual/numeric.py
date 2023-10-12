"""Perpetual numeric solutions"""

from math import exp, log, sqrt
from typing import Callable, Literal, Optional

from ..distributions import CND, CBND, ND
from ..european.black_scholes.analytic import price as bs_price
from ..european.black_scholes.implied_volatility import ivol

from .analytic import price

def delta(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
        return (
             price(is_call, S + dS, X, r, b, v) -
             price(is_call, S - dS, X, r, b, v)
        ) / (2 * dS)

def elasticity(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, X, r, b, v) -
        price(is_call, S - dS, X, r, b, v)
    ) / (2 * dS) * S / price(is_call, S, X, r, b, v)

def gamma(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, X, r, b, v) -
        2 * price(is_call, S, X, r, b, v) +
        price(is_call, S - dS, X, r, b, v)
    ) / dS ** 2

def dgamma_dvol(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S + dS, X, r, b, v + 0.01) -
        2 * price(is_call, S, X, r, b, v + 0.01) +
        price(is_call, S - dS, X, r, b, v + 0.01) -
        price(is_call, S + dS, X, r, b, v - 0.01) +
        2 * price(is_call, S, X, r, b, v - 0.01) -
        price(is_call, S - dS, X, r, b, v - 0.01)
    ) / (2 * 0.01 * dS ** 2) / 100

def gammap(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return S / 100 * (
        price(is_call, S + dS, X, r, b, v) -
        2 * price(is_call, S, X, r, b, v) +
        price(is_call, S - dS, X, r, b, v)
    ) / dS ** 2

def ddelta_dvol(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return 1 / (4 * dS * 0.01) * (
        price(is_call, S + dS, X, r, b, v + 0.01) -
        price(is_call, S + dS, X, r, b, v - 0.01) -
        price(is_call, S - dS, X, r, b, v + 0.01) +
        price(is_call, S - dS, X, r, b, v - 0.01)
    ) / 100

def vega(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X, r, b, v + 0.01) -
        price(is_call, S, X, r, b, v - 0.01)
    ) / 2

def vomma(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X, r, b, v + 0.01) -
        2 * price(is_call, S, X, r, b, v) +
        price(is_call, S, X, r, b, v - 0.01)
    ) / 0.01 ** 2 / 10000

def vegap(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return v / 0.1 * (
            price(is_call, S, X, r, b, v + 0.01) -
            price(is_call, S, X, r, b, v - 0.01)
    ) / 2

def dvega_dvol(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X, r, b, v + 0.01) -
        2 * price(is_call, S, X, r, b, v) +
        price(is_call, S, X, r, b, v - 0.01)
    )

def rho(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X, r + 0.01, b + 0.01, v) -
        price(is_call, S, X, r - 0.01, b - 0.01, v)
    ) / 2

def futures_rho(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X, r + 0.01, b, v) -
        price(is_call, S, X, r - 0.01, b, v)
    ) / 2

def rho2(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X, r, b - 0.01, v) -
        price(is_call, S, X, r, b + 0.01, v)
    ) / 2

def carry(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X, r, b + 0.01, v) -
        price(is_call, S, X, r, b - 0.01, v)
    ) / 2

def speed(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return 1 / dS ** 3 * (
        price(is_call, S + 2 * dS, X, r, b, v) -
        3 * price(is_call, S + dS, X, r, b, v) +
        3 * price(is_call, S, X, r, b, v) -
        price(is_call, S - dS, X, r, b, v)
    )

def strike_delta(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X + dS, r, b, v) -
        price(is_call, S, X - dS, r, b, v)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.01
) -> float:
    return (
        price(is_call, S, X + dS, r, b, v) -
        2 * price(is_call, S, X, r, b, v) +
        price(is_call, S, X - dS, r, b, v)
    ) / dS ** 2
