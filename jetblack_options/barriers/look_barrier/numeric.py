"""Barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, Callable

from ...distributions import CND, CBND
from .partial_fixed_lookback import price as lookback_price

from .analytic import price

def _barrier_price(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float
) -> Optional[float]:
    if is_in:
        if (is_call and S >= H) or (not is_call and S <= H):
            return lookback_price(is_call, S, X, t1, T2, r, b, v)
    else:            
        if (is_call and S >= H) or (not is_call and S <= H):
            return 0

    return None

def delta(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S + dS, X, H, t1, T2, r, b, v) -
        price(is_call, is_in, S - dS, X, H, t1, T2, r, b, v)
    ) / (2 * dS)

def ddelta_dvol(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S + dS, X, H, t1, T2, r, b, v + dv, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S + dS, X, H, t1, T2, r, b, v - dv, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S - dS, X, H, t1, T2, r, b, v + dv, cnd=cnd, cbnd=cbnd)
        + price(is_call, is_in, S - dS, X, H, t1, T2, r, b, v - dv, cnd=cnd, cbnd=cbnd)
    ) / (4 * dS * 0.01) / 100

def gamma(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S + dS, X, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
        - 2 * price(is_call, is_in, S, X, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
        + price(is_call, is_in, S - dS, X, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / (dS ** 2)

def gammap(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return  gammap(is_call, is_in, S + dS, X, H, t1, T2, r, b, v, dS=dS, cnd=cnd, cbnd=cbnd) * S / 100

def dgamma_dvol(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S + dS, X, H, t1, T2, r, b, v + dv, cnd=cnd, cbnd=cbnd)
        - 2 * price(is_call, is_in, S, X, H, t1, T2, r, b, v + dv, cnd=cnd, cbnd=cbnd)
        + price(is_call, is_in, S - dS, X, H, t1, T2, r, b, v + dv, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S + dS, X, H, t1, T2, r, b, v - dv, cnd=cnd, cbnd=cbnd)
        + 2 * price(is_call, is_in, S, X, H, t1, T2, r, b, v - dv, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S - dS, X, H, t1, T2, r, b, v - dv, cnd=cnd, cbnd=cbnd)
    ) / (2 * dv * dS ** 2) / 100

def vega(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S, X, H, t1, T2, r, b, v + dv, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S, X, H, t1, T2, r, b, v - dv, cnd=cnd, cbnd=cbnd)
    ) / 2

def vegap(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return v / 0.1 * vega(is_call, is_in, S + dS, X, H, t1, T2, r, b, v, dv=dv, cnd=cnd, cbnd=cbnd)

def vomma(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dv: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S, X, H, t1, T2, r, b, v + dv, cnd=cnd, cbnd=cbnd)
        - 2 * price(is_call, is_in, S, X, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
        + price(is_call, is_in, S, X, H, t1, T2, r, b, v - dv, cnd=cnd, cbnd=cbnd)
    ) / 0.01 ** 2 / 10000

def rho(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dr: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S, X, H, t1, T2, r + dr, b + dr, v, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S, X, H, t1, T2, r - dr, b - dr, v, cnd=cnd, cbnd=cbnd)
    ) / 2

def rho2(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S, X, H, t1, T2, r, b - db, v, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S, X, H, t1, T2, r, b + db, v, cnd=cnd, cbnd=cbnd)
    ) / 2

def carry(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        db: float = 0.01,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S, X, H, t1, T2, r, b + db, v, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S, X, H, t1, T2, r, b - db, v, cnd=cnd, cbnd=cbnd)
    ) / 2

def theta(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dT: float = 1 / 365,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    if t1 <= dT:
        return (
            price(is_call, is_in, S, X, H, 0.00001, T2 - dT, r, b, v, cnd=cnd, cbnd=cbnd)
            - price(is_call, is_in, S, X, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
        )
    else:
        return (
            price(is_call, is_in, S, X, H, t1 - dT, T2 - dT, r, b, v, cnd=cnd, cbnd=cbnd)
            - price(is_call, is_in, S, X, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
        )

def strike_delta(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S, X + dS, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
        - price(is_call, is_in, S, X - dS, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / (2 * dS)

def strike_gamma(
        is_call: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        cnd: Callable[[float], float] = CND,
        cbnd: Callable[[float, float, float], float] = CBND
) -> float:

    p = _barrier_price(is_call, is_in, S, X, H, t1, T2, r, b, v)
    if p is not None:
        return p

    return (
        price(is_call, is_in, S, X + dS, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
        - 2 * price(is_call, is_in, S, X, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
        + price(is_call, is_in, S, X - dS, H, t1, T2, r, b, v, cnd=cnd, cbnd=cbnd)
    ) / (dS ** 2)
