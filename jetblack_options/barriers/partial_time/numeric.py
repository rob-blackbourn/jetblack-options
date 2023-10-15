"""Barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, Callable

from ...distributions import CND, CBND

from .analytic import price


def _is_at_barrier(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        H: float
) -> bool:
    return (
        (TypeFlag == "cuoA" and S >= H) or
        (TypeFlag == "puoA" and S >= H) or
        (TypeFlag == "cdoA" and S <= H) or
        (TypeFlag == "pdoA" and S <= H)
    )

def delta(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S + dS, K, H, t1, T2, r, b, v)
        - price(TypeFlag, S - dS, K, H, t1, T2, r, b, v)
    ) / (2 * dS)

def ddelta_dvol(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return 1 / (4 * dS * 0.01) * (
        price(TypeFlag, S + dS, K, H, t1, T2, r, b, v + dv)
        - price(TypeFlag, S + dS, K, H, t1, T2, r, b, v - dv)
        - price(TypeFlag, S - dS, K, H, t1, T2, r, b, v + dv)
        + price(TypeFlag, S - dS, K, H, t1, T2, r, b, v - dv)
    ) / 100

def gamma(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: float = 0.0001
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S + dS, K, H, t1, T2, r, b, v)
        - 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v)
        + price(TypeFlag, S - dS, K, H, t1, T2, r, b, v)
    ) / (dS ** 2)

def gammap(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: float = 0.0001
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return S / 100 * gamma(TypeFlag, S + dS, K, H, t1, T2, r, b, v, dS=dS)

def dgamma_dvol(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S + dS, K, H, t1, T2, r, b, v + dv)
        - 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v + dv)
        + price(TypeFlag, S - dS, K, H, t1, T2, r, b, v + dv)
        - price(TypeFlag, S + dS, K, H, t1, T2, r, b, v - dv)
        + 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v - dv)
        - price(TypeFlag, S - dS, K, H, t1, T2, r, b, v - dv)
    ) / (2 * dv * dS ** 2) / 100

def vega(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dv: float = 0.01
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S, K, H, t1, T2, r, b, v + dv)
        - price(TypeFlag, S, K, H, t1, T2, r, b, v - dv)
    ) / 2

def vegap(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: float = 0.0001,
        dv: float = 0.01
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0

    return v / 0.1 * vega(TypeFlag, S + dS, K, H, t1, T2, r, b, v, dv=dv)

def vomma(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dv: float = 0.01
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S, K, H, t1, T2, r, b, v + dv)
        - 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v)
        + price(TypeFlag, S, K, H, t1, T2, r, b, v - dv)
    ) / dv ** 2 / 10000

def rho(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dr: float = 0.01
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S, K, H, t1, T2, r + dr, b + dr, v)
        - price(TypeFlag, S, K, H, t1, T2, r - dr, b - dr, v)
    ) / 2

def futures_rho(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: float = 0.0001
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S + dS, K, H, t1, T2, r, b, v)
        - price(TypeFlag, S - dS, K, H, t1, T2, r, b, v)
    ) / (2 * dS)

def rho2(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        db: float = 0.01
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S, K, H, t1, T2, r, b - db, v)
        - price(TypeFlag, S, K, H, t1, T2, r, b + db, v)
    ) / 2

def carry(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: float = 0.0001
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0

    return (price(TypeFlag, S, K, H, t1, T2, r, b + 0.01, v) - price(TypeFlag, S, K, H, t1, T2, r, b - 0.01, v)) / 2

def theta(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dT: float = 1 / 365
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    if t1 <= dT:
        return (
            price(TypeFlag, S, K, H, 0.00001, T2 - dT, r, b, v)
            - price(TypeFlag, S, K, H, t1, T2, r, b, v))
    else:
        return (
            price(TypeFlag, S, K, H, t1 - dT, T2 - dT, r, b, v)
            - price(TypeFlag, S, K, H, t1, T2, r, b, v)
        )

def strike_delta(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: float = 0.0001
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S, K + dS, H, t1, T2, r, b, v)
        - price(TypeFlag, S, K - dS, H, t1, T2, r, b, v)
    ) / (2 * dS)

def strike_gamma(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: float = 0.0001
) -> float:

    if _is_at_barrier(TypeFlag, S, H):
        return 0
      
    return (
        price(TypeFlag, S, K + dS, H, t1, T2, r, b, v)
        - 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v)
        + price(TypeFlag, S, K - dS, H, t1, T2, r, b, v)
    ) / (dS ** 2)
