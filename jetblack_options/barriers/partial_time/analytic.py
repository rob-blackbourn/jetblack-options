"""Barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, Callable

from ...distributions import CND, CBND
from ...european.black_scholes.analytic import price as bs_price

# Partial-time single asset barrier options
def price(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float
) -> float:
    
    d1 = (log(S / K) + (b + v ** 2 / 2) * T2) / (v * sqrt(T2))
    d2 = d1 - v * sqrt(T2)
    f1 = (log(S / K) + 2 * log(H / S) + (b + v ** 2 / 2) * T2) / (v * sqrt(T2))
    f2 = f1 - v * sqrt(T2)
    e1 = (log(S / H) + (b + v ** 2 / 2) * t1) / (v * sqrt(t1))
    e2 = e1 - v * sqrt(t1)
    e3 = e1 + 2 * log(H / S) / (v * sqrt(t1))
    e4 = e3 - v * sqrt(t1)
    mu = (b - v ** 2 / 2) / v ** 2
    rho = sqrt(t1 / T2)
    g1 = (log(S / H) + (b + v ** 2 / 2) * T2) / (v * sqrt(T2))
    g2 = g1 - v * sqrt(T2)
    g3 = g1 + 2 * log(H / S) / (v * sqrt(T2))
    g4 = g3 - v * sqrt(T2)
    
    z1 = CND(e2) - (H / S) ** (2 * mu) * CND(e4)
    z2 = CND(-e2) - (H / S) ** (2 * mu) * CND(-e4)
    z3 = CBND(g2, e2, rho) - (H / S) ** (2 * mu) * CBND(g4, -e4, -rho)
    z4 = CBND(-g2, -e2, rho) - (H / S) ** (2 * mu) * CBND(-g4, e4, -rho)
    z5 = CND(e1) - (H / S) ** (2 * (mu + 1)) * CND(e3)
    z6 = CND(-e1) - (H / S) ** (2 * (mu + 1)) * CND(-e3)
    z7 = CBND(g1, e1, rho) - (H / S) ** (2 * (mu + 1)) * CBND(g3, -e3, -rho)
    z8 = CBND(-g1, -e1, rho) - (H / S) ** (2 * (mu + 1)) * CBND(-g3, e3, -rho)
    
    if TypeFlag == "cdoA" or TypeFlag == "cuoA": # call down-and out and up-and-out type A
        if TypeFlag == "cdoA":
            eta = 1
        elif TypeFlag == "cuoA":
            eta = -1
        return S * exp((b - r) * T2) * (
            CBND(d1, eta * e1, eta * rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(f1, eta * e3, eta * rho)
        ) - K * exp(-r * T2) * (
            CBND(d2, eta * e2, eta * rho) - (H / S) ** (2 * mu) *
            CBND(f2, eta * e4, eta * rho)
        )
    elif TypeFlag == "cdoB2" and K < H: # call down-and-out type B2
        return S * exp((b - r) * T2) * (
            CBND(g1, e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(g3, -e3, -rho)
        ) - K * exp(-r * T2) * (
            CBND(g2, e2, rho) - (H / S) ** (2 * mu) *
            CBND(g4, -e4, -rho)
        )
    elif TypeFlag == "cdoB2" and K > H:
        return price("coB1", S, K, H, t1, T2, r, b, v)
    elif TypeFlag == "cuoB2" and K < H: # call up-and-out type B2
        return S * exp((b - r) * T2) * (
            CBND(-g1, -e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(-g3, e3, -rho)
        ) - K * exp(-r * T2) * (
            CBND(-g2, -e2, rho) - (H / S) ** (2 * mu) *
            CBND(-g4, e4, -rho)
        ) - S * exp((b - r) * T2) * (
            CBND(-d1, -e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(e3, -f1, -rho)
        ) + K * exp(-r * T2) * (
            CBND(-d2, -e2, rho) - (H / S) ** (2 * mu) *
            CBND(e4, -f2, -rho)
        )
    elif TypeFlag == "coB1" and K > H: # call out type B1
        return S * exp((b - r) * T2) * (
            CBND(d1, e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(f1, -e3, -rho)
        ) - K * exp(-r * T2) * (
            CBND(d2, e2, rho) - (H / S) ** (2 * mu) *
            CBND(f2, -e4, -rho)
        )
    elif TypeFlag == "coB1" and K < H:
        return S * exp((b - r) * T2) * (
            CBND(-g1, -e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(-g3, e3, -rho)
        ) - K * exp(-r * T2) * (
            CBND(-g2, -e2, rho) - (H / S) ** (2 * mu) *
            CBND(-g4, e4, -rho)
        ) - S * exp((b - r) * T2) * (
            CBND(-d1, -e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(-f1, e3, -rho)
        ) + K * exp(-r * T2) * (
            CBND(-d2, -e2, rho) - (H / S) ** (2 * mu) *
            CBND(-f2, e4, -rho)
        ) + S * exp((b - r) * T2) * (
            CBND(g1, e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(g3, -e3, -rho)
        ) - K * exp(-r * T2) * (
            CBND(g2, e2, rho) - (H / S) ** (2 * mu) *
            CBND(g4, -e4, -rho)
        )
    elif TypeFlag == "pdoA": # put down-and out and up-and-out type A
        return (
            price("cdoA", S, K, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z5 +
            K * exp(-r * T2) * z1
        )
    elif TypeFlag == "puoA":
        return (
            price("cuoA", S, K, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z6 +
            K * exp(-r * T2) * z2
        )
    elif TypeFlag == "poB1": # put out type B1
        return (
            price("coB1", S, K, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z8 +
            K * exp(-r * T2) * z4 -
            S * exp((b - r) * T2) * z7 +
            K * exp(-r * T2) * z3
        )
    elif TypeFlag == "pdoB2": # put down-and-out type B2
        return (
            price("cdoB2", S, K, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z7 +
            K * exp(-r * T2) * z3
        )
    elif TypeFlag == "puoB2": # put up-and-out type B2
        return (
            price("cuoB2", S, K, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z8 +
            K * exp(-r * T2) * z4
        )
    else:
        raise ValueError("invalid type flag")
    
def EPartialTimeBarrier(
        OutPutFlag: Literal['p', 'd', 'dddv', 'g', 'gp', 'gv', 'v', 'vp', 'dvdv', 'r', 'f', 'fr', 'b', 't', 'dx', 'dxdx'],
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        K: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.0001
    
    if (
        (TypeFlag == "cuoA" and S >= H) or
        (TypeFlag == "puoA" and S >= H) or
        (TypeFlag == "cdoA" and S <= H) or
        (TypeFlag == "pdoA" and S <= H)
    ):
        return 0
      
    if OutPutFlag == "p": # Value
        return price(TypeFlag, S, K, H, t1, T2, r, b, v)
    elif OutPutFlag == "d": # Delta
        return (price(TypeFlag, S + dS, K, H, t1, T2, r, b, v) - price(TypeFlag, S - dS, K, H, t1, T2, r, b, v)) / (2 * dS)
    elif OutPutFlag == "dddv": # DeltaDVol
        return 1 / (4 * dS * 0.01) * (price(TypeFlag, S + dS, K, H, t1, T2, r, b, v + 0.01) - price(TypeFlag, S + dS, K, H, t1, T2, r, b, v - 0.01) - price(TypeFlag, S - dS, K, H, t1, T2, r, b, v + 0.01) + price(TypeFlag, S - dS, K, H, t1, T2, r, b, v - 0.01)) / 100
    elif OutPutFlag == "g": # Gamma
        return (price(TypeFlag, S + dS, K, H, t1, T2, r, b, v) - 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v) + price(TypeFlag, S - dS, K, H, t1, T2, r, b, v)) / (dS ** 2)
    elif OutPutFlag == "gp": #  GammaP
        return S / 100 * EPartialTimeBarrier("g", TypeFlag, S + dS, K, H, t1, T2, r, b, v)
    elif OutPutFlag == "gv": # DGammaDvol
        return (price(TypeFlag, S + dS, K, H, t1, T2, r, b, v + 0.01) - 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v + 0.01) + price(TypeFlag, S - dS, K, H, t1, T2, r, b, v + 0.01) - price(TypeFlag, S + dS, K, H, t1, T2, r, b, v - 0.01) + 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v - 0.01) - price(TypeFlag, S - dS, K, H, t1, T2, r, b, v - 0.01)) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "v": #  Vega
        return (price(TypeFlag, S, K, H, t1, T2, r, b, v + 0.01) - price(TypeFlag, S, K, H, t1, T2, r, b, v - 0.01)) / 2
    elif OutPutFlag == "vp": #  VegaP
        return v / 0.1 * EPartialTimeBarrier("v", TypeFlag, S + dS, K, H, t1, T2, r, b, v)
    elif OutPutFlag == "dvdv": # DvegaDvol/vomma
        return (price(TypeFlag, S, K, H, t1, T2, r, b, v + 0.01) - 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v) + price(TypeFlag, S, K, H, t1, T2, r, b, v - 0.01)) / 0.01 ** 2 / 10000
    elif OutPutFlag == "r": # Rho
        return (price(TypeFlag, S, K, H, t1, T2, r + 0.01, b + 0.01, v) - price(TypeFlag, S, K, H, t1, T2, r - 0.01, b - 0.01, v)) / 2
    elif OutPutFlag == "fr": # Futures option Rho
        return (price(TypeFlag, S, K, H, t1, T2, r + 0.01, 0, v) - price(TypeFlag, S, K, H, t1, T2, r - 0.01, 0, v)) / 2
    elif OutPutFlag == "f": # Rho2 Phi
        return (price(TypeFlag, S, K, H, t1, T2, r, b - 0.01, v) - price(TypeFlag, S, K, H, t1, T2, r, b + 0.01, v)) / 2
    elif OutPutFlag == "b": #  Carry sensitivity
        return (price(TypeFlag, S, K, H, t1, T2, r, b + 0.01, v) - price(TypeFlag, S, K, H, t1, T2, r, b - 0.01, v)) / 2
    elif OutPutFlag == "t": # Theta
        if t1 <= 1 / 365:
            return (price(TypeFlag, S, K, H, 0.00001, T2 - 1 / 365, r, b, v) - price(TypeFlag, S, K, H, t1, T2, r, b, v))
        else:
            return (price(TypeFlag, S, K, H, t1 - 1 / 365, T2 - 1 / 365, r, b, v) - price(TypeFlag, S, K, H, t1, T2, r, b, v))
    elif OutPutFlag == "dx": # Strike Delta
        return (price(TypeFlag, S, K + dS, H, t1, T2, r, b, v) - price(TypeFlag, S, K - dS, H, t1, T2, r, b, v)) / (2 * dS)
    elif OutPutFlag == "dxdx": # Strike Gamma
        return (price(TypeFlag, S, K + dS, H, t1, T2, r, b, v) - 2 * price(TypeFlag, S, K, H, t1, T2, r, b, v) + price(TypeFlag, S, K - dS, H, t1, T2, r, b, v)) / (dS ** 2)
    else:
        raise ValueError('unhandled output flag')
