"""Barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, cast

from ...distributions import CND

from .analytics import price
    
def delta(
        is_call: bool,
        is_up: bool,
        is_in: bool,
        S: float,
        X: float,
        H: float,
        k: float,
        T: float,
        r: float,
        b: float,
        v: float,
        *,
        dS: float = 0.0001
) -> float:
    
    if (not is_up and not is_in and S <= H) or (is_up and not is_in and S >= H):
        if OutPutFlag == "p":
            return k
        else:
            return 0
    elif (OutInnFlag == "di" and S <= H) or (OutInnFlag == "ui" and S >= H):
        return EGBlackScholes(cast(EGBlackScholes_OutPutFlag, OutPutFlag), CallPutFlag, S, X, T, r, b, v)
      
    if OutPutFlag == "p": # Value
        return price(TypeFlag, S, X, H, k, T, r, b, v)
    elif OutPutFlag == "d": # Delta
        return (
            price(TypeFlag, S + dS, X, H, k, T, r, b, v) -
            price(TypeFlag, S - dS, X, H, k, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dddv": # DeltaDVol
        return 1 / (4 * dS * 0.01) * (
            price(TypeFlag, S + dS, X, H, k, T, r, b, v + 0.01) -
            price(TypeFlag, S + dS, X, H, k, T, r, b, v - 0.01) -
            price(TypeFlag, S - dS, X, H, k, T, r, b, v + 0.01) +
            price(TypeFlag, S - dS, X, H, k, T, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "g": # Gamma
        return (
            price(TypeFlag, S + dS, X, H, k, T, r, b, v) -
            2 * price(TypeFlag, S, X, H, k, T, r, b, v) +
            price(TypeFlag, S - dS, X, H, k, T, r, b, v)
        ) / (dS ** 2)
    elif OutPutFlag == "gp": #  GammaP
        return S / 100 * EStandardBarrier("g", TypeFlag, S + dS, X, H, k, T, r, b, v)
    elif OutPutFlag == "gv": # DGammaDvol
        return (price(TypeFlag, S + dS, X, H, k, T, r, b, v + 0.01) - 2 * price(TypeFlag, S, X, H, k, T, r, b, v + 0.01) + price(TypeFlag, S - dS, X, H, k, T, r, b, v + 0.01) - price(TypeFlag, S + dS, X, H, k, T, r, b, v - 0.01) + 2 * price(TypeFlag, S, X, H, k, T, r, b, v - 0.01) - price(TypeFlag, S - dS, X, H, k, T, r, b, v - 0.01)) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "v": #  Vega
        return (price(TypeFlag, S, X, H, k, T, r, b, v + 0.01) - price(TypeFlag, S, X, H, k, T, r, b, v - 0.01)) / 2
    elif OutPutFlag == "vp": #  VegaP
        return v / 0.1 * EStandardBarrier("v", TypeFlag, S + dS, X, H, k, T, r, b, v)
    elif OutPutFlag == "dvdv": # DvegaDvol/vomma
        return (price(TypeFlag, S, X, H, k, T, r, b, v + 0.01) - 2 * price(TypeFlag, S, X, H, k, T, r, b, v) + price(TypeFlag, S, X, H, k, T, r, b, v - 0.01)) / 0.01 ** 2 / 10000
    elif OutPutFlag == "r": # Rho
        return (price(TypeFlag, S, X, H, k, T, r + 0.01, b + 0.01, v) - price(TypeFlag, S, X, H, k, T, r - 0.01, b - 0.01, v)) / 2
    elif OutPutFlag == "fr": # Futures option Rho
        return (price(TypeFlag, S, X, H, k, T, r + 0.01, 0, v) - price(TypeFlag, S, X, H, k, T, r - 0.01, 0, v)) / 2
    elif OutPutFlag == "f": # Rho2 Phi
        return (price(TypeFlag, S, X, H, k, T, r, b - 0.01, v) - price(TypeFlag, S, X, H, k, T, r, b + 0.01, v)) / 2
    elif OutPutFlag == "b": #  Carry sensitivity
        return (price(TypeFlag, S, X, H, k, T, r, b + 0.01, v) - price(TypeFlag, S, X, H, k, T, r, b - 0.01, v)) / 2
    elif OutPutFlag == "t": # Theta
        if T <= 1 / 365:
            return (
                price(TypeFlag, S, X, H, k, 0.00001, r, b, v) -
                price(TypeFlag, S, X, H, k, T, r, b, v)
            )
        else:
            return (
                price(TypeFlag, S, X, H, k, T - 1 / 365, r, b, v) -
                price(TypeFlag, S, X, H, k, T, r, b, v)
            )
    elif OutPutFlag == "dx": # Strike Delta
        return (price(TypeFlag, S, X + dS, H, k, T, r, b, v) - price(TypeFlag, S, X - dS, H, k, T, r, b, v)) / (2 * dS)
    elif OutPutFlag == "dxdx": # Strike Gamma
        return (price(TypeFlag, S, X + dS, H, k, T, r, b, v) - 2 * price(TypeFlag, S, X, H, k, T, r, b, v) + price(TypeFlag, S, X - dS, H, k, T, r, b, v)) / (dS ** 2)
    else:
        raise ValueError("invalid output flag")
