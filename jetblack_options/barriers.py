"""Barriers"""

from math import exp, log, sqrt
from typing import Literal, Optional, cast

from .distributions import CND, CBND
from .lookback import PartialFixedLB
from .european.plain_vanilla import GBlackScholes, EGBlackScholes, EGBlackScholes_OutPutFlag

# Soft barrier options
def SoftBarrier(TypeFlag: Literal['cdi', 'cdo', 'pui', 'puo'], S: float, X: float, L: float, U: float, T: float, r: float, b: float, v: float) -> float:

    if TypeFlag == "cdi" or TypeFlag == "cdo":
        eta = 1
    else:
        eta = -1
    
    mu = (b + v ** 2 / 2) / v ** 2
    Lambda1 = exp(-1 / 2 * v ** 2 * T * (mu + 0.5) * (mu - 0.5))
    Lambda2 = exp(-1 / 2 * v ** 2 * T * (mu - 0.5) * (mu - 1.5))
    d1 = log(U ** 2 / (S * X)) / (v * sqrt(T)) + mu * v * sqrt(T)
    d2 = d1 - (mu + 0.5) * v * sqrt(T)
    d3 = log(U ** 2 / (S * X)) / (v * sqrt(T)) + (mu - 1) * v * sqrt(T)
    d4 = d3 - (mu - 0.5) * v * sqrt(T)
    e1 = log(L ** 2 / (S * X)) / (v * sqrt(T)) + mu * v * sqrt(T)
    e2 = e1 - (mu + 0.5) * v * sqrt(T)
    e3 = log(L ** 2 / (S * X)) / (v * sqrt(T)) + (mu - 1) * v * sqrt(T)
    e4 = e3 - (mu - 0.5) * v * sqrt(T)
    
    Value = (
        eta * 1 / (U - L) * (S * exp((b - r) * T) * S ** (-2 * mu) * (S * X) ** (mu + 0.5) / (2 * (mu + 0.5)) * ((U ** 2 / (S * X)) ** (mu + 0.5) * CND(eta * d1) - Lambda1 * CND(eta * d2)
    - (L ** 2 / (S * X)) ** (mu + 0.5) * CND(eta * e1) + Lambda1 * CND(eta * e2))
    - X * exp(-r * T) * S ** (-2 * (mu - 1))
    * (S * X) ** (mu - 0.5) / (2 * (mu - 0.5))
    * ((U ** 2 / (S * X)) ** (mu - 0.5) * CND(eta * d3) - Lambda2 * CND(eta * d4)
    - (L ** 2 / (S * X)) ** (mu - 0.5) * CND(eta * e3) + Lambda2 * CND(eta * e4))))
    
    if TypeFlag == "cdi" or TypeFlag == "pui":
        return Value
    elif TypeFlag == "cdo":
        return GBlackScholes("c", S, X, T, r, b, v) - Value
    elif TypeFlag == "puo":
        return GBlackScholes("p", S, X, T, r, b, v) - Value
    else:

        raise ValueError("invalid type flag")
    
def ESoftBarrier(
        OutPutFlag: Literal['p', 'd', 'dddv', 'g', 'gp', 'gv', 'v', 'dvdv', 'vp', 'r', 'f', 'b', 'dx', 'dxdx'],
        TypeFlag: Literal['cdi', 'cdo', 'pui', 'puo'],
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.0001
    
    if OutPutFlag == "p": # Value
        return SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v)
    elif OutPutFlag == "d": #  Delta
        return (
            SoftBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v) -
            SoftBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dddv": #  DeltaDVol
        return 1 / (4 * dS * 0.01) * (
            SoftBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v + 0.01) -
            SoftBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v - 0.01) -
            SoftBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v + 0.01) +
            SoftBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "g": #  Gamma
        return (
            SoftBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v) - 2 *
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v) +
            SoftBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v)
        ) / (dS ** 2)
    elif OutPutFlag == "gp": #  GammaP
        return S / 100 * ESoftBarrier("g", TypeFlag, S + dS, X, L, U, T, r, b, v)
    elif OutPutFlag == "gv": #  DGammaDVol
        return (
            SoftBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v + 0.01) -
            2 * SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v + 0.01) +
            SoftBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v + 0.01) -
            SoftBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v - 0.01) +
            2 * SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v - 0.01) -
            SoftBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "v": #  Vega
        return (
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v + 0.01) -
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": #  DVegaDVol/Vomma
        return (
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v + 0.01) -
            2 * SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v) +
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v - 0.01)
        ) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": #  VegaP
        return v / 0.1 * ESoftBarrier("v", TypeFlag, S + dS, X, L, U, T, r, b, v)
    elif OutPutFlag == "r": #  Rho
        return (
            SoftBarrier(TypeFlag, S, X, L, U, T, r + 0.01, b + 0.01, v) -
            SoftBarrier(TypeFlag, S, X, L, U, T, r - 0.01, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "f": # Rho2/Phi
        return (
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b - 0.01, v) -
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": #  Carry sensitivity
        return (
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b + 0.01, v) -
            SoftBarrier(TypeFlag, S, X, L, U, T, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "t": # Theta
        if T <= 1 / 365:
            return (
                SoftBarrier(TypeFlag, S, X, L, U, 0.00001, r, b, v) -
                SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v)
            )
        else:
            return (
                SoftBarrier(TypeFlag, S, X, L, U, T - 1 / 365, r, b, v) -
                SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v)
            )
    elif OutPutFlag == "dx": # Strike Delta
        return (
            SoftBarrier(TypeFlag, S, X + dS, L, U, T, r, b, v) -
            SoftBarrier(TypeFlag, S, X - dS, L, U, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": #  Strike Gamma
        return (
            SoftBarrier(TypeFlag, S, X + dS, L, U, T, r, b, v) -
            2 * SoftBarrier(TypeFlag, S, X, L, U, T, r, b, v) +
            SoftBarrier(TypeFlag, S, X - dS, L, U, T, r, b, v)
        ) / (dS ** 2)
    else:
        raise ValueError('unhandled output type')

# Look-barrier options
def LookBarrier(TypeFlag: Literal['cuo', 'cui', 'pdo', 'pdi'], S: float, X: float, H: float, t1: float, T2: float, r: float, b: float, v: float) -> float:

    hh = log(H / S)
    k = log(X / S)
    mu1 = b - v ** 2 / 2
    mu2 = b + v ** 2 / 2
    rho = sqrt(t1 / T2)
    
    if TypeFlag == "cuo" or TypeFlag == "cui":
        eta = 1
        m = min(hh, k)
    elif TypeFlag == "pdo" or TypeFlag == "pdi":
        eta = -1
        m = max(hh, k)
    else:
        raise ValueError("invalid type flag")
    
    g1 = (
        CND(eta * (hh - mu2 * t1) / (v * sqrt(t1))) -
        exp(2 * mu2 * hh / v ** 2) * CND(eta * (-hh - mu2 * t1) / (v * sqrt(t1)))
    ) - (
        CND(eta * (m - mu2 * t1) / (v * sqrt(t1))) -
        exp(2 * mu2 * hh / v ** 2) * CND(eta * (m - 2 * hh - mu2 * t1) / (v * sqrt(t1)))
    )
    g2 = (
        CND(eta * (hh - mu1 * t1) / (v * sqrt(t1))) -
        exp(2 * mu1 * hh / v ** 2) * CND(eta * (-hh - mu1 * t1) / (v * sqrt(t1)))
    ) - (
        CND(eta * (m - mu1 * t1) / (v * sqrt(t1))) -
        exp(2 * mu1 * hh / v ** 2) * CND(eta * (m - 2 * hh - mu1 * t1) / (v * sqrt(t1)))
    )

    part1 = (
        S *
        exp((b - r) * T2) *
        (1 + v ** 2 / (2 * b)) *
        (
            CBND(
                eta * (m - mu2 * t1) / (v * sqrt(t1)),
                eta * (-k + mu2 * T2) / (v * sqrt(T2)),
                -rho
            ) -
            exp(2 * mu2 * hh / v ** 2) *
            CBND(
                eta * (m - 2 * hh - mu2 * t1) / (v * sqrt(t1)),
                eta * (2 * hh - k + mu2 * T2) / (v * sqrt(T2)),
                -rho
            )
        )
    )
    part2 = (
        -exp(-r * T2) * X * (
            CBND(
                eta * (m - mu1 * t1) / (v * sqrt(t1)),
                eta * (-k + mu1 * T2) / (v * sqrt(T2)),
                -rho
            ) -
            exp(2 * mu1 * hh / v ** 2) *
            CBND(
                eta * (m - 2 * hh - mu1 * t1) / (v * sqrt(t1)),
                eta * (2 * hh - k + mu1 * T2) / (v * sqrt(T2)),
                -rho
            )
        )
    )
    part3 = (
        -exp(-r * T2) * v ** 2 / (2 * b) * (
            S * (S / X) ** (-2 * b / v ** 2) *
            CBND(
                eta * (m + mu1 * t1) / (v * sqrt(t1)),
                eta * (-k - mu1 * T2) / (v * sqrt(T2)),
                -rho
            ) -
            H * (H / X) ** (-2 * b / v ** 2) *
            CBND(
                eta * (m - 2 * hh + mu1 * t1) / (v * sqrt(t1)),
                eta * (2 * hh - k - mu1 * T2) / (v * sqrt(T2)),
                -rho
            )
        )
    )
    part4 = (
        S * exp((b - r) * T2) * (
            (1 + v ** 2 / (2 * b)) *
            CND(eta * mu2 * (T2 - t1) / (v * sqrt(T2 - t1))) +
            exp(-b * (T2 - t1)) * (1 - v ** 2 / (2 * b)) *
            CND(eta * (-mu1 * (T2 - t1)) / (v * sqrt(T2 - t1)))
        ) * g1 - exp(-r * T2) * X * g2)
    OutValue = eta * (part1 + part2 + part3 + part4)

    if TypeFlag == "cuo" or TypeFlag == "pdo":
        return OutValue
    elif TypeFlag == "cui":
        return PartialFixedLB("c", S, X, t1, T2, r, b, v) - OutValue
    elif TypeFlag == "pdi":
        return PartialFixedLB("p", S, X, t1, T2, r, b, v) - OutValue
    else:
        raise ValueError("invalid type flag")
    
def ELookBarrier(
        OutPutFlag: Literal['p', 'd', 'dddv', 'g', 'gp', 'gv', 'v', 'vp'],
        TypeFlag: Literal['cuo', 'cui', 'pdo', 'pdi'],
        S: float,
        X: float,
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
     
    CallPutFlag = cast(Literal['c', 'p'], TypeFlag[0])
    
    if (TypeFlag == "cuo" and S >= H) or (TypeFlag == "pdo" and S <= H):
        return 0
    elif (TypeFlag == "cui" and S >= H) or (TypeFlag == "pdi" and S <= H):
        return PartialFixedLB(CallPutFlag, S, X, t1, T2, r, b, v)
      
    if OutPutFlag == "p": # Value
        return LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v)
    elif OutPutFlag == "d": # Delta
        return (
            LookBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v) -
            LookBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dddv": # DeltaDVol
        return 1 / (4 * dS * 0.01) * (
            LookBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v + 0.01) -
            LookBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v - 0.01) -
            LookBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v + 0.01) +
            LookBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "g": # Gamma
        return (
            LookBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v) -
            2 * LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v) +
            LookBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v)
        ) / (dS ** 2)
    elif OutPutFlag == "gp": #  GammaP
        return S / 100 * ELookBarrier("g", TypeFlag, S + dS, X, H, t1, T2, r, b, v)
    elif OutPutFlag == "gv": # DGammaDvol
        return (
            LookBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v + 0.01) -
            2 * LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v + 0.01) +
            LookBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v + 0.01) -
            LookBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v - 0.01) +
            2 * LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v - 0.01) -
            LookBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "v": #  Vega
        return (
            LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v + 0.01) -
            LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "vp": #  VegaP
        return v / 0.1 * ELookBarrier("v", TypeFlag, S + dS, X, H, t1, T2, r, b, v)
    elif OutPutFlag == "dvdv": # DvegaDvol/vomma
        return (LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v + 0.01) - 2 * LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v) + LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v - 0.01)) / 0.01 ** 2 / 10000
    elif OutPutFlag == "r": # Rho
        return (LookBarrier(TypeFlag, S, X, H, t1, T2, r + 0.01, b + 0.01, v) - LookBarrier(TypeFlag, S, X, H, t1, T2, r - 0.01, b - 0.01, v)) / 2
    elif OutPutFlag == "f": # Rho2 Phi
        return (LookBarrier(TypeFlag, S, X, H, t1, T2, r, b - 0.01, v) - LookBarrier(TypeFlag, S, X, H, t1, T2, r, b + 0.01, v)) / 2
    elif OutPutFlag == "b": #  Carry sensitivity
        return (LookBarrier(TypeFlag, S, X, H, t1, T2, r, b + 0.01, v) - LookBarrier(TypeFlag, S, X, H, t1, T2, r, b - 0.01, v)) / 2
    elif OutPutFlag == "t": # Theta
        if t1 <= 1 / 365:
            return (
                LookBarrier(TypeFlag, S, X, H, 0.00001, T2 - 1 / 365, r, b, v) -
                LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v)
            )
        else:
            return (
                LookBarrier(TypeFlag, S, X, H, t1 - 1 / 365, T2 - 1 / 365, r, b, v) -
                LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v)
            )
    elif OutPutFlag == "dx": # Strike Delta
        return (
            LookBarrier(TypeFlag, S, X + dS, H, t1, T2, r, b, v) -
            LookBarrier(TypeFlag, S, X - dS, H, t1, T2, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": # Strike Gamma
        return (
            LookBarrier(TypeFlag, S, X + dS, H, t1, T2, r, b, v) -
            2 * LookBarrier(TypeFlag, S, X, H, t1, T2, r, b, v) +
            LookBarrier(TypeFlag, S, X - dS, H, t1, T2, r, b, v)
        ) / (dS ** 2)
    else:
        raise ValueError('unhandled output flag')


# Partial-time singel asset barrier options
def PartialTimeBarrier(
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        X: float,
        H: float,
        t1: float,
        T2: float,
        r: float,
        b: float,
        v: float
) -> float:
    
    d1 = (log(S / X) + (b + v ** 2 / 2) * T2) / (v * sqrt(T2))
    d2 = d1 - v * sqrt(T2)
    f1 = (log(S / X) + 2 * log(H / S) + (b + v ** 2 / 2) * T2) / (v * sqrt(T2))
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
        ) - X * exp(-r * T2) * (
            CBND(d2, eta * e2, eta * rho) - (H / S) ** (2 * mu) *
            CBND(f2, eta * e4, eta * rho)
        )
    elif TypeFlag == "cdoB2" and X < H: # call down-and-out type B2
        return S * exp((b - r) * T2) * (
            CBND(g1, e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(g3, -e3, -rho)
        ) - X * exp(-r * T2) * (
            CBND(g2, e2, rho) - (H / S) ** (2 * mu) *
            CBND(g4, -e4, -rho)
        )
    elif TypeFlag == "cdoB2" and X > H:
        return PartialTimeBarrier("coB1", S, X, H, t1, T2, r, b, v)
    elif TypeFlag == "cuoB2" and X < H: # call up-and-out type B2
        return S * exp((b - r) * T2) * (
            CBND(-g1, -e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(-g3, e3, -rho)
        ) - X * exp(-r * T2) * (
            CBND(-g2, -e2, rho) - (H / S) ** (2 * mu) *
            CBND(-g4, e4, -rho)
        ) - S * exp((b - r) * T2) * (
            CBND(-d1, -e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(e3, -f1, -rho)
        ) + X * exp(-r * T2) * (
            CBND(-d2, -e2, rho) - (H / S) ** (2 * mu) *
            CBND(e4, -f2, -rho)
        )
    elif TypeFlag == "coB1" and X > H: # call out type B1
        return S * exp((b - r) * T2) * (
            CBND(d1, e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(f1, -e3, -rho)
        ) - X * exp(-r * T2) * (
            CBND(d2, e2, rho) - (H / S) ** (2 * mu) *
            CBND(f2, -e4, -rho)
        )
    elif TypeFlag == "coB1" and X < H:
        return S * exp((b - r) * T2) * (
            CBND(-g1, -e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(-g3, e3, -rho)
        ) - X * exp(-r * T2) * (
            CBND(-g2, -e2, rho) - (H / S) ** (2 * mu) *
            CBND(-g4, e4, -rho)
        ) - S * exp((b - r) * T2) * (
            CBND(-d1, -e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(-f1, e3, -rho)
        ) + X * exp(-r * T2) * (
            CBND(-d2, -e2, rho) - (H / S) ** (2 * mu) *
            CBND(-f2, e4, -rho)
        ) + S * exp((b - r) * T2) * (
            CBND(g1, e1, rho) - (H / S) ** (2 * (mu + 1)) *
            CBND(g3, -e3, -rho)
        ) - X * exp(-r * T2) * (
            CBND(g2, e2, rho) - (H / S) ** (2 * mu) *
            CBND(g4, -e4, -rho)
        )
    elif TypeFlag == "pdoA": # put down-and out and up-and-out type A
        return (
            PartialTimeBarrier("cdoA", S, X, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z5 +
            X * exp(-r * T2) * z1
        )
    elif TypeFlag == "puoA":
        return (
            PartialTimeBarrier("cuoA", S, X, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z6 +
            X * exp(-r * T2) * z2
        )
    elif TypeFlag == "poB1": # put out type B1
        return (
            PartialTimeBarrier("coB1", S, X, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z8 +
            X * exp(-r * T2) * z4 -
            S * exp((b - r) * T2) * z7 +
            X * exp(-r * T2) * z3
        )
    elif TypeFlag == "pdoB2": # put down-and-out type B2
        return (
            PartialTimeBarrier("cdoB2", S, X, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z7 +
            X * exp(-r * T2) * z3
        )
    elif TypeFlag == "puoB2": # put up-and-out type B2
        return (
            PartialTimeBarrier("cuoB2", S, X, H, t1, T2, r, b, v) -
            S * exp((b - r) * T2) * z8 +
            X * exp(-r * T2) * z4
        )
    else:
        raise ValueError("invalid type flag")
    
def EPartialTimeBarrier(
        OutPutFlag: Literal['p', 'd', 'dddv', 'g', 'gp', 'gv', 'v', 'vp', 'dvdv', 'r', 'f', 'fr', 'b', 't', 'dx', 'dxdx'],
        TypeFlag: Literal['cdoA', 'cuoA', 'coB1', 'coB2', 'cuoB2', 'cdoB2', 'pdoA', 'puoA', 'poB1', 'poB2'],
        S: float,
        X: float,
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
        return PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v)
    elif OutPutFlag == "d": # Delta
        return (PartialTimeBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v) - PartialTimeBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v)) / (2 * dS)
    elif OutPutFlag == "dddv": # DeltaDVol
        return 1 / (4 * dS * 0.01) * (PartialTimeBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v + 0.01) - PartialTimeBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v - 0.01) - PartialTimeBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v + 0.01) + PartialTimeBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v - 0.01)) / 100
    elif OutPutFlag == "g": # Gamma
        return (PartialTimeBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v) - 2 * PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v) + PartialTimeBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v)) / (dS ** 2)
    elif OutPutFlag == "gp": #  GammaP
        return S / 100 * EPartialTimeBarrier("g", TypeFlag, S + dS, X, H, t1, T2, r, b, v)
    elif OutPutFlag == "gv": # DGammaDvol
        return (PartialTimeBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v + 0.01) - 2 * PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v + 0.01) + PartialTimeBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v + 0.01) - PartialTimeBarrier(TypeFlag, S + dS, X, H, t1, T2, r, b, v - 0.01) + 2 * PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v - 0.01) - PartialTimeBarrier(TypeFlag, S - dS, X, H, t1, T2, r, b, v - 0.01)) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "v": #  Vega
        return (PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v + 0.01) - PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v - 0.01)) / 2
    elif OutPutFlag == "vp": #  VegaP
        return v / 0.1 * EPartialTimeBarrier("v", TypeFlag, S + dS, X, H, t1, T2, r, b, v)
    elif OutPutFlag == "dvdv": # DvegaDvol/vomma
        return (PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v + 0.01) - 2 * PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v) + PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v - 0.01)) / 0.01 ** 2 / 10000
    elif OutPutFlag == "r": # Rho
        return (PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r + 0.01, b + 0.01, v) - PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r - 0.01, b - 0.01, v)) / 2
    elif OutPutFlag == "fr": # Futures option Rho
        return (PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r + 0.01, 0, v) - PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r - 0.01, 0, v)) / 2
    elif OutPutFlag == "f": # Rho2 Phi
        return (PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b - 0.01, v) - PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b + 0.01, v)) / 2
    elif OutPutFlag == "b": #  Carry sensitivity
        return (PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b + 0.01, v) - PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b - 0.01, v)) / 2
    elif OutPutFlag == "t": # Theta
        if t1 <= 1 / 365:
            return (PartialTimeBarrier(TypeFlag, S, X, H, 0.00001, T2 - 1 / 365, r, b, v) - PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v))
        else:
            return (PartialTimeBarrier(TypeFlag, S, X, H, t1 - 1 / 365, T2 - 1 / 365, r, b, v) - PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v))
    elif OutPutFlag == "dx": # Strike Delta
        return (PartialTimeBarrier(TypeFlag, S, X + dS, H, t1, T2, r, b, v) - PartialTimeBarrier(TypeFlag, S, X - dS, H, t1, T2, r, b, v)) / (2 * dS)
    elif OutPutFlag == "dxdx": # Strike Gamma
        return (PartialTimeBarrier(TypeFlag, S, X + dS, H, t1, T2, r, b, v) - 2 * PartialTimeBarrier(TypeFlag, S, X, H, t1, T2, r, b, v) + PartialTimeBarrier(TypeFlag, S, X - dS, H, t1, T2, r, b, v)) / (dS ** 2)
    else:
        raise ValueError('unhandled output flag')


# Double barrier options
def DoubleBarrier(TypeFlag: Literal['co', 'ci', 'po', 'pi'], S: float, X: float, L: float, U: float, T: float, r: float, b: float, v: float, delta1: float, delta2: float) -> float:
    
    F = U * exp(delta1 * T)
    E = L * exp(delta2 * T)
    Sum1 = 0
    Sum2 = 0
    
    if TypeFlag == "co" or TypeFlag == "ci":
        for n in range(-5, 5+1):
            d1 = (log(S * U ** (2 * n) / (X * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d2 = (log(S * U ** (2 * n) / (F * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d3 = (log(L ** (2 * n + 2) / (X * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d4 = (log(L ** (2 * n + 2) / (F * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / v ** 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / v ** 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / v ** 2 + 1
            Sum1 = Sum1 + (U ** n / L ** n) ** mu1 * (L / S) ** mu2 * (CND(d1) - CND(d2)) - (L ** (n + 1) / (U ** n * S)) ** mu3 * (CND(d3) - CND(d4))
            Sum2 = Sum2 + (U ** n / L ** n) ** (mu1 - 2) * (L / S) ** mu2 * (CND(d1 - v * sqrt(T)) - CND(d2 - v * sqrt(T))) - (L ** (n + 1) / (U ** n * S)) ** (mu3 - 2) * (CND(d3 - v * sqrt(T)) - CND(d4 - v * sqrt(T)))
        OutValue = S * exp((b - r) * T) * Sum1 - X * exp(-r * T) * Sum2
    elif TypeFlag == "po" or TypeFlag == "pi":
        for n in range(-5, 5+1):
            d1 = (log(S * U ** (2 * n) / (E * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d2 = (log(S * U ** (2 * n) / (X * L ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d3 = (log(L ** (2 * n + 2) / (E * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            d4 = (log(L ** (2 * n + 2) / (X * S * U ** (2 * n))) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
            mu1 = 2 * (b - delta2 - n * (delta1 - delta2)) / v ** 2 + 1
            mu2 = 2 * n * (delta1 - delta2) / v ** 2
            mu3 = 2 * (b - delta2 + n * (delta1 - delta2)) / v ** 2 + 1
            Sum1 = Sum1 + (U ** n / L ** n) ** mu1 * (L / S) ** mu2 * (CND(d1) - CND(d2)) - (L ** (n + 1) / (U ** n * S)) ** mu3 * (CND(d3) - CND(d4))
            Sum2 = Sum2 + (U ** n / L ** n) ** (mu1 - 2) * (L / S) ** mu2 * (CND(d1 - v * sqrt(T)) - CND(d2 - v * sqrt(T))) - (L ** (n + 1) / (U ** n * S)) ** (mu3 - 2) * (CND(d3 - v * sqrt(T)) - CND(d4 - v * sqrt(T)))

        OutValue = X * exp(-r * T) * Sum2 - S * exp((b - r) * T) * Sum1

    if TypeFlag == "co" or TypeFlag == "po":
        return OutValue
    elif TypeFlag == "ci":
        return GBlackScholes("c", S, X, T, r, b, v) - OutValue
    elif TypeFlag == "pi":
        return GBlackScholes("p", S, X, T, r, b, v) - OutValue
    else:
        raise ValueError('invalid type flag')

def EDoubleBarrier(
        OutPutFlag: Literal['p', 'd', 'dddv', 'g', 'gp', 'gv', 'v', 'dvdv', 'vp', 'r', 'fr', 'f', 'b', 't', 'dx', 'dxdx'],
        TypeFlag: Literal['co', 'ci', 'po', 'pi'],
        S: float,
        X: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.0001
    
    OutInnFlag = cast(Literal['i', 'o'], TypeFlag[-1])
    CallPutFlag = cast(Literal['c', 'p'], TypeFlag[0])
    
    if OutInnFlag == "o" and (S <= L or S >= U):
        return 0
    elif OutInnFlag == "i" and (S <= L or S >= U):
        return EGBlackScholes(cast(EGBlackScholes_OutPutFlag, OutPutFlag), CallPutFlag, S, X, T, r, b, v)
    
    if OutPutFlag == "p": # Value
        return DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v, delta1, delta2)
    elif OutPutFlag == "d": #' Delta
        return (
            DoubleBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v, delta1, delta2) -
            DoubleBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v, delta1, delta2)
        ) / (2 * dS)
    elif OutPutFlag == "dddv": #' DeltaDVol
        return 1 / (4 * dS * 0.01) * (
            DoubleBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v + 0.01, delta1, delta2) -
            DoubleBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v - 0.01, delta1, delta2) -
            DoubleBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v + 0.01, delta1, delta2) +
            DoubleBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v - 0.01, delta1, delta2)
            ) / 100
    elif OutPutFlag == "g": #' Gamma
        return (DoubleBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v, delta1, delta2) - 2 * DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v, delta1, delta2) + DoubleBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v, delta1, delta2)) / (dS ** 2)
    elif OutPutFlag == "gp": #' GammaP
        return S / 100 * EDoubleBarrier("g", TypeFlag, S + dS, X, L, U, T, r, b, v, delta1, delta2)
    elif OutPutFlag == "gv": #' DGammaDVol
        return (
            DoubleBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v + 0.01, delta1, delta2) -
            2 * DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v + 0.01, delta1, delta2) +
            DoubleBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v + 0.01, delta1, delta2) -
            DoubleBarrier(TypeFlag, S + dS, X, L, U, T, r, b, v - 0.01, delta1, delta2) +
            2 * DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v - 0.01, delta1, delta2) -
            DoubleBarrier(TypeFlag, S - dS, X, L, U, T, r, b, v - 0.01, delta1, delta2)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "v": #' Vega
        return (
            DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v + 0.01, delta1, delta2) -
            DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v - 0.01, delta1, delta2)
        ) / 2
    elif OutPutFlag == "dvdv": #' DVegaDVol/Vomma
        return (
            DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v + 0.01, delta1, delta2) -
            2 * DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v, delta1, delta2) +
            DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v - 0.01, delta1, delta2)
        ) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": #' VegaP
        return v / 0.1 * EDoubleBarrier("v", TypeFlag, S + dS, X, L, U, T, r, b, v, delta1, delta2)
    elif OutPutFlag == "r": #' Rho
        return (
            DoubleBarrier(TypeFlag, S, X, L, U, T, r + 0.01, b + 0.01, v, delta1, delta2) -
            DoubleBarrier(TypeFlag, S, X, L, U, T, r - 0.01, b - 0.01, v, delta1, delta2)
        ) / 2
    elif OutPutFlag == "fr": #' Futures option Rho
        return (DoubleBarrier(TypeFlag, S, X, L, U, T, r + 0.01, 0, v, delta1, delta2) - DoubleBarrier(TypeFlag, S, X, L, U, T, r - 0.01, 0, v, delta1, delta2)) / 2
    elif OutPutFlag == "f": #Rho2/Phi
        return (DoubleBarrier(TypeFlag, S, X, L, U, T, r, b - 0.01, v, delta1, delta2) - DoubleBarrier(TypeFlag, S, X, L, U, T, r, b + 0.01, v, delta1, delta2)) / 2
    elif OutPutFlag == "b": #' Carry sensitivity
        return (DoubleBarrier(TypeFlag, S, X, L, U, T, r, b + 0.01, v, delta1, delta2) - DoubleBarrier(TypeFlag, S, X, L, U, T, r, b - 0.01, v, delta1, delta2)) / 2
    elif OutPutFlag == "t": #Theta
        if T <= 1 / 365:
            return DoubleBarrier(TypeFlag, S, X, L, U, 0.00001, r, b, v, delta1, delta2) - DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v, delta1, delta2)
        else:
            return DoubleBarrier(TypeFlag, S, X, L, U, T - 1 / 365, r, b, v, delta1, delta2) - DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v, delta1, delta2)
    elif OutPutFlag == "dx": #Strike Delta
        return (DoubleBarrier(TypeFlag, S, X + dS, L, U, T, r, b, v, delta1, delta2) - DoubleBarrier(TypeFlag, S, X - dS, L, U, T, r, b, v, delta1, delta2)) / (2 * dS)
    elif OutPutFlag == "dxdx": #' Strike Gamma
        return (DoubleBarrier(TypeFlag, S, X + dS, L, U, T, r, b, v, delta1, delta2) - 2 * DoubleBarrier(TypeFlag, S, X, L, U, T, r, b, v, delta1, delta2) + DoubleBarrier(TypeFlag, S, X - dS, L, U, T, r, b, v, delta1, delta2)) / (dS ** 2)
    else:
        raise ValueError("invalid output flag")
    
# Standard barrier options
def StandardBarrier(
        TypeFlag: Literal['cdi', 'cui', 'pdi', 'pui', 'cdo', 'cuo','pdo', 'puo'],
        S: float,
        X: float,
        H: float,
        k: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:

    # TypeFlag:      The "TypeFlag" gives you 8 different standard barrier options:
    #                1) "cdi"=Down-and-in call,    2) "cui"=Up-and-in call
    #                3) "pdi"=Down-and-in put,     4) "pui"=Up-and-in put
    #                5) "cdo"=Down-and-out call,   6) "cuo"=Up-out-in call
    #                7) "pdo"=Down-and-out put,    8) "puo"=Up-out-in put
    
    # Dim mu: float
    # Dim lambda_ As Double
    # Dim X1 As Double, X2 As Double
    # Dim y1 As Double, y2 As Double
    # Dim z As Double
    
    # Dim eta As Integer    'Binary variable that can take the value of 1 or -1
    # Dim phi As Integer    'Binary variable that can take the value of 1 or -1
    
    # Dim f1 As Double    'Equal to formula "A" in the book
    # Dim f2 As Double    'Equal to formula "B" in the book
    # Dim f3 As Double    'Equal to formula "C" in the book
    # Dim f4 As Double    'Equal to formula "D" in the book
    # Dim f5 As Double    'Equal to formula "E" in the book
    # Dim f6 As Double    'Equal to formula "F" in the book

    mu = (b - v ** 2 / 2) / v ** 2
    lambda_ = sqrt(mu ** 2 + 2 * r / v ** 2)
    X1 = log(S / X) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T)
    X2 = log(S / H) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T)
    y1 = log(H ** 2 / (S * X)) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T)
    y2 = log(H / S) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T)
    z = log(H / S) / (v * sqrt(T)) + lambda_ * v * sqrt(T)
    
    if TypeFlag == "cdi" or TypeFlag == "cdo":
        eta = 1
        phi = 1
    elif TypeFlag == "cui" or TypeFlag == "cuo":
        eta = -1
        phi = 1
    elif TypeFlag == "pdi" or TypeFlag == "pdo":
        eta = 1
        phi = -1
    elif TypeFlag == "pui" or TypeFlag == "puo":
        eta = -1
        phi = -1
    
    f1 = phi * S * exp((b - r) * T) * CND(phi * X1) - phi * X * exp(-r * T) * CND(phi * X1 - phi * v * sqrt(T))
    f2 = phi * S * exp((b - r) * T) * CND(phi * X2) - phi * X * exp(-r * T) * CND(phi * X2 - phi * v * sqrt(T))
    f3 = phi * S * exp((b - r) * T) * (H / S) ** (2 * (mu + 1)) * CND(eta * y1) - phi * X * exp(-r * T) * (H / S) ** (2 * mu) * CND(eta * y1 - eta * v * sqrt(T))
    f4 = phi * S * exp((b - r) * T) * (H / S) ** (2 * (mu + 1)) * CND(eta * y2) - phi * X * exp(-r * T) * (H / S) ** (2 * mu) * CND(eta * y2 - eta * v * sqrt(T))
    f5 = k * exp(-r * T) * (CND(eta * X2 - eta * v * sqrt(T)) - (H / S) ** (2 * mu) * CND(eta * y2 - eta * v * sqrt(T)))
    f6 = k * ((H / S) ** (mu + lambda_) * CND(eta * z) + (H / S) ** (mu - lambda_) * CND(eta * z - 2 * eta * lambda_ * v * sqrt(T)))
    
    if X > H:
        if TypeFlag == "cdi": # 1a) cdi
            return f3 + f5
        elif TypeFlag == "cui": # 2a) cui
            return f1 + f5
        elif TypeFlag == "pdi": # 3a) pdi
            return f2 - f3 + f4 + f5
        elif TypeFlag == "pui": # 4a) pui
            return f1 - f2 + f4 + f5
        elif TypeFlag == "cdo": # 5a) cdo
            return f1 - f3 + f6
        elif TypeFlag == "cuo": # 6a) cuo
            return f6
        elif TypeFlag == "pdo": # 7a) pdo
            return f1 - f2 + f3 - f4 + f6
        elif TypeFlag == "puo": # 8a) puo
            return f2 - f4 + f6
        else:
            raise ValueError("unhandled type flag")
    elif X < H:
        if TypeFlag == "cdi": # 1b) cdi
            return f1 - f2 + f4 + f5
        elif TypeFlag == "cui":  # 2b) cui
            return f2 - f3 + f4 + f5
        elif TypeFlag == "pdi": # 3b) pdi
            return f1 + f5
        elif TypeFlag == "pui":   # 4b) pui
            return f3 + f5
        elif TypeFlag == "cdo": # 5b) cdo
            return f2 + f6 - f4
        elif TypeFlag == "cuo": # 6b) cuo
            return f1 - f2 + f3 - f4 + f6
        elif TypeFlag == "pdo":   # 7b) pdo
            return f6
        elif TypeFlag == "puo":  # 8b) puo
            return f1 - f3 + f6
        else:
            raise ValueError("unhandled type flag")
    else:
        raise ValueError("no solution")
    
def EStandardBarrier(
        OutPutFlag: Literal['p', 'd', 'g', 'gp', 'gv', 'v', 'vp', 'dvdv', 'r', 'f', 'fr', 'b', 't'],
        TypeFlag: Literal['cdi', 'cui', 'pdi', 'pui', 'cdo', 'cuo','pdo', 'puo'],
        S: float,
        X: float,
        H: float,
        k: float,
        T: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.0001
    
    OutInnFlag = cast(Literal['i', 'o'], TypeFlag[-1])
    CallPutFlag = cast(Literal['c', 'p'], TypeFlag[0])
    
    if (OutInnFlag == "do" and S <= H) or (OutInnFlag == "uo" and S >= H):
        if OutPutFlag == "p":
            return k
        else:
            return 0
    elif (OutInnFlag == "di" and S <= H) or (OutInnFlag == "ui" and S >= H):
        return EGBlackScholes(cast(EGBlackScholes_OutPutFlag, OutPutFlag), CallPutFlag, S, X, T, r, b, v)
      
    if OutPutFlag == "p": # Value
        return StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v)
    elif OutPutFlag == "d": # Delta
        return (
            StandardBarrier(TypeFlag, S + dS, X, H, k, T, r, b, v) -
            StandardBarrier(TypeFlag, S - dS, X, H, k, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dddv": # DeltaDVol
        return 1 / (4 * dS * 0.01) * (
            StandardBarrier(TypeFlag, S + dS, X, H, k, T, r, b, v + 0.01) -
            StandardBarrier(TypeFlag, S + dS, X, H, k, T, r, b, v - 0.01) -
            StandardBarrier(TypeFlag, S - dS, X, H, k, T, r, b, v + 0.01) +
            StandardBarrier(TypeFlag, S - dS, X, H, k, T, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "g": # Gamma
        return (
            StandardBarrier(TypeFlag, S + dS, X, H, k, T, r, b, v) -
            2 * StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v) +
            StandardBarrier(TypeFlag, S - dS, X, H, k, T, r, b, v)
        ) / (dS ** 2)
    elif OutPutFlag == "gp": #  GammaP
        return S / 100 * EStandardBarrier("g", TypeFlag, S + dS, X, H, k, T, r, b, v)
    elif OutPutFlag == "gv": # DGammaDvol
        return (StandardBarrier(TypeFlag, S + dS, X, H, k, T, r, b, v + 0.01) - 2 * StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v + 0.01) + StandardBarrier(TypeFlag, S - dS, X, H, k, T, r, b, v + 0.01) - StandardBarrier(TypeFlag, S + dS, X, H, k, T, r, b, v - 0.01) + 2 * StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v - 0.01) - StandardBarrier(TypeFlag, S - dS, X, H, k, T, r, b, v - 0.01)) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "v": #  Vega
        return (StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v + 0.01) - StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v - 0.01)) / 2
    elif OutPutFlag == "vp": #  VegaP
        return v / 0.1 * EStandardBarrier("v", TypeFlag, S + dS, X, H, k, T, r, b, v)
    elif OutPutFlag == "dvdv": # DvegaDvol/vomma
        return (StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v + 0.01) - 2 * StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v) + StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v - 0.01)) / 0.01 ** 2 / 10000
    elif OutPutFlag == "r": # Rho
        return (StandardBarrier(TypeFlag, S, X, H, k, T, r + 0.01, b + 0.01, v) - StandardBarrier(TypeFlag, S, X, H, k, T, r - 0.01, b - 0.01, v)) / 2
    elif OutPutFlag == "fr": # Futures option Rho
        return (StandardBarrier(TypeFlag, S, X, H, k, T, r + 0.01, 0, v) - StandardBarrier(TypeFlag, S, X, H, k, T, r - 0.01, 0, v)) / 2
    elif OutPutFlag == "f": # Rho2 Phi
        return (StandardBarrier(TypeFlag, S, X, H, k, T, r, b - 0.01, v) - StandardBarrier(TypeFlag, S, X, H, k, T, r, b + 0.01, v)) / 2
    elif OutPutFlag == "b": #  Carry sensitivity
        return (StandardBarrier(TypeFlag, S, X, H, k, T, r, b + 0.01, v) - StandardBarrier(TypeFlag, S, X, H, k, T, r, b - 0.01, v)) / 2
    elif OutPutFlag == "t": # Theta
        if T <= 1 / 365:
            return (
                StandardBarrier(TypeFlag, S, X, H, k, 0.00001, r, b, v) -
                StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v)
            )
        else:
            return (
                StandardBarrier(TypeFlag, S, X, H, k, T - 1 / 365, r, b, v) -
                StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v)
            )
    elif OutPutFlag == "dx": # Strike Delta
        return (StandardBarrier(TypeFlag, S, X + dS, H, k, T, r, b, v) - StandardBarrier(TypeFlag, S, X - dS, H, k, T, r, b, v)) / (2 * dS)
    elif OutPutFlag == "dxdx": # Strike Gamma
        return (StandardBarrier(TypeFlag, S, X + dS, H, k, T, r, b, v) - 2 * StandardBarrier(TypeFlag, S, X, H, k, T, r, b, v) + StandardBarrier(TypeFlag, S, X - dS, H, k, T, r, b, v)) / (dS ** 2)
    else:
        raise ValueError("invalid output flag")


# Discrete barrier monitoring adjustment
def DiscreteAdjustedBarrier(S: float, H: float, v: float, dt: float) -> float:
    if H > S:
        return H * exp(0.5826 * v * sqrt(dt))
    elif H < S:
        return H * exp(-0.5826 * v * sqrt(dt))
    else:
        raise ValueError('no solution')
