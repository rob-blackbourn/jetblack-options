"""Plain Vanilla"""

from math import exp, log, pi, sqrt
from typing import Literal, Optional

from ..distributions import CND, ND, CNDEV, CHIINV

# The generalized Black and Scholes formula
def GBlackScholes(
        CallPutFlag: Literal['c', 'p'],
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
)-> float:

    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if CallPutFlag == "c":
        return S * exp((b - r) * T) * CND(d1) - X * exp(-r * T) * CND(d2)
    elif CallPutFlag == "p":
        return X * exp(-r * T) * CND(-d2) - S * exp((b - r) * T) * CND(-d1)
    else:
        raise ValueError(f'Invalid CallPutFlag "{CallPutFlag}"')
    
def GBlackScholesNGreeks(
        OutPutFlag: Literal['p', 'd', 'dP', 'e', 'g', 'gv', 'gp', 'dddv', 'v', 'vp', 'dvdv', 't', 'r', 'fr', 'f', 'b', 's', 'dx', 'dxdx'],
        CallPutFlag: Literal['c', 'p'],
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.01

    if OutPutFlag == "p": # Value
        return GBlackScholes(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "d": # Delta
        return (
            GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v) -
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dP": # Delta %
        return (
            GBlackScholes(CallPutFlag, S * (1 + dS), X, T, r, b, v) -
            GBlackScholes(CallPutFlag, S * (1 - dS), X, T, r, b, v)
        ) * 2 / S
    elif OutPutFlag == "e": # Elasticity
        return (
            GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v) -
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v)
        ) / (2 * dS) * S / GBlackScholes(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "g": # Gamma
        return (
            GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * GBlackScholes(CallPutFlag, S, X, T, r, b, v) +
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": # DGammaDVol
        return (
            GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            2 * GBlackScholes(CallPutFlag, S, X, T, r, b, v + 0.01) +
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v + 0.01) -
            GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v - 0.01) +
            2 * GBlackScholes(CallPutFlag, S, X, T, r, b, v - 0.01) -
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": # GammaP
        return S / 100 * (
            GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * GBlackScholes(CallPutFlag, S, X, T, r, b, v) +
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "dddv": # DDeltaDvol
        return 1 / (4 * dS * 0.01) * (
            GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v - 0.01) -
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v + 0.01) +
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "v": # Vega
        return (
            GBlackScholes(CallPutFlag, S, X, T, r, b, v + 0.01) -
            GBlackScholes(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "vp": # VegaP
        return v / 0.1 * (
            GBlackScholes(CallPutFlag, S, X, T, r, b, v + 0.01) -
            GBlackScholes(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": # DvegaDvol
        return (
            GBlackScholes(CallPutFlag, S, X, T, r, b, v + 0.01) -
            2 * GBlackScholes(CallPutFlag, S, X, T, r, b, v) +
            GBlackScholes(CallPutFlag, S, X, T, r, b, v - 0.01)
        )
    elif OutPutFlag == "t": # Theta
        if T <= 1 / 365:
            return (
                GBlackScholes(CallPutFlag, S, X, 0.00001, r, b, v) -
                GBlackScholes(CallPutFlag, S, X, T, r, b, v)
            )
        else:
            return (
                GBlackScholes(CallPutFlag, S, X, T - 1 / 365, r, b, v) -
                GBlackScholes(CallPutFlag, S, X, T, r, b, v)
            )
    elif OutPutFlag == "r": # Rho
        return (
            GBlackScholes(CallPutFlag, S, X, T, r + 0.01, b + 0.01, v) -
            GBlackScholes(CallPutFlag, S, X, T, r - 0.01, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "fr": # Futures options rho
        return (
            GBlackScholes(CallPutFlag, S, X, T, r + 0.01, 0, v) -
            GBlackScholes(CallPutFlag, S, X, T, r - 0.01, 0, v)
        ) / 2
    elif OutPutFlag == "f": # Rho2
        return (
            GBlackScholes(CallPutFlag, S, X, T, r, b - 0.01, v) -
            GBlackScholes(CallPutFlag, S, X, T, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": # Carry
        return (
            GBlackScholes(CallPutFlag, S, X, T, r, b + 0.01, v) -
            GBlackScholes(CallPutFlag, S, X, T, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "s": # Speed
        return 1 / dS ** 3 * (
            GBlackScholes(CallPutFlag, S + 2 * dS, X, T, r, b, v) -
            3 * GBlackScholes(CallPutFlag, S + dS, X, T, r, b, v) +
            3 * GBlackScholes(CallPutFlag, S, X, T, r, b, v) -
            GBlackScholes(CallPutFlag, S - dS, X, T, r, b, v)
        )
    elif OutPutFlag == "dx": # Strike Delta
        return (
            GBlackScholes(CallPutFlag, S, X + dS, T, r, b, v) -
            GBlackScholes(CallPutFlag, S, X - dS, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": # Strike Gamma
        return (
            GBlackScholes(CallPutFlag, S, X + dS, T, r, b, v) -
            2 * GBlackScholes(CallPutFlag, S, X, T, r, b, v) +
            GBlackScholes(CallPutFlag, S, X - dS, T, r, b, v)
        ) / dS ** 2
    else:
        raise ValueError("invalid output flag")

# Black (1976) Options on futures/forwards
def Black76(CallPutFlag: Literal['c', 'p'], F: float, X: float, T: float, r: float, v: float) -> float:

    d1 = (log(F / X) + (v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if CallPutFlag == "c":
        return exp(-r * T) * (F * CND(d1) - X * CND(d2))
    elif CallPutFlag == "p":
        return exp(-r * T) * (X * CND(-d2) - F * CND(-d1))
    else:
        raise ValueError("invalid call put flag")

# Garman and Kolhagen (1983) Currency options
def GarmanKolhagen(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, rf: float, v: float) -> float:
                
    d1 = (log(S / X) + (r - rf + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if CallPutFlag == "c":
        return S * exp(-rf * T) * CND(d1) - X * exp(-r * T) * CND(d2)
    elif CallPutFlag == "p":
        return X * exp(-r * T) * CND(-d2) - S * exp(-rf * T) * CND(-d1)
    else:
        raise ValueError("invalid call put flag")


# Delta for the generalized Black and Scholes formula
def GDelta(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if CallPutFlag == "c":
        return exp((b - r) * T) * CND(d1)
    elif CallPutFlag == "p":
        return -exp((b - r) * T) * CND(-d1)
    else:
        raise ValueError("invalid call put flag")


# Forward delta for the generalized Black and Scholes formula
def GForwardDelta(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
                
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    
    if CallPutFlag == "c":
        return exp(-r * T) * CND(d1)
    elif CallPutFlag == "p":
        return exp(-r * T) * (CND(d1) - 1)
    else:
        raise ValueError("invalid call put flag")

# DDeltaDvol also known as vanna
def GDdeltaDvol(S: float, X: float, T: float, r: float, b: float, v: float) -> float:

    d1 = (log(S / X) + (b + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return -exp((b - r) * T) * d2 / v * ND(d1)

# DDeltaDvolDvol also known as DVannaDvol
def GDdeltaDvolDvol(S: float, X: float, T: float, r: float, b: float, v: float) -> float:

    d1 = (log(S / X) + (b + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return GDdeltaDvol(S, X, T, r, b, v) * 1 / v * (d1 * d2 - d1 / d2 - 1)

# DdeltaDtime/Charm for the generalized Black and Scholes formula
def GDdeltaDtime(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:

    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    
    if CallPutFlag == "c":
        return -exp((b - r) * T) * (ND(d1) * (b / (v * sqrt(T)) - d2 / (2 * T)) + (b - r) * CND(d1))
    elif CallPutFlag == "p":
        return -exp((b - r) * T) * (ND(d1) * (b / (v * sqrt(T)) - d2 / (2 * T)) - (b - r) * CND(-d1))
    else:
        raise ValueError("invalid call put flag")


# SaddleGamma for the generalized Black and Scholes formula
def GSaddleGamma(X: float, T: float, r: float, b: float, v: float)-> float:
    return sqrt(exp(1) / pi) * sqrt((2 * b - r) / v ** 2 + 1) / X


# Gamma for the generalized Black and Scholes formula
def GGamma(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return exp((b - r) * T) * ND(d1) / (S * v * sqrt(T))


# DgammaDspot/Speed for the generalized Black and Scholes formula
def GDgammaDspot(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -GGamma(S, X, T, r, b, v) * (1 + d1 / (v * sqrt(T))) / S


# DgammaDvol/Zomma for the generalized Black and Scholes formula
def GDgammaDvol(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return GGamma(S, X, T, r, b, v) * ((d1 * d2 - 1) / v)


# GGammaDtime for the generalized Black and Scholes formula
def GDgammaDtime(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return GGamma(S, X, T, r, b, v) * (r - b + b * d1 / (v * sqrt(T)) + (1 - d1 * d2) / (2 * T))


# GammaP for the generalized Black and Scholes formula
def GGammaP(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    return S * GGamma(S, X, T, r, b, v) / 100


# DgammaPDspot/SpeedP for the generalized Black and Scholes formula
def GDgammaPDspot(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -GGamma(S, X, T, r, b, v) * (d1) / (100 * v * sqrt(T))

# DgammaPDvol for the generalized Black and Scholes formula
def GDgammaPDvol(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S / 100 * GGamma(S, X, T, r, b, v) * ((d1 * d2 - 1) / v)

# GGammaPDtime for the generalized Black and Scholes formula
def GDgammaPDtime(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        GGammaP(S, X, T, r, b, v) * 
        (r - b + b * d1 / (v * sqrt(T)) + (1 - d1 * d2) / (2 * T))
    )


# Vega for the generalized Black and Scholes formula
def GVega(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * exp((b - r) * T) * ND(d1) * sqrt(T)


# DvegaDtime for the generalized Black and Scholes formula
def GDvegaDtime(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        GVega(S, X, T, r, b, v) * 
        (r - b + b * d1 / (v * sqrt(T)) - (1 + d1 * d2) / (2 * T))
    )

# DvegaDvol/Vomma for the generalized Black and Scholes formula
def GDvegaDvol(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return GVega(S, X, T, r, b, v) * d1 * d2 / v


# DVommaDVol for the generalized Black and Scholes formula
def GDvommaDvol(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        GDvegaDvol(S, X, T, r, b, v) *
        1 / v * (d1 * d2 - d1 / d2 - d2 / d1 - 1)
    )


# VegaP for the generalized Black and Scholes formula
def GVegaP(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    return v / 10 * GVega(S, X, T, r, b, v)


# DvegaPDvol/VommaP for the generalized Black and Scholes formula
def GDvegaPDvol(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return GVegaP(S, X, T, r, b, v) * d1 * d2 / v


# Vega for the generalized Black and Scholes formula
def GVegaLeverage(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    return (
        GVega(S, X, T, r, b, v) * v /
        GBlackScholes(CallPutFlag, S, X, T, r, b, v)
    )

# Variance-vega for the generalized Black and Scholes formula
def GVarianceVega(S: float, x: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / x) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * exp((b - r) * T) * ND(d1) * sqrt(T) / (2 * v)


# Variance-delta for the generalized Black and Scholes formula
def GVarianceDelta(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S * exp((b - r) * T) * ND(d1) * (-d2) / (2 * v ** 2)


# Variance-vomma for the generalized Black and Scholes formula
def GVarianceVomma(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S * exp((b - r) * T) * sqrt(T) / (4 * v ** 3) * ND(d1) * (d1 * d2 - 1)


# Variance-ultima for the generalized Black and Scholes formula
def GVarianceUltima(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        S * exp((b - r) * T) * sqrt(T) / 
        (8 * v ** 5) * ND(d1) * 
        ((d1 * d2 - 1) * 
         (d1 * d2 - 3) - 
         (d1 ** 2 + d2 ** 2))
    )


# Theta for the generalized Black and Scholes formula
def GTheta(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if CallPutFlag == "c":
        return (
            -S * exp((b - r) * T) * ND(d1) * v / (2 * sqrt(T)) - 
            (b - r) * S * exp((b - r) * T) * CND(d1) -
            r * X * exp(-r * T) * CND(d2)
        )
    elif CallPutFlag == "p":
        return (
            -S * exp((b - r) * T) * ND(d1) * v / (2 * sqrt(T)) +
            (b - r) * S * exp((b - r) * T) * CND(-d1) +
            r * X * exp(-r * T) * CND(-d2)
        )
    else:
        raise ValueError("invalid call put flag")


# Drift-less Theta for the generalized Black and Scholes formula
def GThetaDriftLess(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -S * exp((b - r) * T) * ND(d1) * v / (2 * sqrt(T))



# Rho for the generalized Black and Scholes formula for all options except futures
def GRho(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if CallPutFlag == "c":
        return T * X * exp(-r * T) * CND(d2)
    elif CallPutFlag == "p":
        return -T * X * exp(-r * T) * CND(-d2)
    else:
        raise ValueError("invalid call put flag")



# Rho for the generalized Black and Scholes formula for Futures option
def GRhoFO(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    return -T * GBlackScholes(CallPutFlag, S, X, T, r, 0, v)


# Carry rho sensitivity for the generalized Black and Scholes formula
def GCarry(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if CallPutFlag == "c":
        return T * S * exp((b - r) * T) * CND(d1)
    elif CallPutFlag == "p":
        return -T * S * exp((b - r) * T) * CND(-d1)
    else:
        raise ValueError("invalid call put flag")


# Rho2/Phi for the generalized Black and Scholes formula
def GPhi(CallPutFlag: Literal['c', 'p'], S: float, x: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / x) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if CallPutFlag == "c":
        return -T * S * exp((b - r) * T) * CND(d1)
    elif CallPutFlag == "p":
        return T * S * exp((b - r) * T) * CND(-d1)
    else:
        raise ValueError("invalid call put flag")


# DZetaDvol for the generalized Black and Scholes formula
def GDzetaDvol(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if CallPutFlag == "c":
        return -ND(d2) * d1 / v
    elif CallPutFlag == "p":
        return ND(d2) * d1 / v
    else:
        raise ValueError("invalid call put flag")



# DZetaDtime for the generalized Black and Scholes formula
def GDzetaDtime(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if CallPutFlag == "c":
        return ND(d2) * (b / (v * sqrt(T)) - d1 / (2 * T))
    elif CallPutFlag == "p":
        return -ND(d2) * (b / (v * sqrt(T)) - d1 / (2 * T))
    else:
        raise ValueError("invalid call put flag")


# Risk neutral break even probability for the generalized Black and Scholes formula
def GBreakEvenProbability(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    if CallPutFlag == "c":
        X = X + GBlackScholes("c", S, X, T, r, b, v) * exp(r * T)
        d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
        return CND(d2)
    elif CallPutFlag == "p":
        X = X - GBlackScholes("p", S, X, T, r, b, v) * exp(r * T)
        d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
        return CND(-d2)
    else:
        raise ValueError("invalid call put flag")


# StrikeDelta for the generalized Black and Scholes formula
def GStrikeDelta(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    if CallPutFlag == "c":
        return -exp(-r * T) * CND(d2)
    elif CallPutFlag == "p":
        return exp(-r * T) * CND(-d2)
    else:
        raise ValueError("invalid call put flag")


# Risk Neutral Denisty for the generalized Black and Scholes formula
def GRiskNeutralDensity(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    return exp(-r * T) * ND(d2) / (X * v * sqrt(T))


# Gamma from delta
def GGammaFromDelta(S: float, T: float, r: float, b: float, v: float, delta: float) -> float:
    return exp((b - r) * T) * ND(CNDEV(exp((r - b) * T) * abs(delta))) / (S * v * sqrt(T))


# GammaP from delta
def GGammaPFromDelta(S: float, T: float, r: float, b: float, v: float, delta: float) -> float:
    return S / 100 * GGammaFromDelta(S, T, r, b, v, delta)


# Vega from delta
def GVegaFromDelta(S: float, T: float, r: float, b: float, delta: float) -> float:
    return S * exp((b - r) * T) * sqrt(T) * ND(CNDEV(exp((r - b) * T) * abs(delta)))


# VegaP from delta
def GVegaPFromDelta(S: float, T: float, r: float, b: float, v: float, delta: float) -> float:
    return v / 10 * GVegaFromDelta(S, T, r, b, delta)


# Closed form solution to find strike given the delta
def GStrikeFromDelta(CallPutFlag: Literal['c', 'p'], S: float, T: float, r: float, b: float, v: float, delta: float) -> float:
    if CallPutFlag == "c":
        return S * exp(-CNDEV(delta * exp((r - b) * T)) * v * sqrt(T) + (b + v * v / 2) * T)
    elif CallPutFlag == "p":
        return S * exp(CNDEV(-delta * exp((r - b) * T)) * v * sqrt(T) + (b + v * v / 2) * T)
    else:
        raise ValueError("invalid call put flag")



# Closed form solution to find in-the-money risk-neutral probaility given the delta
def InTheMoneyProbFromDelta(CallPutFlag: Literal['c', 'p'], S: float, T: float, r: float, b: float, v: float, delta: float) -> float:

    if CallPutFlag == "c":
        return CND(CNDEV(delta / exp((b - r) * T)) - v * sqrt(T))
    elif CallPutFlag == "p":
        return CND(CNDEV(-delta / exp((b - r) * T)) + v * sqrt(T))
    else:
        raise ValueError("invalid call put flag")


# Closed form solution to find strike given the in-the-money risk neutral probability
def GStrikeFromInTheMoneyProb(CallPutFlag: Literal['c', 'p'], S: float, v: float, T: float, b: float, InTheMoneyProb: float) -> float:
    if CallPutFlag == "c":
        return S * exp(-CNDEV(InTheMoneyProb) * v * sqrt(T) + (b - v * v / 2) * T)
    elif CallPutFlag == "p":
        return S * exp(CNDEV(InTheMoneyProb) * v * sqrt(T) + (b - v * v / 2) * T)
    else:
        raise ValueError("invalid call put flag")


# Risk Neutral Density from in-the-money probability
def GRNDFromInTheMoneyProb(X: float, T: float, r: float, v: float, Probability: float) -> float:
    return exp(-r * T) * ND(CNDEV(Probability)) / (X * v * sqrt(T))




# Closed form solution to find in-the-money risk-neutral probaility given the delta
def GDeltaFromInTheMoneyProb(CallPutFlag: Literal['c', 'p'], S: float, T: float, r: float, b: float, v: float, InTheMoneyProb: float) -> float:
    if CallPutFlag == "c":
        return CND(CNDEV(InTheMoneyProb * exp((b - r) * T)) - v * sqrt(T))
    elif CallPutFlag == "p":
        return -CND(CNDEV(InTheMoneyProb * exp((b - r) * T)) + v * sqrt(T))
    else:
        raise ValueError("invalid call put flag")


EGBlackScholes_OutPutFlag = Literal['p', 'd', 'df', 'dddv', 'dvv', 'dt', 'dmx', 'e', 'sg', 'g', 's', 'gv', 'gt', 'gp', 'gps', 'gpv', 'gpt', 'v', 'vt', 'dvdv' 'vvv', 'vp', 'vpv', 'vl', 'varvega', 'vardelta', 'varvar', 'varult', 't', 'dlt', 'r', 'fr', 'b', 'f', 'z', 'zv', 'zt', 'bp', 'dx', 'dxdx', 'xfip', 'RNDfip', 'dfip', 'd1', 'd2', 'nd1', 'nd2', 'CNDd1', 'CNDd2']

# This is the generalized Black-Scholes-Merton formula including all greeeks
# This function is simply calling all the other functions
def EGBlackScholes(
        OutPutFlag: EGBlackScholes_OutPutFlag,
        CallPutFlag: Literal['c', 'p'],
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta: Optional[float]= None,
        InTheMoneyProb: Optional[float] = None,
        ThetaDays: Optional[float] = None
) -> float:

    output = 0
                    
    if OutPutFlag == "p": # Value
        return GBlackScholes(CallPutFlag, S, X, T, r, b, v)

    # DELTA GREEKS
    elif OutPutFlag == "d": # Delta
        return GDelta(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "df": # Forward Delta
        return GForwardDelta(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "dddv": # DDeltaDvol
        return GDdeltaDvol(S, X, T, r, b, v) / 100
    elif OutPutFlag == "dvv": # DDeltaDvolDvol
        return GDdeltaDvolDvol(S, X, T, r, b, v) / 10000
    elif OutPutFlag == "dt": # DDeltaDtime/Charm
        return GDdeltaDtime(CallPutFlag, S, X, T, r, b, v) / 365
    elif OutPutFlag == "dmx":
        return S ** 2 / X * exp((2 * b + v ** 2) * T)
    elif OutPutFlag == "e": # Elasticity
        return GElasticity(CallPutFlag, S, X, T, r, b, v)

    # GAMMA GREEKS
    elif OutPutFlag == "sg": # SaddleGamma
        return GSaddleGamma(X, T, r, b, v)
    elif OutPutFlag == "g": # Gamma
        return GGamma(S, X, T, r, b, v)
    elif OutPutFlag == "s": # DgammaDspot/speed
        return GDgammaDspot(S, X, T, r, b, v)
    elif OutPutFlag == "gv": # DgammaDvol/Zomma
        return GDgammaDvol(S, X, T, r, b, v) / 100
    elif OutPutFlag == "gt": # DgammaDtime
        return GDgammaDtime(S, X, T, r, b, v) / 365
    elif OutPutFlag == "gp": # GammaP
        return GGammaP(S, X, T, r, b, v)
    elif OutPutFlag == "gps": # DgammaPDspot
        return GDgammaPDspot(S, X, T, r, b, v)
    elif OutPutFlag == "gpv": # DgammaDvol/Zomma
        return GDgammaPDvol(S, X, T, r, b, v) / 100
    elif OutPutFlag == "gpt": # DgammaPDtime
        return GDgammaPDtime(S, X, T, r, b, v) / 365
    
    # VEGA GREEKS
    elif OutPutFlag == "v": # Vega
        return GVega(S, X, T, r, b, v) / 100
    elif OutPutFlag == "vt": # DvegaDtime
        return GDvegaDtime(S, X, T, r, b, v) / 365
    elif OutPutFlag == "dvdv": # DvegaDvol/Vomma
        return GDvegaDvol(S, X, T, r, b, v) / 10000
    elif OutPutFlag == "vvv": # DvommaDvol
        return GDvommaDvol(S, X, T, r, b, v) / 1000000
    elif OutPutFlag == "vp": # VegaP
        return GVegaP(S, X, T, r, b, v)
    elif OutPutFlag == "vpv": # DvegaPDvol/VommaP
        return GDvegaPDvol(S, X, T, r, b, v) / 100
    elif OutPutFlag == "vl": # Vega Leverage
        return GVegaLeverage(CallPutFlag, S, X, T, r, b, v)

    # VARIANCE GREEKS
    elif OutPutFlag == "varvega": # Variance-Vega
        return GVarianceVega(S, X, T, r, b, v) / 100
    elif OutPutFlag == "vardelta": # Variance-delta
        return GVarianceDelta(S, X, T, r, b, v) / 100
    elif OutPutFlag == "varvar": # Variance-vomma
        return GVarianceVomma(S, X, T, r, b, v) / 10000
    elif OutPutFlag == "varult": # Variance-ultima
        return GVarianceUltima(S, X, T, r, b, v) / 1000000

    # THETA GREEKS
    elif OutPutFlag == "t": # Theta
        return GTheta(CallPutFlag, S, X, T, r, b, v) / 365
    elif OutPutFlag == "Dlt": # Drift-less Theta
        return GThetaDriftLess(S, X, T, r, b, v) / 365

    # RATE/CARRY GREEKS
    elif OutPutFlag == "r": # Rho
        return GRho(CallPutFlag, S, X, T, r, b, v) / 100
    elif OutPutFlag == "fr": # Rho futures option
        return GRhoFO(CallPutFlag, S, X, T, r, b, v) / 100
    elif OutPutFlag == "b": # Carry Rho
        return GCarry(CallPutFlag, S, X, T, r, b, v) / 100
    elif OutPutFlag == "f": # Phi/Rho2
        return GPhi(CallPutFlag, S, X, T, r, b, v) / 100

    # PROB GREEKS
    elif OutPutFlag == "z": # Zeta/In-the-money risk neutral probability
        return GInTheMoneyProbability(CallPutFlag, S, X, T, b, v)
    elif OutPutFlag == "zv": # DzetaDvol
        return GDzetaDvol(CallPutFlag, S, X, T, r, b, v) / 100
    elif OutPutFlag == "zt": # DzetaDtime
        return GDzetaDtime(CallPutFlag, S, X, T, r, b, v) / 365
    elif OutPutFlag == "bp": # Brak even probability
        return GBreakEvenProbability(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "dx": # StrikeDelta
        return GStrikeDelta(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "dxdx": # Risk Neutral Density
        return GRiskNeutralDensity(S, X, T, r, b, v)

    # FROM DELTA GREEKS
    elif OutPutFlag == "gfd": # Gamma from delta
        assert delta is not None
        return GGammaFromDelta(S, T, r, b, v, delta)
    elif OutPutFlag == "gpfd": # GammaP from delta
        assert delta is not None
        return GGammaPFromDelta(S, T, r, b, v, delta)
    elif OutPutFlag == "vfd": # Vega from delta
        assert delta is not None
        return GVegaFromDelta(S, T, r, b, delta) / 100
    elif OutPutFlag == "vpfd": # VegaP from delta
        assert delta is not None
        return GVegaPFromDelta(S, T, r, b, v, delta)
    elif OutPutFlag == "xfd": # Strike from delta
        assert delta is not None
        return GStrikeFromDelta(CallPutFlag, S, T, r, b, v, delta)
    elif OutPutFlag == "ipfd": # In-the-money probability from delta
        assert delta is not None
        return InTheMoneyProbFromDelta(CallPutFlag, S, T, r, b, v, delta)

    # FROM IN-THE GREEKS
    elif OutPutFlag == "xfip": # Strike from in-the-money probability
        assert InTheMoneyProb is not None
        return GStrikeFromInTheMoneyProb(CallPutFlag, S, v, T, b, InTheMoneyProb)
    elif OutPutFlag == "RNDfip": # Risk Neutral Density from in-the-money probability
        assert InTheMoneyProb is not None
        return GRNDFromInTheMoneyProb(X, T, r, v, InTheMoneyProb)
    elif OutPutFlag == "dfip": # Strike from in-the-money probability
        assert InTheMoneyProb is not None
        return GDeltaFromInTheMoneyProb(CallPutFlag, S, T, r, b, v, InTheMoneyProb)

    # CALCULATIONS
    elif OutPutFlag == "d1": # d1
        return (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    elif OutPutFlag == "d2": # d2
        return (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    elif OutPutFlag == "nd1": # n(d1)
        return ND((log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T)))
    elif OutPutFlag == "nd2": # n(d2)
        return ND((log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T)))
    elif OutPutFlag == "CNDd1": # N(d1)
        return CND((log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T)))
    elif OutPutFlag == "CNDd2": # N(d2)
        return CND((log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T)))
    else:
        raise ValueError("invalid output flag")
    

# The generalized Black and Scholes formula on variance form
def GBlackScholesVariance(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v / 2) * T) / sqrt(v * T)
    d2 = d1 - sqrt(v * T)

    if CallPutFlag == "c":
        return S * exp((b - r) * T) * CND(d1) - X * exp(-r * T) * CND(d2)
    elif CallPutFlag == "p":
        return X * exp(-r * T) * CND(-d2) - S * exp((b - r) * T) * CND(-d1)
    else:
        raise ValueError('invalid call put flag')


def GBlackScholesVarianceNGreeks(
        OutPutFlag: Literal['p', 'd', 'e', 'g', 'gv', 'gp', 'dddv', 'v', 'vp', 'dvdv', 't', 'r', 'fr', 'f', 'b', 'g', 'dx', 'dxdx'],
        CallPutFlag: Literal['c', 'p'],
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.01
    
    
    if OutPutFlag == "p": # Value
        return GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "d": # Delta
        return (
             GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v) -
             GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "e": # Elasticity
         return (
             GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v) -
             GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v)
            ) / (2 * dS) * S / GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "g": # Gamma
        return (
            GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v) +
            GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": # DGammaDvariance
        return (
            GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            2 * GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v + 0.01) +
            GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v + 0.01) -
            GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v - 0.01) +
            2 * GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v - 0.01) -
            GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": # GammaP
        return S / 100 * (
            GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v) +
            GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "dddv": # DDeltaDvariance
        return 1 / (4 * dS * 0.01) * (
            GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v - 0.01) -
            GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v + 0.01) +
            GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "v": # Variance Vega
        return (
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v + 0.01) -
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "vp": # Variance VegaP
        return v / 0.1 * (
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v + 0.01) -
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": # Variance Dvegavariance
        return (
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v + 0.01) -
            2 * GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v) +
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v - 0.01)
        )
    elif OutPutFlag == "t": # Theta
        if T <= 1 / 365:
            return (
                GBlackScholesVariance(CallPutFlag, S, X, 0.00001, r, b, v) -
                GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v)
            )
        else:
            return (
                GBlackScholesVariance(CallPutFlag, S, X, T - 1 / 365, r, b, v) -
                GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v)
            )
    elif OutPutFlag == "r": # Rho
        return (
            GBlackScholesVariance(CallPutFlag, S, X, T, r + 0.01, b + 0.01, v) -
            GBlackScholesVariance(CallPutFlag, S, X, T, r - 0.01, b - 0.01, v)
        ) / (2)
    elif OutPutFlag == "fr": # Futures options rho
        return (
            GBlackScholesVariance(CallPutFlag, S, X, T, r + 0.01, 0, v) -
            GBlackScholesVariance(CallPutFlag, S, X, T, r - 0.01, 0, v)
        ) / 2
    elif OutPutFlag == "f": # Rho2
        return (
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b - 0.01, v) -
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": # Carry
        return (
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b + 0.01, v) -
            GBlackScholesVariance(CallPutFlag, S, X, T, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "s": # Speed
        return 1 / dS ** 3 * (
            GBlackScholesVariance(CallPutFlag, S + 2 * dS, X, T, r, b, v) -
            3 * GBlackScholesVariance(CallPutFlag, S + dS, X, T, r, b, v) +
            3 * GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v) -
            GBlackScholesVariance(CallPutFlag, S - dS, X, T, r, b, v)
        )
    elif OutPutFlag == "dx": # Strike Delta
        return (
            GBlackScholesVariance(CallPutFlag, S, X + dS, T, r, b, v) -
            GBlackScholesVariance(CallPutFlag, S, X - dS, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": # Gamma
        return (
            GBlackScholesVariance(CallPutFlag, S, X + dS, T, r, b, v) -
            2 * GBlackScholesVariance(CallPutFlag, S, X, T, r, b, v) +
            GBlackScholesVariance(CallPutFlag, S, X - dS, T, r, b, v)
        ) / dS ** 2
    else:
        raise ValueError("invalid output flag")



# What asset price that gives maximum DdeltaDvol
def MaxDdeltaDvolAsset(UpperLowerFlag: Literal['l', 'u'], x: float, T: float, b: float, v: float) -> float:
    # UpperLowerFlag"l" gives lower asset level that gives max DdeltaDvol
    # UpperLowerFlag"l" gives upper asset level that gives max DdeltaDvol
    
    if UpperLowerFlag == "l":
        return x * exp(-b * T - v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)
    elif UpperLowerFlag == "u":
        return x * exp(-b * T + v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)
    else:
        raise ValueError("invalid upper lower flag")

# What strike price that gives maximum DdeltaDvol
def MaxDdeltaDvolStrike(UpperLowerFlag: Literal['l', 'u'], S: float, T: float, b: float, v: float) -> float:
    
    # UpperLowerFlag"l" gives lower strike level that gives max DdeltaDvol
    # UpperLowerFlag"l" gives upper strike level that gives max DdeltaDvol

    if UpperLowerFlag == "l":
        return S * exp(b * T - v * sqrt(T) * sqrt(4 + T * v * 2) / 2)
    elif UpperLowerFlag == "u":
        return S * exp(b * T + v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)
    else:
        raise ValueError("invalid upper lower flag")

# What strike price that gives maximum gamma and vega
def GMaxGammaVegaatX(S: float, b: float, T: float, v: float) -> float:
    return S * exp((b + v * v / 2) * T)

# What asset price that gives maximum gamma
def GMaxGammaatS(x: float, b: float, T: float, v: float) -> float:
    return x * exp((-b - 3 * v * v / 2) * T)


# What asset price that gives maximum vega
def GMaxVegaatS(X: float, b: float, T: float, v: float) -> float:
    return X * exp((-b + v * v / 2) * T)

# Delta for the generalized Black and Scholes formula
def GInTheMoneyProbability(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, b: float, v: float) -> float:
    d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    
    if CallPutFlag == "c":
        return CND(d2)
    elif CallPutFlag == "p":
        return CND(-d2)
    else:
        raise ValueError("invalid call put flag")

# MirrorDeltaStrike, delta neutral straddle strike in the BSM formula
def GDeltaMirrorStrike(S: float, T: float, b: float, v: float) -> float:
    return S * exp((b + v ** 2 / 2) * T)

# MirrorProbabilityStrike, probability neutral straddle strike in the BSM formula
def GProbabilityMirrorStrike(S: float, T: float, b: float, v: float) -> float:
    return S * exp((b - v ** 2 / 2) * T)

# MirrorDeltaStrike, general delta symmmetric strike in the BSM formula
def GDeltaMirrorCallPutStrike(S: float, x: float, T: float, b: float, v: float) -> float:
    return S ** 2 / x * exp((2 * b + v ** 2) * T)

# Elasticity for the generalized Black and Scholes formula
def GElasticity(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    return (
        GDelta(CallPutFlag, S, X, T, r, b, v) * S /
        GBlackScholes(CallPutFlag, S, X, T, r, b, v)
    )

# Volatility estimate confidence interval
def GConfidenceIntervalVolatility(Alfa: float, n: int, VolatilityEstimate: float, UpperLower: Literal['u', 'l']) -> float:
    # UpperLower     ="L" gives the lower cofidence interval
    #                ="U" gives the upper cofidence interval
    # n: number of observations
    if UpperLower == "L":
        return VolatilityEstimate * sqrt((n - 1) / (CHIINV(Alfa / 2, n - 1)))
    elif UpperLower == "U":
        return VolatilityEstimate * sqrt((n - 1) / (CHIINV(1 - Alfa / 2, n - 1)))
    else:
        raise ValueError("invalid upper lower flag")


# Profit/Loss STD for the generalized Black and Scholes formula
def GProfitLossSTD(TypeFlag: Literal['a', 'p'], CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float, NHedges: int) -> float:
    
    if TypeFlag == "a": # in dollars
        return sqrt(pi / 4) * GVega(S, X, T, r, b, v) * v / sqrt(NHedges)
    elif TypeFlag == "p": # in percent
        return sqrt(pi / 4) * GVega(S, X, T, r, b, v) * v / sqrt(NHedges) / GBlackScholes(CallPutFlag, S, X, T, r, b, v)
    else:
        raise ValueError('invalid type flag')
