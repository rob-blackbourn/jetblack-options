"""American"""

from math import exp, log, sqrt
from typing import Literal, Optional

from .black_scholes_merton import ImpliedVolGBlackScholes
from .distributions import CND, CBND, ND
from .european.plain_vanilla import GBlackScholes

def phi(S: float, T: float, gamma: float, h: float, i: float, r: float, b: float, v: float) -> float:
    lambda_ = (-r + gamma * b + 0.5 * gamma * (gamma - 1) * v ** 2) * T
    d = -(log(S / h) + (b + (gamma - 0.5) * v ** 2) * T) / (v * sqrt(T))
    kappa = 2 * b / v ** 2 + 2 * gamma - 1
    return exp(lambda_) * S ** gamma * (CND(d) - (i / S) ** kappa * CND(d - 2 * log(i / S) / (v * sqrt(T))))

def ksi(S: float, T2: float, gamma: float, h: float, I2: float, I1: float, t1: float, r: float, b: float, v: float) -> float:
    e1 = (log(S / I1) + (b + (gamma - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e2 = (log(I2 ** 2 / (S * I1)) + (b + (gamma - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e3 = (log(S / I1) - (b + (gamma - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    e4 = (log(I2 ** 2 / (S * I1)) - (b + (gamma - 0.5) * v ** 2) * t1) / (v * sqrt(t1))
    
    f1 = (log(S / h) + (b + (gamma - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f2 = (log(I2 ** 2 / (S * h)) + (b + (gamma - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f3 = (log(I1 ** 2 / (S * h)) + (b + (gamma - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    f4 = (log(S * I1 ** 2 / (h * I2 ** 2)) + (b + (gamma - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    
    rho = sqrt(t1 / T2)
    lambda_ = -r + gamma * b + 0.5 * gamma * (gamma - 1) * v ** 2
    kappa = 2 * b / (v ** 2) + (2 * gamma - 1)
    
    return (
        exp(lambda_ * T2) *
        S ** gamma * 
        (
            CBND(-e1, -f1, rho) - 
            (I2 / S) ** kappa * CBND(-e2, -f2, rho) - 
            (I1 / S) ** kappa * CBND(-e3, -f3, -rho) + 
            (I1 / I2) ** kappa * CBND(-e4, -f4, -rho)
            )
            )

def BSAmericanCallApprox2002(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    
    t1 = 1 / 2 * (sqrt(5) - 1) * T
    
    if b >= r: # Never optimal to exersice before maturity
            return GBlackScholes("c", S, X, T, r, b, v)
    else:
        Beta = (1 / 2 - b / v ** 2) + sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
        BInfinity = Beta / (Beta - 1) * X
        B0 = max(X, r / (r - b) * X)
        
        ht1 = -(b * t1 + 2 * v * sqrt(t1)) * X ** 2 / ((BInfinity - B0) * B0)
        ht2 = -(b * T + 2 * v * sqrt(T)) * X ** 2 / ((BInfinity - B0) * B0)
        I1 = B0 + (BInfinity - B0) * (1 - exp(ht1))
        I2 = B0 + (BInfinity - B0) * (1 - exp(ht2))
        alfa1 = (I1 - X) * I1 ** (-Beta)
        alfa2 = (I2 - X) * I2 ** (-Beta)
    
        if S >= I2:
            return S - X
        else:
            return (
                alfa2 * S ** Beta - alfa2 * phi(S, t1, Beta, I2, I2, r, b, v)
                + phi(S, t1, 1, I2, I2, r, b, v) - phi(S, t1, 1, I1, I2, r, b, v)
                - X * phi(S, t1, 0, I2, I2, r, b, v) + X * phi(S, t1, 0, I1, I2, r, b, v)
                + alfa1 * phi(S, t1, Beta, I1, I2, r, b, v) - alfa1 * ksi(S, T, Beta, I1, I2, I1, t1, r, b, v)
                + ksi(S, T, 1, I1, I2, I1, t1, r, b, v) - ksi(S, T, 1, X, I2, I1, t1, r, b, v)
                - X * ksi(S, T, 0, I1, I2, I1, t1, r, b, v) + X * ksi(S, T, 0, X, I2, I1, t1, r, b, v))


# The Bjerksund and Stensland (2002) American approximation
def BSAmericanApprox2002(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    if CallPutFlag == "c":
        return BSAmericanCallApprox2002(S, X, T, r, b, v)
    elif CallPutFlag == "p": # Use the Bjerksund and Stensland put-call transformation
        return BSAmericanCallApprox2002(X, S, T, r - b, -b, v)

def EBSAmericanApprox2002(
        OutPutFlag: Literal['p'],
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
        return BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "d": # Delta
        return (
             BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v) -
             BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "e": # Elasticity
        return (
            BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v) -
            BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v)
        ) / (2 * dS) * S / BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "g": # Gamma
        return (
            BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v) +
            BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": # DGammaDVol
        return (
            BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            2 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v + 0.01) +
            BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v + 0.01) -
            BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v - 0.01) +
            2 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v - 0.01) -
            BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": # GammaP
        return S / 100 * (
            BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v) +
            BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "tg": # time Gamma
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T + 1 / 365, r, b, v) -
            2 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v) +
            BSAmericanApprox2002(CallPutFlag, S, X, T - 1 / 365, r, b, v)
        ) / (1 / 365) ** 2
    elif OutPutFlag == "dddv": # DDeltaDvol
        return 1 / (4 * dS * 0.01) * (
            BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v - 0.01) -
            BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v + 0.01) +
            BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "v": # Vega
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v + 0.01) -
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "vv": # DvegaDvol/vomma
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v + 0.01) -
            2 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v) +
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": # VegaP
        return v / 0.1 * (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v + 0.01) -
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": # DvegaDvol
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v + 0.01) -
            2 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v) +
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v - 0.01)
        )
    elif OutPutFlag == "t": # Theta
        if T <= 1 / 365:
            return (
                BSAmericanApprox2002(CallPutFlag, S, X, 0.00001, r, b, v) -
                BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v)
            )
        else:
            return (
                BSAmericanApprox2002(CallPutFlag, S, X, T - 1 / 365, r, b, v) -
                BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v)
            )
    elif OutPutFlag == "r": # Rho
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r + 0.01, b + 0.01, v) -
            BSAmericanApprox2002(CallPutFlag, S, X, T, r - 0.01, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "fr": # Futures options rho
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r + 0.01, b, v) -
            BSAmericanApprox2002(CallPutFlag, S, X, T, r - 0.01, b, v)
        ) / 2
    elif OutPutFlag == "f": # Rho2
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b - 0.01, v) -
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": # Carry
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b + 0.01, v) -
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "s": # Speed
        return 1 / dS ** 3 * (
            BSAmericanApprox2002(CallPutFlag, S + 2 * dS, X, T, r, b, v) -
            3 * BSAmericanApprox2002(CallPutFlag, S + dS, X, T, r, b, v) +
            3 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v) -
            BSAmericanApprox2002(CallPutFlag, S - dS, X, T, r, b, v)
        )
    elif OutPutFlag == "dx": # Strike Delta
        return (
            BSAmericanApprox2002(CallPutFlag, S, X + dS, T, r, b, v) -
            BSAmericanApprox2002(CallPutFlag, S, X - dS, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": # Strike Gamma
        return (
            BSAmericanApprox2002(CallPutFlag, S, X + dS, T, r, b, v) -
            2 * BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v) +
            BSAmericanApprox2002(CallPutFlag, S, X - dS, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "di": # Difference in value between BS Approx and Black-Scholes Merton value
        return (
            BSAmericanApprox2002(CallPutFlag, S, X, T, r, b, v) -
            GBlackScholes(CallPutFlag, S, X, T, r, b, v)
        )
    else:
        raise ValueError('invalid output flag')
    

def BSAmericanCallApprox(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    
    if b >= r: # Never optimal to exercise before maturity
        return GBlackScholes("c", S, X, T, r, b, v)
    else:
        Beta = (1 / 2 - b / v ** 2) + sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
        BInfinity = Beta / (Beta - 1) * X
        B0 = max(X, r / (r - b) * X)
        ht = -(b * T + 2 * v * sqrt(T)) * B0 / (BInfinity - B0)
        i = B0 + (BInfinity - B0) * (1 - exp(ht))
        Alpha = (i - X) * i ** (-Beta)
        if S >= i:
            return S - X
        else:
            return (
                Alpha * S ** Beta -
                Alpha * phi(S, T, Beta, i, i, r, b, v) +
                phi(S, T, 1, i, i, r, b, v) -
                phi(S, T, 1, X, i, r, b, v) -
                X * phi(S, T, 0, i, i, r, b, v) +
                X * phi(S, T, 0, X, i, r, b, v)
            )

# The Bjerksund and Stensland (1993) American approximation
def BSAmericanApprox(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    if CallPutFlag == "c":
        return BSAmericanCallApprox(S, X, T, r, b, v)
    elif CallPutFlag == "p": # Use the Bjerksund and Stensland put-call transformation
        return BSAmericanCallApprox(X, S, T, r - b, -b, v)
    else:
        raise ValueError('invalid call put flag')
    
def EBSAmericanApprox(
        OutPutFlag: Literal['p'],
        CallPutFlag: Literal['c', 'p'],
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float]
) -> float:
            
    if dS is None:
        dS = 0.01
    
    
    if OutPutFlag == "p": # Value
        return BSAmericanApprox(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "d": #Delta
        return (
            BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) -
            BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "e": #Elasticity
         return (
             BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) -
             BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)
            ) / (2 * dS) * S / BSAmericanApprox(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "g": #Gamma
        return (
            BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v) +
            BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": #DGammaDVol
        return (
            BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            2 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) +
            BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v + 0.01) -
            BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v - 0.01) +
            2 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01) -
            BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": #GammaP
        return S / 100 * (BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) - 2 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v) + BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)) / dS ** 2
    elif OutPutFlag == "tg": #time Gamma
        return (BSAmericanApprox(CallPutFlag, S, X, T + 1 / 365, r, b, v) - 2 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v) + BSAmericanApprox(CallPutFlag, S, X, T - 1 / 365, r, b, v)) / (1 / 365) ** 2
    elif OutPutFlag == "dddv": #DDeltaDvol
        return 1 / (4 * dS * 0.01) * (BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v + 0.01) - BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v - 0.01) - BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v + 0.01) + BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v - 0.01)) / 100
    elif OutPutFlag == "v": #Vega
         return (BSAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) - BSAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01)) / 2
    elif OutPutFlag == "vv": #DvegaDvol/vomma
        return (BSAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) - 2 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v) + BSAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01)) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": #VegaP
        return v / 0.1 * (BSAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) - BSAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01)) / 2
    elif OutPutFlag == "dvdv": #DvegaDvol
        return (BSAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) - 2 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v) + BSAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01))
    elif OutPutFlag == "t": #Theta
        if T <= 1 / 365:
            return BSAmericanApprox(CallPutFlag, S, X, 0.00001, r, b, v) - BSAmericanApprox(CallPutFlag, S, X, T, r, b, v)
        else:
            return BSAmericanApprox(CallPutFlag, S, X, T - 1 / 365, r, b, v) - BSAmericanApprox(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "r": #Rho
        return (BSAmericanApprox(CallPutFlag, S, X, T, r + 0.01, b + 0.01, v) - BSAmericanApprox(CallPutFlag, S, X, T, r - 0.01, b - 0.01, v)) / (2)
    elif OutPutFlag == "fr": #Futures options rho
        return (BSAmericanApprox(CallPutFlag, S, X, T, r + 0.01, b, v) - BSAmericanApprox(CallPutFlag, S, X, T, r - 0.01, b, v)) / (2)
    elif OutPutFlag == "f": #Rho2
        return (BSAmericanApprox(CallPutFlag, S, X, T, r, b - 0.01, v) - BSAmericanApprox(CallPutFlag, S, X, T, r, b + 0.01, v)) / (2)
    elif OutPutFlag == "b": #Carry
        return (BSAmericanApprox(CallPutFlag, S, X, T, r, b + 0.01, v) - BSAmericanApprox(CallPutFlag, S, X, T, r, b - 0.01, v)) / (2)
    elif OutPutFlag == "s": #Speed
        return 1 / dS ** 3 * (BSAmericanApprox(CallPutFlag, S + 2 * dS, X, T, r, b, v) - 3 * BSAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) + 3 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v) - BSAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v))
    elif OutPutFlag == "dx": #Strike Delta
        return (BSAmericanApprox(CallPutFlag, S, X + dS, T, r, b, v) - BSAmericanApprox(CallPutFlag, S, X - dS, T, r, b, v)) / (2 * dS)
    elif OutPutFlag == "dxdx": #Strike Gamma
        return (BSAmericanApprox(CallPutFlag, S, X + dS, T, r, b, v) - 2 * BSAmericanApprox(CallPutFlag, S, X, T, r, b, v) + BSAmericanApprox(CallPutFlag, S, X - dS, T, r, b, v)) / dS ** 2
    elif OutPutFlag == "di": #Difference in value between BS Approx and Black-Scholes Merton value
        return BSAmericanApprox(CallPutFlag, S, X, T, r, b, v) - GBlackScholes(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "BSIVol": #Equivalent Black-Scholes-Merton implied volatility
        if S >= S * exp(b * T):
            CallPutFlag = "c"
        else:
            CallPutFlag = "p"
        return ImpliedVolGBlackScholes(CallPutFlag, S, X, T, r, b, BSAmericanApprox(CallPutFlag, S, X, T, r, b, v))

def PerpetualOption(CallPutFlag: Literal['c', 'p'], S: float, X: float, r: float, b: float, v: float) -> float:

    y1 = 1 / 2 - b / v ** 2 + sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
    y2 = 1 / 2 - b / v ** 2 - sqrt((b / v ** 2 - 1 / 2) ** 2 + 2 * r / v ** 2)
    if CallPutFlag == "c":
        return X / (y1 - 1) * ((y1 - 1) / y1 * S / X) ** y1
    elif CallPutFlag == "p":
        return X / (1 - y2) * ((y2 - 1) / y2 * S / X) ** y2
    else:
        raise ValueError("invalid call put flag")

def EPerpetualOption(
        OutPutFlag: Literal['p', 'd'],
        CallPutFlag: Literal['c', 'p'],
        S: float,
        X: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.01
    
    
    if OutPutFlag == "p": # Value
        return PerpetualOption(CallPutFlag, S, X, r, b, v)
    elif OutPutFlag == "d": #Delta
        return (
             PerpetualOption(CallPutFlag, S + dS, X, r, b, v) -
             PerpetualOption(CallPutFlag, S - dS, X, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "e": #Elasticity
        return (
            PerpetualOption(CallPutFlag, S + dS, X, r, b, v) -
            PerpetualOption(CallPutFlag, S - dS, X, r, b, v)
        ) / (2 * dS) * S / PerpetualOption(CallPutFlag, S, X, r, b, v)
    elif OutPutFlag == "g": #Gamma
        return (
            PerpetualOption(CallPutFlag, S + dS, X, r, b, v) -
            2 * PerpetualOption(CallPutFlag, S, X, r, b, v) +
            PerpetualOption(CallPutFlag, S - dS, X, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": #DGammaDVol
        return (
            PerpetualOption(CallPutFlag, S + dS, X, r, b, v + 0.01) -
            2 * PerpetualOption(CallPutFlag, S, X, r, b, v + 0.01) +
            PerpetualOption(CallPutFlag, S - dS, X, r, b, v + 0.01) -
            PerpetualOption(CallPutFlag, S + dS, X, r, b, v - 0.01) +
            2 * PerpetualOption(CallPutFlag, S, X, r, b, v - 0.01) -
            PerpetualOption(CallPutFlag, S - dS, X, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": #GammaP
        return S / 100 * (
            PerpetualOption(CallPutFlag, S + dS, X, r, b, v) -
            2 * PerpetualOption(CallPutFlag, S, X, r, b, v) +
            PerpetualOption(CallPutFlag, S - dS, X, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "dddv": #DDeltaDvol
        return 1 / (4 * dS * 0.01) * (
            PerpetualOption(CallPutFlag, S + dS, X, r, b, v + 0.01) -
            PerpetualOption(CallPutFlag, S + dS, X, r, b, v - 0.01) -
            PerpetualOption(CallPutFlag, S - dS, X, r, b, v + 0.01) +
            PerpetualOption(CallPutFlag, S - dS, X, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "v": #Vega
        return (
            PerpetualOption(CallPutFlag, S, X, r, b, v + 0.01) -
            PerpetualOption(CallPutFlag, S, X, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "vv": #DvegaDvol/vomma
        return (
            PerpetualOption(CallPutFlag, S, X, r, b, v + 0.01) -
            2 * PerpetualOption(CallPutFlag, S, X, r, b, v) +
            PerpetualOption(CallPutFlag, S, X, r, b, v - 0.01)
        ) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": #VegaP
        return v / 0.1 * (
             PerpetualOption(CallPutFlag, S, X, r, b, v + 0.01) -
             PerpetualOption(CallPutFlag, S, X, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": #DvegaDvol
        return (
            PerpetualOption(CallPutFlag, S, X, r, b, v + 0.01) -
            2 * PerpetualOption(CallPutFlag, S, X, r, b, v) +
            PerpetualOption(CallPutFlag, S, X, r, b, v - 0.01)
        )
    elif OutPutFlag == "r": #Rho
        return (
            PerpetualOption(CallPutFlag, S, X, r + 0.01, b + 0.01, v) -
            PerpetualOption(CallPutFlag, S, X, r - 0.01, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "fr": #Futures options rho
        return (
            PerpetualOption(CallPutFlag, S, X, r + 0.01, b, v) -
            PerpetualOption(CallPutFlag, S, X, r - 0.01, b, v)
        ) / 2
    elif OutPutFlag == "f": #Rho2
        return (
            PerpetualOption(CallPutFlag, S, X, r, b - 0.01, v) -
            PerpetualOption(CallPutFlag, S, X, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": #Carry
        return (
            PerpetualOption(CallPutFlag, S, X, r, b + 0.01, v) -
            PerpetualOption(CallPutFlag, S, X, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "s": #Speed
        return 1 / dS ** 3 * (
            PerpetualOption(CallPutFlag, S + 2 * dS, X, r, b, v) -
            3 * PerpetualOption(CallPutFlag, S + dS, X, r, b, v) +
            3 * PerpetualOption(CallPutFlag, S, X, r, b, v) -
            PerpetualOption(CallPutFlag, S - dS, X, r, b, v)
        )
    elif OutPutFlag == "dx": #Strike Delta
        return (
            PerpetualOption(CallPutFlag, S, X + dS, r, b, v) -
            PerpetualOption(CallPutFlag, S, X - dS, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": #Strike Gamma
        return (
            PerpetualOption(CallPutFlag, S, X + dS, r, b, v) -
            2 * PerpetualOption(CallPutFlag, S, X, r, b, v) +
            PerpetualOption(CallPutFlag, S, X - dS, r, b, v)
        ) / dS ** 2
    else:
        raise ValueError("invalid output flag")


# Newton Raphson algorithm to solve for the critical commodity price for a Call
def Kc(X: float, T: float, r: float, b: float, v: float) -> float:

    # Calculation of seed value, Si
    N = 2 * b / v ** 2
    m = 2 * r / v ** 2
    q2u = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * m)) / 2
    su = X / (1 - 1 / q2u)
    h2 = -(b * T + 2 * v * sqrt(T)) * X / (su - X)
    Si = X + (su - X) * (1 - exp(h2))

    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Si / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q2 = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * k)) / 2
    LHS = Si - X
    RHS = GBlackScholes("c", Si, X, T, r, b, v) + (1 - exp((b - r) * T) * CND(d1)) * Si / Q2
    bi = exp((b - r) * T) * CND(d1) * (1 - 1 / Q2) + (1 - exp((b - r) * T) * CND(d1) / (v * sqrt(T))) / Q2
    E = 0.000001
    # Newton Raphson algorithm for finding critical price Si
    while abs(LHS - RHS) / X > E:
        Si = (X + RHS - bi * Si) / (1 - bi)
        d1 = (log(Si / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        LHS = Si - X
        RHS = GBlackScholes("c", Si, X, T, r, b, v) + (1 - exp((b - r) * T) * CND(d1)) * Si / Q2
        bi = exp((b - r) * T) * CND(d1) * (1 - 1 / Q2) + (1 - exp((b - r) * T) * ND(d1) / (v * sqrt(T))) / Q2

    return Si


# American call
def BAWAmericanCallApprox(S: float, X: float, T: float, r: float, b: float, v: float) -> float:

    if b >= r:
        return GBlackScholes("c", S, X, T, r, b, v)
    else:
        Sk = Kc(X, T, r, b, v)
        N = 2 * b / v ** 2
        k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
        d1 = (log(Sk / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        Q2 = (-(N - 1) + sqrt((N - 1) ** 2 + 4 * k)) / 2
        a2 = (Sk / Q2) * (1 - exp((b - r) * T) * CND(d1))
        if S < Sk:
            return GBlackScholes("c", S, X, T, r, b, v) + a2 * (S / Sk) ** Q2
        else:
            return S - X



# Newton Raphson algorithm to solve for the critical commodity price for a Put
def Kp(X: float, T: float, r: float, b: float, v: float) -> float:
    
    # Calculation of seed value, Si
    N = 2 * b / v ** 2
    m = 2 * r / v ** 2
    q1u = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * m)) / 2
    su = X / (1 - 1 / q1u)
    h1 = (b * T - 2 * v * sqrt(T)) * X / (X - su)
    Si = su + (X - su) * exp(h1)

    k = 2 * r / (v * 2 * (1 - exp(-r * T)))
    d1 = (log(Si / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q1 = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * k)) / 2
    LHS = X - Si
    RHS = GBlackScholes("p", Si, X, T, r, b, v) - (1 - exp((b - r) * T) * CND(-d1)) * Si / Q1
    bi = -exp((b - r) * T) * CND(-d1) * (1 - 1 / Q1) - (1 + exp((b - r) * T) * ND(-d1) / (v * sqrt(T))) / Q1
    E = 0.000001
    # Newton Raphson algorithm for finding critical price Si
    while abs(LHS - RHS) / X > E:
        Si = (X - RHS + bi * Si) / (1 + bi)
        d1 = (log(Si / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
        LHS = X - Si
        RHS = GBlackScholes("p", Si, X, T, r, b, v) - (1 - exp((b - r) * T) * CND(-d1)) * Si / Q1
        bi = -exp((b - r) * T) * CND(-d1) * (1 - 1 / Q1) - (1 + exp((b - r) * T) * CND(-d1) / (v * sqrt(T))) / Q1

    return Si

# American put
def BAWAmericanPutApprox(S: float, X: float, T: float, r: float, b: float, v: float) -> float:

    Sk = Kp(X, T, r, b, v)
    N = 2 * b / v ** 2
    k = 2 * r / (v ** 2 * (1 - exp(-r * T)))
    d1 = (log(Sk / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    Q1 = (-(N - 1) - sqrt((N - 1) ** 2 + 4 * k)) / 2
    a1 = -(Sk / Q1) * (1 - exp((b - r) * T) * CND(-d1))

    if S > Sk:
        return GBlackScholes("p", S, X, T, r, b, v) + a1 * (S / Sk) ** Q1
    else:
        return X - S

# The Barone-Adesi and Whaley (1987) American approximation
def BAWAmericanApprox(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    if CallPutFlag == "c":
        return BAWAmericanCallApprox(S, X, T, r, b, v)
    elif CallPutFlag == "p":
        return BAWAmericanPutApprox(S, X, T, r, b, v)
    else:
        raise ValueError("invalid call put flag")

def EBAWAmericanApprox(
        OutPutFlag: Literal['p', 'd'],
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
        return BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "d": #Delta
        return (
             BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) -
             BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "e": #Elasticity
        return (
            BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) -
            BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)
        ) / (2 * dS) * S / BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v)
    elif OutPutFlag == "g": #Gamma
        return (
            BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v) +
            BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": #DGammaDVol
        return (
            BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            2 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) +
            BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v + 0.01) -
            BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v - 0.01) +
            2 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01) -
            BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": #GammaP
        return S / 100 * (
            BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) -
            2 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v) +
            BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "tg": #time Gamma
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T + 1 / 365, r, b, v) -
            2 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v) +
            BAWAmericanApprox(CallPutFlag, S, X, T - 1 / 365, r, b, v)
        ) / (1 / 365) ** 2
    elif OutPutFlag == "dddv": #DDeltaDvol
        return 1 / (4 * dS * 0.01) * (
            BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v + 0.01) -
            BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v - 0.01) -
            BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v + 0.01) +
            BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "v": #Vega
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) -
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "vv": #DvegaDvol/vomma
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) -
            2 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v) +
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": #VegaP
        return v / 0.1 * (
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) -
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": #DvegaDvol
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v + 0.01) -
            2 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v) +
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v - 0.01)
        )
    elif OutPutFlag == "t": #Theta
        if T <= 1 / 365:
            return (
                BAWAmericanApprox(CallPutFlag, S, X, 0.00001, r, b, v) -
                BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v)
            )
        else:
            return (
                BAWAmericanApprox(CallPutFlag, S, X, T - 1 / 365, r, b, v) -
                BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v)
            )
    elif OutPutFlag == "r": #Rho
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T, r + 0.01, b + 0.01, v) -
            BAWAmericanApprox(CallPutFlag, S, X, T, r - 0.01, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "fr": #Futures options rho
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T, r + 0.01, b, v) -
            BAWAmericanApprox(CallPutFlag, S, X, T, r - 0.01, b, v)
        ) / 2
    elif OutPutFlag == "f": #Rho2
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b - 0.01, v) -
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": #Carry
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b + 0.01, v) -
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "s": #Speed
        return 1 / dS ** 3 * (
            BAWAmericanApprox(CallPutFlag, S + 2 * dS, X, T, r, b, v) -
            3 * BAWAmericanApprox(CallPutFlag, S + dS, X, T, r, b, v) +
            3 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v) -
            BAWAmericanApprox(CallPutFlag, S - dS, X, T, r, b, v)
        )
    elif OutPutFlag == "dx": #Strike Delta
        return (
            BAWAmericanApprox(CallPutFlag, S, X + dS, T, r, b, v) -
            BAWAmericanApprox(CallPutFlag, S, X - dS, T, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": #Strike Gamma
        return (
            BAWAmericanApprox(CallPutFlag, S, X + dS, T, r, b, v) -
            2 * BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v) +
            BAWAmericanApprox(CallPutFlag, S, X - dS, T, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "di": #Difference in value between BS Approx and Black-Scholes Merton value
        return (
            BAWAmericanApprox(CallPutFlag, S, X, T, r, b, v) -
            GBlackScholes(CallPutFlag, S, X, T, r, b, v)
        )

# Muligens forskjellig fra phi i Bjerksun Stensland 1993
def phi2(S: float, T2: float, gamma: float, h: float, i: float, r: float, b: float, v: float) -> float:
    
    lambda_ = -r + gamma * b + 0.5 * gamma * (gamma - 1) * v ** 2
    kappa = 2 * b / v ** 2 + 2 * gamma - 1
    
    d = (log(S / h) + (b + (gamma - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    d2 = (log(i ** 2 / (S * h)) + (b + (gamma - 0.5) * v ** 2) * T2) / (v * sqrt(T2))
    
    return exp(lambda_ * T2) * S ** gamma * (CND(-d) - (i / S) ** kappa * CND(-d2))
