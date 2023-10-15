"""Digitals"""

from math import exp, log, sqrt
from typing import Literal, Optional

from .distributions import CND

# Binary barrier options
def BinaryBarrier(
        TypeFlag: int,
        S: float,
        x: float,
        h: float,
        k: float,
        T: float,
        r: float,
        b: float,
        v: float,
        eta: int,
        phi: int
) -> float:

    # TypeFlag:  Value 1 to 28 dependent on binary option type,
    #            look in the book for specifications.
    
    mu = (b - v ** 2 / 2) / v ** 2
    lambda_ = sqrt(mu ** 2 + 2 * r / v ** 2)
    X1 = log(S / x) / (v * sqrt(T)) + (mu + 1) * v * sqrt(T)
    X2 = log(S / h) / (v * sqrt(T)) + (mu + 1) * v * sqrt(T)
    y1 = log(h ** 2 / (S * x)) / (v * sqrt(T)) + (mu + 1) * v * sqrt(T)
    y2 = log(h / S) / (v * sqrt(T)) + (mu + 1) * v * sqrt(T)
    Z = log(h / S) / (v * sqrt(T)) + lambda_ * v * sqrt(T)
    
    a1 = S * exp((b - r) * T) * CND(phi * X1)
    b1 = k * exp(-r * T) * CND(phi * X1 - phi * v * sqrt(T))
    a2 = S * exp((b - r) * T) * CND(phi * X2)
    b2 = k * exp(-r * T) * CND(phi * X2 - phi * v * sqrt(T))
    a3 = S * exp((b - r) * T) * (h / S) ** (2 * (mu + 1)) * CND(eta * y1)
    b3 = k * exp(-r * T) * (h / S) ** (2 * mu) * CND(eta * y1 - eta * v * sqrt(T))
    a4 = S * exp((b - r) * T) * (h / S) ** (2 * (mu + 1)) * CND(eta * y2)
    b4 = k * exp(-r * T) * (h / S) ** (2 * mu) * CND(eta * y2 - eta * v * sqrt(T))
    a5 = k * ((h / S) ** (mu + lambda_) * CND(eta * Z) + (h / S) ** (mu - lambda_) * CND(eta * Z - 2 * eta * lambda_ * v * sqrt(T)))
    
    if x > h:
        if TypeFlag < 5:
            return a5
        elif TypeFlag < 7:
            return b2 + b4
        elif TypeFlag < 9:
            return a2 + a4
        elif TypeFlag < 11:
            return b2 - b4
        elif TypeFlag < 13:
            return a2 - a4
        elif TypeFlag == 13:
            return b3
        elif TypeFlag == 14:
            return b3
        elif TypeFlag == 15:
            return a3
        elif TypeFlag == 16:
            return a1
        elif TypeFlag == 17:
            return b2 - b3 + b4
        elif TypeFlag == 18:
            return b1 - b2 + b4
        elif TypeFlag == 19:
            return a2 - a3 + a4
        elif TypeFlag == 20:
            return a1 - a2 + a3
        elif TypeFlag == 21:
            return b1 - b3
        elif TypeFlag == 22:
            return 0
        elif TypeFlag == 23:
            return a1 - a3
        elif TypeFlag == 24:
           return 0
        elif TypeFlag == 25:
            return b1 - b2 + b3 - b4
        elif TypeFlag == 26:
            return b2 - b4
        elif TypeFlag == 27:
            return a1 - a2 + a3 - a4
        elif TypeFlag == 28:
            return a2 - a4
        else:
            raise ValueError('unhandled type flag')
    elif x < h:
        if TypeFlag < 5:
            return a5
        elif TypeFlag < 7:
            return b2 + b4
        elif TypeFlag < 9:
            return a2 + a4
        elif TypeFlag < 11:
            return b2 - b4
        elif TypeFlag < 13:
            return a2 - a4
        elif TypeFlag == 13:
            return b1 - b2 + b4
        elif TypeFlag == 14:
            return b2 - b3 + b4
        elif TypeFlag == 15:
            return a1 - a2 + a4
        elif TypeFlag == 16:
            return a2 - a3 + a4
        elif TypeFlag == 17:
            return b1
        elif TypeFlag == 18:
            return b3
        elif TypeFlag == 19:
            return a1
        elif TypeFlag == 20:
            return a3
        elif TypeFlag == 21:
            return b2 - b4
        elif TypeFlag == 22:
            return b1 - b2 + b3 - b4
        elif TypeFlag == 23:
            return a2 - a4
        elif TypeFlag == 24:
            return a1 - a2 + a3 - a4
        elif TypeFlag == 25:
            return 0
        elif TypeFlag == 26:
            return b1 - b3
        elif TypeFlag == 27:
            return 0
        elif TypeFlag == 28:
            return a1 - a3
        else:
            raise ValueError("unhandled type flag")
    else:
        raise ValueError("no solution")
    
def EBinaryBarrier(
        OutPutFlag: Literal['p', 'd'],
        TypeFlag: int,
        S: float,
        x: float,
        h: float,
        k: float,
        T: float,
        r: float,
        b: float,
        v: float,
        eta: int,
        phi: int,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.01
   
    if S >= h and TypeFlag / 2 - int(TypeFlag / 2) == 0:
        if OutPutFlag == "p":
            if TypeFlag == 2:
                return k
            elif TypeFlag == 4:
                return S
            elif TypeFlag == 6:
                return exp(-r * T) * k
            elif TypeFlag == 8:
                return exp(-r * T) * S
            elif TypeFlag == 10:
                return 0
            elif TypeFlag == 12:
                return 0
            elif TypeFlag == 14:
                if S >= x:
                    return exp(-r * T) * k
                else:
                    return 0
            elif TypeFlag == 16:
                if S >= x:
                    return exp(-r * T) * S
                else:
                    return 0
            elif TypeFlag == 18:
                if S <= x:
                    return exp(-r * T) * k
                else:
                    return 0
            elif TypeFlag == 20:
                if S <= x:
                    return exp(-r * T) * S
                else:
                    return 0
            elif TypeFlag == 22:
                return 0
            elif TypeFlag == 24:
                return 0
            elif TypeFlag == 26:
                return 0
            elif TypeFlag == 28:
                return 0
            else:
                raise ValueError("invalid type flag")
        elif OutPutFlag == "d":
            if TypeFlag == 4:
                return 1
            elif TypeFlag == 8 or TypeFlag == 16 or TypeFlag == 20:
                return exp((b - r) * T) * 1
        else:
            return 0
    
    if S <= h and TypeFlag / 2 - int(TypeFlag / 2) != 0:
        if OutPutFlag == "p":
            if TypeFlag == 1:
                return k
            elif TypeFlag == 3:
                return S
            elif TypeFlag == 5:
                return exp(-r * T) * k
            elif TypeFlag == 7:
                return exp(-r * T) * S
            elif TypeFlag == 9:
                return 0
            elif TypeFlag == 11:
                return 0
            elif TypeFlag == 13:
                if S >= x:
                    return exp(-r * T) * k
                else:
                    return 0
            elif TypeFlag == 15:
                if S >= x:
                    return exp(-r * T) * S
                else:
                    return 0
            elif TypeFlag == 17:
                if S <= x:
                    return exp(-r * T) * k
                else:
                    return 0
            elif TypeFlag == 19:
                if S <= x:
                    return exp(-r * T) * S
                else:
                    return 0
            elif TypeFlag == 21:
                return 0
            elif TypeFlag == 23:
                return 0
            elif TypeFlag == 25:
                return 0
            elif TypeFlag == 27:
                return 0
        elif OutPutFlag == "d":
            if TypeFlag == 3:
                return 1
            elif TypeFlag == 7 or TypeFlag == 15 or TypeFlag == 19:
                return exp((b - r) * T) * 1
            else:
                raise ValueError("invalid type flag")
        else:
            return 0
    
    
    if OutPutFlag == "p": # Value
        return BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v, eta, phi)
    elif OutPutFlag == "d": # Delta
         return (BinaryBarrier(TypeFlag, S + dS, x, h, k, T, r, b, v, eta, phi) - BinaryBarrier(TypeFlag, S - dS, x, h, k, T, r, b, v, eta, phi)) / (2 * dS)
    elif OutPutFlag == "g": # Gamma
        return (BinaryBarrier(TypeFlag, S + dS, x, h, k, T, r, b, v, eta, phi) - 2 * BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v, eta, phi) + BinaryBarrier(TypeFlag, S - dS, x, h, k, T, r, b, v, eta, phi)) / dS ** 2
    elif OutPutFlag == "gv": # DGammaDVol
        return (BinaryBarrier(TypeFlag, S + dS, x, h, k, T, r, b, v + 0.01, eta, phi) - 2 * BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v + 0.01, eta, phi) + BinaryBarrier(TypeFlag, S - dS, x, h, k, T, r, b, v + 0.01, eta, phi) - BinaryBarrier(TypeFlag, S + dS, x, h, k, T, r, b, v - 0.01, eta, phi) + 2 * BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v - 0.01, eta, phi) - BinaryBarrier(TypeFlag, S - dS, x, h, k, T, r, b, v - 0.01, eta, phi)) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "dddv": # DDeltaDvol
        return 1 / (4 * dS * 0.01) * (BinaryBarrier(TypeFlag, S + dS, x, h, k, T, r, b, v + 0.01, eta, phi) - BinaryBarrier(TypeFlag, S + dS, x, h, k, T, r, b, v - 0.01, eta, phi) - BinaryBarrier(TypeFlag, S - dS, x, h, k, T, r, b, v + 0.01, eta, phi) + BinaryBarrier(TypeFlag, S - dS, x, h, k, T, r, b, v - 0.01, eta, phi)) / 100
    elif OutPutFlag == "v": # Vega
         return (BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v + 0.01, eta, phi) - BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v - 0.01, eta, phi)) / 2
    elif OutPutFlag == "vp": # VegaP
         return v / 0.1 * (BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v + 0.01, eta, phi) - BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v - 0.01, eta, phi)) / 2
    elif OutPutFlag == "dvdv": # DvegaDvol
        return (BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v + 0.01, eta, phi) - 2 * BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v, eta, phi) + BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v - 0.01, eta, phi))
    elif OutPutFlag == "t": # Theta
        if T <= 1 / 365:
            return BinaryBarrier(TypeFlag, S, x, h, k, 0.00001, r, b, v, eta, phi) - BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v, eta, phi)
        else:
            return BinaryBarrier(TypeFlag, S, x, h, k, T - 1 / 365, r, b, v, eta, phi) - BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v, eta, phi)
    elif OutPutFlag == "r": #Rho
        return (BinaryBarrier(TypeFlag, S, x, h, k, T, r + 0.01, b + 0.01, v, eta, phi) - BinaryBarrier(TypeFlag, S, x, h, k, T, r - 0.01, b - 0.01, v, eta, phi)) / (2)
    elif OutPutFlag == "fr": # Futures Rho
         return (BinaryBarrier(TypeFlag, S, x, h, k, T, r + 0.01, b, v, eta, phi) - BinaryBarrier(TypeFlag, S, x, h, k, T, r - 0.01, b, v, eta, phi)) / (2)
    elif OutPutFlag == "f": # Rho2
         return (BinaryBarrier(TypeFlag, S, x, h, k, T, r, b - 0.01, v, eta, phi) - BinaryBarrier(TypeFlag, S, x, h, k, T, r, b + 0.01, v, eta, phi)) / (2)
    elif OutPutFlag == "b": # Carry
        return (BinaryBarrier(TypeFlag, S, x, h, k, T, r, b + 0.01, v, eta, phi) - BinaryBarrier(TypeFlag, S, x, h, k, T, r, b - 0.01, v, eta, phi)) / (2)
    elif OutPutFlag == "s": # Speed
        return 1 / dS ** 3 * (BinaryBarrier(TypeFlag, S + 2 * dS, x, h, k, T, r, b, v, eta, phi) - 3 * BinaryBarrier(TypeFlag, S + dS, x, h, k, T, r, b, v, eta, phi) + 3 * BinaryBarrier(TypeFlag, S, x, h, k, T, r, b, v, eta, phi) - BinaryBarrier(TypeFlag, S - dS, x, h, k, T, r, b, v, eta, phi))
    else:
        raise ValueError("invalid output flag")


# Cash-or-nothing options
def CashOrNothing(CallPutFlag: Literal['c', 'p'], S: float, x: float, k: float, T: float, r: float, b: float, v: float) -> float:

    d = (log(S / x) + (b - v ** 2 / 2) * T) / (v * sqrt(T))

    if CallPutFlag == "c":
        return k * exp(-r * T) * CND(d)
    elif CallPutFlag == "p":
        return k * exp(-r * T) * CND(-d)
    else:
        raise ValueError("invalid call put flag")
