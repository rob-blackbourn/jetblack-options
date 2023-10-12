"""Plain Vanilla"""

from math import exp, log, pi, sqrt
from typing import Literal, Optional

from ..distributions import CND, ND, CNDEV, CHIINV


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
