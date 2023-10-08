"""Black Scholes Merton"""

from typing import Literal

from .plain_vanilla import GBlackScholes

def ImpliedVolGBlackScholes(CallPutFlag: Literal['c', 'p'], S: float, X: float, T: float, r: float, b: float, cm: float) -> float:

    vLow = 0.005
    vHigh = 4
    epsilon = 0.00000001
    cLow = GBlackScholes(CallPutFlag, S, X, T, r, b, vLow)
    cHigh = GBlackScholes(CallPutFlag, S, X, T, r, b, vHigh)
    N = 0
    vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow)
    while abs(cm - GBlackScholes(CallPutFlag, S, X, T, r, b, vi)) > epsilon:
        N = N + 1
        if N > 20:
            break
        
        if GBlackScholes(CallPutFlag, S, X, T, r, b, vi) < cm:
            vLow = vi
        else:
            vHigh = vi

        cLow = GBlackScholes(CallPutFlag, S, X, T, r, b, vLow)
        cHigh = GBlackScholes(CallPutFlag, S, X, T, r, b, vHigh)
        vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow)

    return vi
