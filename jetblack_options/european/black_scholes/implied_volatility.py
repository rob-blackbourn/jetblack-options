"""Black Scholes Merton"""

from typing import Callable

from ...distributions import CND

from .analytic import price

def ivol(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        cm: float,
        *,
        cnd: Callable[[float], float] = CND
) -> float:

    vLow = 0.005
    vHigh = 4
    epsilon = 0.00000001
    cLow = price(is_call, S, X, T, r, b, vLow, cnd=cnd)
    cHigh = price(is_call, S, X, T, r, b, vHigh, cnd=cnd)
    N = 0
    vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow)
    while abs(cm - price(is_call, S, X, T, r, b, vi, cnd=cnd)) > epsilon:
        N = N + 1
        if N > 20:
            break
        
        if price(is_call, S, X, T, r, b, vi, cnd=cnd) < cm:
            vLow = vi
        else:
            vHigh = vi

        cLow = price(is_call, S, X, T, r, b, vLow, cnd=cnd)
        cHigh = price(is_call, S, X, T, r, b, vHigh, cnd=cnd)
        vi = vLow + (cm - cLow) * (vHigh - vLow) / (cHigh - cLow)

    return vi
