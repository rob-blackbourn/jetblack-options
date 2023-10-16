"""Double barriers"""

from typing import Callable

from ...european.black_scholes_merton import price as bs_price
from ...distributions import CDF
from ...numeric_greeks import NumericGreeks

from .analytics import price

def delta(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dS: float = 0.0001,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.delta(is_call, S, K, T, r, b, v, dS=dS)
    
    return (
        price(is_call, is_in, S + dS, K, L, U, T, r, b, v, delta1, delta2)
        - price(is_call, is_in, S - dS, K, L, U, T, r, b, v, delta1, delta2)
    ) / (2 * dS)


def gamma(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dS: float = 0.0001,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.gamma(is_call, S, K, T, r, b, v, dS=dS)
    
    return (
        price(is_call, is_in, S + dS, K, L, U, T, r, b, v, delta1, delta2)
        - 2 * price(is_call, is_in, S, K, L, U, T, r, b, v, delta1, delta2)
        + price(is_call, is_in, S - dS, K, L, U, T, r, b, v, delta1, delta2)
    ) / (dS ** 2)


def ddelta_dvol(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.ddelta_dvol(is_call, S, K, T, r, b, v, dS=dS, dv=dv)
    
    return 1 / (4 * dS * 0.01) * (
        price(is_call, is_in, S + dS, K, L, U, T, r, b, v + 0.01, delta1, delta2) -
        price(is_call, is_in, S + dS, K, L, U, T, r, b, v - 0.01, delta1, delta2) -
        price(is_call, is_in, S - dS, K, L, U, T, r, b, v + 0.01, delta1, delta2) +
        price(is_call, is_in, S - dS, K, L, U, T, r, b, v - 0.01, delta1, delta2)
        ) / 100

def gammap(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dS: float = 0.0001,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.gammap(is_call, S, K, T, r, b, v, dS=dS)
    
    return S / 100 * gamma(is_call, is_in, S + dS, K, L, U, T, r, b, v, delta1, delta2)

def dgamma_dvol(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.dgamma_dvol(is_call, S, K, T, r, b, v, dS=dS, dv=dv)
    
    return (
        price(is_call, is_in, S + dS, K, L, U, T, r, b, v + dv, delta1, delta2) -
        2 * price(is_call, is_in, S, K, L, U, T, r, b, v + dv, delta1, delta2) +
        price(is_call, is_in, S - dS, K, L, U, T, r, b, v + dv, delta1, delta2) -
        price(is_call, is_in, S + dS, K, L, U, T, r, b, v - dv, delta1, delta2) +
        2 * price(is_call, is_in, S, K, L, U, T, r, b, v - dv, delta1, delta2) -
        price(is_call, is_in, S - dS, K, L, U, T, r, b, v - dv, delta1, delta2)
    ) / (2 * dv * dS ** 2) / 100

def vega(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.vega(is_call, S, K, T, r, b, v, dv=dv)
    
    return (
        price(is_call, is_in, S, K, L, U, T, r, b, v + dv, delta1, delta2) -
        price(is_call, is_in, S, K, L, U, T, r, b, v - dv, delta1, delta2)
    ) / 2

def vomma(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.vomma(is_call, S, K, T, r, b, v, dv=dv)
    
    return (
        price(is_call, is_in, S, K, L, U, T, r, b, v + dv, delta1, delta2) -
        2 * price(is_call, is_in, S, K, L, U, T, r, b, v, delta1, delta2) +
        price(is_call, is_in, S, K, L, U, T, r, b, v - dv, delta1, delta2)
    ) / dv ** 2 / 10000

def vegap(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dS: float = 0.0001,
        dv: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.vegap(is_call, S, K, T, r, b, v, dv=dv)
    
    return v / 0.1 * vega(is_call, is_in, S + dS, K, L, U, T, r, b, v, delta1, delta2, dv=dv)

def rho(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dr: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.rho(is_call, S, K, T, r, b, v, dr=dr)
    
    return (
        price(is_call, is_in, S, K, L, U, T, r + dr, b + dr, v, delta1, delta2) -
        price(is_call, is_in, S, K, L, U, T, r - dr, b - dr, v, delta1, delta2)
    ) / 2

def futures_rho(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dr: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.futures_rho(is_call, S, K, T, r, b, v, dr=dr)
    
    return (
        price(is_call, is_in, S, K, L, U, T, r + dr, 0, v, delta1, delta2)
        - price(is_call, is_in, S, K, L, U, T, r - dr, 0, v, delta1, delta2)
    ) / 2

def rho2(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        db: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.carry(is_call, S, K, T, r, b, v, db=db)
    
    return (
        price(is_call, is_in, S, K, L, U, T, r, b - db, v, delta1, delta2)
        - price(is_call, is_in, S, K, L, U, T, r, b + db, v, delta1, delta2)
    ) / 2

def carry(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        db: float = 0.01,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.carry(is_call, S, K, T, r, b, v, db=db)
    
    return (
        price(is_call, is_in, S, K, L, U, T, r, b + db, v, delta1, delta2)
        - price(is_call, is_in, S, K, L, U, T, r, b - db, v, delta1, delta2)
    ) / 2

def theta(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dT: float = 1 / 365,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.theta(is_call, S, K, T, r, b, v, dT=dT)
    
    if T <= dT:
        return (
            price(is_call, is_in, S, K, L, U, 0.00001, r, b, v, delta1, delta2)
            - price(is_call, is_in, S, K, L, U, T, r, b, v, delta1, delta2)
        )
    else:
        return (
            price(is_call, is_in, S, K, L, U, T - dT, r, b, v, delta1, delta2)
            - price(is_call, is_in, S, K, L, U, T, r, b, v, delta1, delta2)
        )

def strike_delta(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dX: float = 0.0001,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.strike_delta(is_call, S, K, T, r, b, v, dX=dX)
    
    return (
        price(is_call, is_in, S, K + dX, L, U, T, r, b, v, delta1, delta2)
        - price(is_call, is_in, S, K - dX, L, U, T, r, b, v, delta1, delta2)
    ) / (2 * dX)

def strike_gamma(
        is_call: bool,
        is_in: bool,
        S: float,
        K: float,
        L: float,
        U: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta1: float,
        delta2: float,
        *,
        dX: float = 0.0001,
        cdf: Callable[[float], float] = CDF
) -> float:
    
    if not is_in and (S <= L or S >= U):
        return 0
    elif is_in and (S <= L or S >= U):
        bs = NumericGreeks(
            lambda is_call, S, K, T, r, b, v: bs_price(is_call, S, K, T, r, b, v, cdf=cdf)
        )
        return bs.strike_gamma(is_call, S, K, T, r, b, v, dX=dX)
    
    return (
        price(is_call, is_in, S, K + dX, L, U, T, r, b, v, delta1, delta2)
        - 2 * price(is_call, is_in, S, K, L, U, T, r, b, v, delta1, delta2)
        + price(is_call, is_in, S, K - dX, L, U, T, r, b, v, delta1, delta2)
    ) / (dX ** 2)
    