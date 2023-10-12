"""Asian options"""

from math import exp, log, sqrt
from typing import Literal, Optional

from .distributions import CND
from .european.black_scholes.plain_vanilla import GBlackScholes

def AsianCurranApprox(CallPutFlag: Literal['c', 'p'], S: float, SA: float, X: float, t1: float, T: float, n: int, m: int, r: float, b: float, v: float) -> float:

    z = 1
    if CallPutFlag == "p":
        z = -1

    dt = (T - t1) / (n - 1)
    
    if b == 0:
        EA = S
    else:
        EA = S / n * exp(b * t1) * (1 - exp(b * dt * n)) / (1 - exp(b * dt))

    if m > 0:
        if SA > n / m * X:
            # Exercise is certain for call, put must be out-of-the-money:
            if CallPutFlag == "p":
                return 0
            elif CallPutFlag == "c":
                SA = SA * m / n + EA * (n - m) / n
                return (SA - X) * exp(-r * T)

    if m == n - 1:
        # Only one fix left use Black-Scholes weighted with time
        X = n * X - (n - 1) * SA
        return GBlackScholes(CallPutFlag, S, X, T, r, b, v) * 1 / n

    if m > 0:
        X = n / (n - m) * X - m / (n - m) * SA

    vx = v * sqrt(t1 + dt * (n - 1) * (2 * n - 1) / (6 * n))
    my = log(S) + (b - v * v * 0.5) * (t1 + (n - 1) * dt / 2)

    sum1 = 0
    for i in range(1, 1+n):
    
        ti = dt * i + t1 - dt
        vi = v * sqrt(t1 + (i - 1) * dt)
        vxi = v * v * (t1 + dt * ((i - 1) - i * (i - 1) / (2 * n)))
        myi = log(S) + (b - v * v * 0.5) * ti
        sum1 = sum1 + exp(myi + vxi / (vx * vx) * (log(X) - my) + (vi * vi - vxi * vxi / (vx * vx)) * 0.5)
    Km = 2 * X - 1 / n * sum1
    sum2 = 0

    for i in range(1, 1+n):
    
        ti = dt * i + t1 - dt
        vi = v * sqrt(t1 + (i - 1) * dt)
        vxi = v * v * (t1 + dt * ((i - 1) - i * (i - 1) / (2 * n)))
        myi = log(S) + (b - v * v * 0.5) * ti
        sum2 = sum2 + exp(myi + vi * vi * 0.5) * CND(z * ((my - log(Km)) / vx + vxi / vx))

    return exp(-r * T) * z * (1 / n * sum2 - X * CND(z * (my - log(Km)) / vx)) * (n - m) / n

def EAsianCurranApprox(
        OutPutFlag: Literal['p', 'd'],
        CallPutFlag: Literal['c', 'p'],
        S: float,
        SA: float,
        X: float,
        t1: float,
        T: float,
        n: int,
        m: int,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.01

    if OutPutFlag == "p": # Value
        return AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v)
    elif OutPutFlag == "d": #Delta
        return (
            AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) -
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "e": #Elasticity
        return (
            AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) -
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        ) / (2 * dS) * S / AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v)
    elif OutPutFlag == "g": #Gamma
        return (
            AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) -
            2 * AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": #DGammaDVol
        return (
            AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v + 0.01) -
            2 * AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) +
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v + 0.01) -
            AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v - 0.01) +
            2 * AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01) -
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": #GammaP
        return S / 100 * (
            AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) -
            2 * AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "dddv": #DDeltaDvol
        return 1 / (4 * dS * 0.01) * (
            AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v + 0.01) -
            AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v - 0.01) -
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v + 0.01) +
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "v": #Vega
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) -
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "vv": #DvegaDvol/vomma
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) -
            2 * AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": #VegaP
        return v / 0.1 * (
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) -
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": #DvegaDvol
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) -
            2 * AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01)
        )
    elif OutPutFlag == "t": #Theta
        if t1 > 1 / 365 and T > 1 / 365:
            return (
                AsianCurranApprox(CallPutFlag, S, SA, X, t1 - 1 / 365, T - 1 / 365, n, m, r, b, v) -
                AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v)
            )
        else:
            raise ValueError("must have at least d day")
    elif OutPutFlag == "r": #Rho
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r + 0.01, b + 0.01, v) -
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r - 0.01, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "fr": #Futures options rho
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r + 0.01, b, v) -
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r - 0.01, b, v)
        ) / 2
    elif OutPutFlag == "f": #Rho2
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b - 0.01, v) -
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": #Carry
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b + 0.01, v) -
            AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "s": #Speed
        return 1 / dS ** 3 * (
            AsianCurranApprox(CallPutFlag, S + 2 * dS, SA, X, t1, T, n, m, r, b, v) -
            3 * AsianCurranApprox(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) +
            3 * AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) -
            AsianCurranApprox(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        )
    elif OutPutFlag == "dx": #Strike Delta
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X + dS, t1, T, n, m, r, b, v) -
            AsianCurranApprox(CallPutFlag, S, SA, X - dS, t1, T, n, m, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": #Gamma
        return (
            AsianCurranApprox(CallPutFlag, S, SA, X + dS, t1, T, n, m, r, b, v) -
            2 * AsianCurranApprox(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            AsianCurranApprox(CallPutFlag, S, SA, X - dS, t1, T, n, m, r, b, v)
        ) / dS ** 2
    else:
        raise ValueError("unknown output flag")

def DiscreteAsianHHM(CallPutFlag: Literal['c', 'p'], S: float, SA: float, X: float, t1: float, T: float, n: float, m: float, r: float, b: float, v: float) -> float:

    # This is a modified version of the Levy formula, this is the formula published in "Asian Pyramid Power" By
    #  Haug, Haug and Margrabe

    h = (T - t1) / (n - 1)

    if b == 0:
        EA = S
    else:
        EA = S / n * exp(b * t1) * (1 - exp(b * h * n)) / (1 - exp(b * h))
   
    if m > 0:
        if SA > n / m * X: #  Exercise is certain for call, put must be out-of-the-money
        
            if CallPutFlag == "p":
                return 0
            elif CallPutFlag == "c":
                SA = SA * m / n + EA * (n - m) / n
                return (SA - X) * exp(-r * T)
            else:
                raise ValueError("invalid call put flag")

    if m == n - 1: # Only one fix left use Black-Scholes weighted with time
   
        X = n * X - (n - 1) * SA
        return GBlackScholes(CallPutFlag, S, X, T, r, b, v) * 1 / n

    if b == 0:
        EA2 = (
            S * S * exp(v * v * t1) / (n * n) * 
            (
                (1 - exp(v * v * h * n)) / (1 - exp(v * v * h)) +
                2 / (1 - exp(v * v * h)) * (n - (1 - exp(v * v * h * n)) / (1 - exp(v * v * h)))
            )
        )
    else:
        EA2 = (
            S * S * exp((2 * b + v * v) * t1) / (n * n) *
            (
                (1 - exp((2 * b + v * v) * h * n)) / (1 - exp((2 * b + v * v) * h)) +
                2 / (1 - exp((b + v * v) * h)) *
                (
                    (1 - exp(b * h * n)) / (1 - exp(b * h)) -
                    (1 - exp((2 * b + v * v) * h * n)) / (1 - exp((2 * b + v * v) * h))
                )
            )
        )

    vA = sqrt((log(EA2) - 2 * log(EA)) / T)

    OptionValue = 0
    
    if m > 0:
        X = n / (n - m) * X - m / (n - m) * SA
    
    d1 = (log(EA / X) + vA ** 2 / 2 * T) / (vA * sqrt(T))
    d2 = d1 - vA * sqrt(T)

    if CallPutFlag == "c":
        OptionValue = exp(-r * T) * (EA * CND(d1) - X * CND(d2))
    elif CallPutFlag == "p":
        OptionValue = exp(-r * T) * (X * CND(-d2) - EA * CND(-d1))
    else:
        raise ValueError("invalid call put flag")

    return OptionValue * (n - m) / n

def EDiscreteAsianHHM(
        OutPutFlag: Literal['p', 'd'],
        CallPutFlag: Literal['c', 'p'],
        S: float,
        SA: float,
        X: float,
        t1: float,
        T: float,
        n: float,
        m: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.01
        
    if OutPutFlag == "p": # Value
        return DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v)
    elif OutPutFlag == "d": #Delta
        return (
            DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) -
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "e": #Elasticity
        return (
            DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) -
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        ) / (2 * dS) * S / DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v)
    elif OutPutFlag == "g": #Gamma
        return (
            DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) -
            2 * DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": #DGammaDVol
        return (
            DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v + 0.01) -
            2 * DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) +
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v + 0.01) -
            DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v - 0.01) +
            2 * DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01) -
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": #GammaP
        return S / 100 * (
            DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) -
            2 * DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "dddv": #DDeltaDvol
        return 1 / (4 * dS * 0.01) * (
            DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v + 0.01) -
            DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v - 0.01) -
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v + 0.01) +
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "v": #Vega
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) -
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "vv": #DvegaDvol/vomma
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) -
            2 * DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": #VegaP
        return v / 0.1 * (
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) -
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": #DvegaDvol
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v + 0.01) -
            2 * DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v - 0.01)
        )
    elif OutPutFlag == "t": #Theta
        if t1 > 1 / 365 and T > 1 / 365:
            return (
                DiscreteAsianHHM(CallPutFlag, S, SA, X, t1 - 1 / 365, T - 1 / 365, n, m, r, b, v) -
                DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v)
            )
        else:
            raise ValueError("must have at least one day to expiry")
    elif OutPutFlag == "r": #Rho
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r + 0.01, b + 0.01, v) -
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r - 0.01, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "fr": #Futures options rho
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r + 0.01, b, v) -
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r - 0.01, b, v)
        ) / 2
    elif OutPutFlag == "f": #Rho2
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b - 0.01, v) -
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": #Carry
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b + 0.01, v) -
            DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "s": #Speed
        return 1 / dS ** 3 * (
            DiscreteAsianHHM(CallPutFlag, S + 2 * dS, SA, X, t1, T, n, m, r, b, v) -
            3 * DiscreteAsianHHM(CallPutFlag, S + dS, SA, X, t1, T, n, m, r, b, v) +
            3 * DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) -
            DiscreteAsianHHM(CallPutFlag, S - dS, SA, X, t1, T, n, m, r, b, v)
        )
    elif OutPutFlag == "dx": #Strike Delta
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X + dS, t1, T, n, m, r, b, v) -
            DiscreteAsianHHM(CallPutFlag, S, SA, X - dS, t1, T, n, m, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": #Gamma
        return (
            DiscreteAsianHHM(CallPutFlag, S, SA, X + dS, t1, T, n, m, r, b, v) -
            2 * DiscreteAsianHHM(CallPutFlag, S, SA, X, t1, T, n, m, r, b, v) +
            DiscreteAsianHHM(CallPutFlag, S, SA, X - dS, t1, T, n, m, r, b, v)
        ) / dS ** 2
    else:
        raise ValueError("unknown output flag")
    

# Arithmetic average rate option
def TurnbullWakemanAsian(CallPutFlag: Literal['c', 'p'], S: float, SA: float, X: float, T: float, T2: float, r: float, b: float, v: float) -> float:

    # CallPutFlag = "c" for call and "p" for put option
    # S = Asset price
    # SA= Realized average so far
    # X = Strike price
    # t1 = Time to start of average period in years
    # T =  Time to maturity in years of option  T
    # T2 = Original time in average period in years, constant over life of option
    # r = risk-free rate
    # b = cost of carry underlying asset can be positive and negative
    # v = annualized volatility of asset price

    # tau: reminding time of average perios

    t1 = max(0, T - T2)
    tau = T2 - T
   
    if b == 0:    
        M1 = 1
    else:
        M1 = (exp(b * T) - exp(b * t1)) / (b * (T - t1))

    # Take into account when option wil  be exercised
    if tau > 0:
    
        if T2 / T * X - tau / T * SA < 0:
    
            if CallPutFlag == "c":
                SA = SA * (T2 - T) / T2 + S * M1 * T / T2 # Expected average at maturity
                return max(0, SA - X) * exp(-r * T)
            elif CallPutFlag == 'p':
                return 0
            else:
                raise ValueError('invalid call put flag')

    if b == 0: # Extended to hold for options on futures 16 May 1999 Espen Haug
       M2 = (
           2 * exp(v * v * T) / (v ** 4 * (T - t1) ** 2) -
           2 * exp(v * v * t1) * (1 + v * v * (T - t1)) / (v ** 4 * (T - t1) ** 2)
        )
    else:
        M2 = (
            2 * exp((2 * b + v * v) * T) / ((b + v * v) * (2 * b + v * v) * (T - t1) ** 2) +
            2 * exp((2 * b + v * v) * t1) / (b * (T - t1) ** 2) * (1 / (2 * b + v * v) - exp(b * (T - t1)) / (b + v * v))
        )

    bA = log(M1) / T
    vA = sqrt(log(M2) / T - 2 * bA)
    if tau > 0:
        X = T2 / T * X - tau / T * SA
        return GBlackScholes(CallPutFlag, S, X, T, r, bA, vA) * T / T2
    else:
        return GBlackScholes(CallPutFlag, S, X, T, r, bA, vA)


def ETurnbullWakemanAsian(
        OutPutFlag: Literal['p', 'd'],
        CallPutFlag: Literal['c', 'p'],
        S: float,
        SA: float,
        X: float,
        T: float,
        T2: float,
        r: float,
        b: float,
        v: float,
        dS: Optional[float] = None
) -> float:

    if dS is None:
        dS = 0.01

    if OutPutFlag == "p": # Value
        return TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v)
    elif OutPutFlag == "d": #Delta
        return (
            TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v) -
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "e": #Elasticity
        return (
            TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v) -
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v)
        ) / (2 * dS) * S / TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v)
    elif OutPutFlag == "g": #Gamma
        return (
            TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v) -
            2 * TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v) +
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "gv": #DGammaDVol
        return (
            TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v + 0.01) -
            2 * TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v + 0.01) +
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v + 0.01) -
            TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v - 0.01) +
            2 * TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v - 0.01) -
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v - 0.01)
        ) / (2 * 0.01 * dS ** 2) / 100
    elif OutPutFlag == "gp": #GammaP
        return S / 100 * (
            TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v) -
            2 * TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v) +
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v)
        ) / dS ** 2
    elif OutPutFlag == "dddv": #DDeltaDvol
        return 1 / (4 * dS * 0.01) * (
            TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v + 0.01) - 
            TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v - 0.01) - 
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v + 0.01) + 
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v - 0.01)
        ) / 100
    elif OutPutFlag == "v": #Vega
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v + 0.01) -
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v - 0.01)
            ) / 2
    elif OutPutFlag == "vv": #DvegaDvol/vomma
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v + 0.01) -
            2 * TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v) +
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v - 0.01)
        ) / 0.01 ** 2 / 10000
    elif OutPutFlag == "vp": #VegaP
        return v / 0.1 * (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v + 0.01) -
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v - 0.01)
        ) / 2
    elif OutPutFlag == "dvdv": #DvegaDvol
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v + 0.01) -
            2 * TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v) +
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v - 0.01)
        )
    elif OutPutFlag == "t": #Theta
        if T <= 1 / 365:
            return (
                TurnbullWakemanAsian(CallPutFlag, S, SA, X, 0.00001, T2, r, b, v) -
                TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v)
            )
        else:
            return (
                TurnbullWakemanAsian(CallPutFlag, S, SA, X, T - 1 / 365, T2, r, b, v) -
                TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v)
            )
    elif OutPutFlag == "r": #Rho
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r + 0.01, b + 0.01, v) -
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r - 0.01, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "fr": #Futures options rho
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r + 0.01, b, v) -
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r - 0.01, b, v)
        ) / 2
    elif OutPutFlag == "f": #Rho2
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b - 0.01, v) -
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b + 0.01, v)
        ) / 2
    elif OutPutFlag == "b": #Carry
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b + 0.01, v) -
            TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b - 0.01, v)
        ) / 2
    elif OutPutFlag == "s": #Speed
        return 1 / dS ** 3 * (
            TurnbullWakemanAsian(CallPutFlag, S + 2 * dS, SA, X, T, T2, r, b, v) -
            3 * TurnbullWakemanAsian(CallPutFlag, S + dS, SA, X, T, T2, r, b, v) +
            3 * TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v) -
            TurnbullWakemanAsian(CallPutFlag, S - dS, SA, X, T, T2, r, b, v)
        )
    elif OutPutFlag == "dx": #Strike Delta
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X + dS, T, T2, r, b, v) -
            TurnbullWakemanAsian(CallPutFlag, S, SA, X - dS, T, T2, r, b, v)
        ) / (2 * dS)
    elif OutPutFlag == "dxdx": #Gamma
        return (
            TurnbullWakemanAsian(CallPutFlag, S, SA, X + dS, T, T2, r, b, v) -
            2 * TurnbullWakemanAsian(CallPutFlag, S, SA, X, T, T2, r, b, v) +
            TurnbullWakemanAsian(CallPutFlag, S, SA, X - dS, T, T2, r, b, v)
        ) / dS ** 2
    else:
        raise ValueError("invalid output flag")
