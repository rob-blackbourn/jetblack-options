"""Plain Vanilla"""

from math import exp, log, pi, sqrt
from typing import Literal, Optional

from ...distributions import CND, ND, CNDEV, CHIINV


def price(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:

    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return S * exp((b - r) * T) * CND(d1) - X * exp(-r * T) * CND(d2)
    else:
        return X * exp(-r * T) * CND(-d2) - S * exp((b - r) * T) * CND(-d1)


def delta(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:

    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))

    if is_call:
        return exp((b - r) * T) * CND(d1)
    else:
        return -exp((b - r) * T) * CND(-d1)


def forward_delta(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:

    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))

    if is_call:
        return exp(-r * T) * CND(d1)
    else:
        return exp(-r * T) * (CND(d1) - 1)

# DDeltaDvol also known as vanna


def ddelta_dvol(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return -exp((b - r) * T) * d2 / v * ND(d1)

# DDeltaDvolDvol also known as DVannaDvol


def ddelta_dvol_dvol(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v * v / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return ddelta_dvol(S, X, T, r, b, v) * 1 / v * (d1 * d2 - d1 / d2 - 1)

# DdeltaDtime/Charm for the generalized Black and Scholes formula


def ddelta_dtime(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:

    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return -exp((b - r) * T) * (
            ND(d1) * (b / (v * sqrt(T)) - d2 / (2 * T)) + (b - r) * CND(d1)
        )
    else:
        return -exp((b - r) * T) * (
            ND(d1) * (b / (v * sqrt(T)) - d2 / (2 * T)) - (b - r) * CND(-d1)
        )


# SaddleGamma for the generalized Black and Scholes formula
def saddle_gamma(X: float, T: float, r: float, b: float, v: float) -> float:
    return sqrt(exp(1) / pi) * sqrt((2 * b - r) / v ** 2 + 1) / X


# Gamma for the generalized Black and Scholes formula
def gamma(S: float, X: float, T: float, r: float, b: float, v: float) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return exp((b - r) * T) * ND(d1) / (S * v * sqrt(T))


# DgammaDspot/Speed for the generalized Black and Scholes formula
def dgamma_dspot(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -gamma(S, X, T, r, b, v) * (1 + d1 / (v * sqrt(T))) / S


# DgammaDvol/Zomma for the generalized Black and Scholes formula
def dgamma_dvol(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return gamma(S, X, T, r, b, v) * ((d1 * d2 - 1) / v)


# GGammaDtime for the generalized Black and Scholes formula
def dgamma_dtime(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return gamma(S, X, T, r, b, v) * (
        r - b + b * d1 / (v * sqrt(T)) +
        (1 - d1 * d2) / (2 * T)
    )


# GammaP for the generalized Black and Scholes formula
def gammap(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    return S * gamma(S, X, T, r, b, v) / 100


# DgammaPDspot/SpeedP for the generalized Black and Scholes formula
def dgammap_dspot(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -gamma(S, X, T, r, b, v) * (d1) / (100 * v * sqrt(T))

# DgammaPDvol for the generalized Black and Scholes formula


def dgammap_dvol(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S / 100 * gamma(S, X, T, r, b, v) * ((d1 * d2 - 1) / v)

# GGammaPDtime for the generalized Black and Scholes formula


def dgammap_dtime(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        gammap(S, X, T, r, b, v) *
        (r - b + b * d1 / (v * sqrt(T)) + (1 - d1 * d2) / (2 * T))
    )


# Vega for the generalized Black and Scholes formula
def vega(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * exp((b - r) * T) * ND(d1) * sqrt(T)


# DvegaDtime for the generalized Black and Scholes formula
def dvega_dtime(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        vega(S, X, T, r, b, v) *
        (r - b + b * d1 / (v * sqrt(T)) - (1 + d1 * d2) / (2 * T))
    )

# DvegaDvol/Vomma for the generalized Black and Scholes formula


def dvega_dvol(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return vega(S, X, T, r, b, v) * d1 * d2 / v


# DVommaDVol for the generalized Black and Scholes formula
def dvomma_dvol(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return (
        dvega_dvol(S, X, T, r, b, v) *
        1 / v * (d1 * d2 - d1 / d2 - d2 / d1 - 1)
    )


# VegaP for the generalized Black and Scholes formula
def vegap(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    return v / 10 * vega(S, X, T, r, b, v)


# DvegaPDvol/VommaP for the generalized Black and Scholes formula
def dvegap_dvol(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return vegap(S, X, T, r, b, v) * d1 * d2 / v


# Vega for the generalized Black and Scholes formula
def vega_leverage(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    return (
        vega(S, X, T, r, b, v) * v /
        price(is_call, S, X, T, r, b, v)
    )

# Variance-vega for the generalized Black and Scholes formula


def variance_vega(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return S * exp((b - r) * T) * ND(d1) * sqrt(T) / (2 * v)


# Variance-delta for the generalized Black and Scholes formula
def variance_delta(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S * exp((b - r) * T) * ND(d1) * (-d2) / (2 * v ** 2)


# Variance-vomma for the generalized Black and Scholes formula
def variance_vomma(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    return S * exp((b - r) * T) * sqrt(T) / (4 * v ** 3) * ND(d1) * (d1 * d2 - 1)


# Variance-ultima for the generalized Black and Scholes formula
def variance_ultima(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
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
def theta(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    if is_call:
        return (
            -S * exp((b - r) * T) * ND(d1) * v / (2 * sqrt(T)) -
            (b - r) * S * exp((b - r) * T) * CND(d1) -
            r * X * exp(-r * T) * CND(d2)
        )
    else:
        return (
            -S * exp((b - r) * T) * ND(d1) * v / (2 * sqrt(T)) +
            (b - r) * S * exp((b - r) * T) * CND(-d1) +
            r * X * exp(-r * T) * CND(-d2)
        )


# Drift-less Theta for the generalized Black and Scholes formula
def theta_driftless(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    return -S * exp((b - r) * T) * ND(d1) * v / (2 * sqrt(T))


# Rho for the generalized Black and Scholes formula for all options except futures
def rho(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return T * X * exp(-r * T) * CND(d2)
    else:
        return -T * X * exp(-r * T) * CND(-d2)


# Rho for the generalized Black and Scholes formula for Futures option
def futures_rho(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    return -T * price(is_call, S, X, T, r, 0, v)


# Carry rho sensitivity for the generalized Black and Scholes formula
def carry(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return T * S * exp((b - r) * T) * CND(d1)
    else:
        return -T * S * exp((b - r) * T) * CND(-d1)


# Rho2/Phi for the generalized Black and Scholes formula
def phi(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return -T * S * exp((b - r) * T) * CND(d1)
    else:
        return T * S * exp((b - r) * T) * CND(-d1)


# DZetaDvol for the generalized Black and Scholes formula
def dzeta_dvol(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return -ND(d2) * d1 / v
    else:
        return ND(d2) * d1 / v


# DZetaDtime for the generalized Black and Scholes formula
def dzeta_dtime(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d1 = (log(S / X) + (b + v ** 2 / 2) * T) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)
    if is_call:
        return ND(d2) * (b / (v * sqrt(T)) - d1 / (2 * T))
    else:
        return -ND(d2) * (b / (v * sqrt(T)) - d1 / (2 * T))


# Risk neutral break even probability for the generalized Black and Scholes formula
def break_even_probability(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    if is_call:
        X = X + price(True, S, X, T, r, b, v) * exp(r * T)
        d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
        return CND(d2)
    else:
        X = X - price(False, S, X, T, r, b, v) * exp(r * T)
        d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
        return CND(-d2)


# StrikeDelta for the generalized Black and Scholes formula
def strike_delta(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    if is_call:
        return -exp(-r * T) * CND(d2)
    else:
        return exp(-r * T) * CND(-d2)


# Risk Neutral Density for the generalized Black and Scholes formula
def risk_neutral_density(
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))
    return exp(-r * T) * ND(d2) / (X * v * sqrt(T))


# Gamma from delta
def gamma_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta: float
) -> float:
    return exp((b - r) * T) * ND(
        CNDEV(exp((r - b) * T) * abs(delta))
    ) / (S * v * sqrt(T))


# GammaP from delta
def gammap_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta: float
) -> float:
    return S / 100 * gamma_from_delta(S, T, r, b, v, delta)


# Vega from delta
def vega_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        delta: float
) -> float:
    return S * exp((b - r) * T) * sqrt(T) * ND(
        CNDEV(exp((r - b) * T) * abs(delta))
    )


# VegaP from delta
def vegap_from_delta(
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta: float
) -> float:
    return v / 10 * vega_from_delta(S, T, r, b, delta)


# Closed form solution to find strike given the delta
def strike_from_delta(
        is_call: bool,
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta: float
) -> float:
    if is_call:
        return S * exp(
            -CNDEV(delta * exp((r - b) * T)) *
            v * sqrt(T) + (b + v * v / 2) * T
        )
    else:
        return S * exp(
            CNDEV(-delta * exp((r - b) * T)) *
            v * sqrt(T) + (b + v * v / 2) * T
        )


# Closed form solution to find in-the-money risk-neutral probaility given the delta
def in_the_money_prob_from_delta(
        is_call: bool,
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        delta: float
) -> float:

    if is_call:
        return CND(CNDEV(delta / exp((b - r) * T)) - v * sqrt(T))
    else:
        return CND(CNDEV(-delta / exp((b - r) * T)) + v * sqrt(T))


# Closed form solution to find strike given the in-the-money risk neutral probability
def strike_from_in_the_money_prob(
        is_call: bool,
        S: float,
        v: float,
        T: float,
        b: float,
        in_the_money_prob: float
) -> float:
    if is_call:
        return S * exp(
            -CNDEV(in_the_money_prob) * v * sqrt(T) + (b - v * v / 2) * T
        )
    else:
        return S * exp(
            CNDEV(in_the_money_prob) * v * sqrt(T) + (b - v * v / 2) * T
        )


# Risk Neutral Density from in-the-money probability
def rnd_from_in_the_money_prob(
        X: float,
        T: float,
        r: float,
        v: float,
        in_the_money_prob: float
) -> float:
    return exp(-r * T) * ND(CNDEV(in_the_money_prob)) / (X * v * sqrt(T))


# Closed form solution to find in-the-money risk-neutral probaility given the delta
def delta_from_in_the_money_prob(
        is_call: bool,
        S: float,
        T: float,
        r: float,
        b: float,
        v: float,
        in_the_money_prob: float
) -> float:
    if is_call:
        return CND(CNDEV(in_the_money_prob * exp((b - r) * T)) - v * sqrt(T))
    else:
        return -CND(CNDEV(in_the_money_prob * exp((b - r) * T)) + v * sqrt(T))


# What asset price that gives maximum DdeltaDvol
def max_ddelta_dvol_asset(is_lower: bool, x: float, T: float, b: float, v: float) -> float:
    # UpperLowerFlag"l" gives lower asset level that gives max DdeltaDvol
    # UpperLowerFlag"l" gives upper asset level that gives max DdeltaDvol

    if is_lower:
        return x * exp(-b * T - v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)
    else:
        return x * exp(-b * T + v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)

# What strike price that gives maximum DdeltaDvol


def max_ddelta_dvol_strike(is_lower: bool, S: float, T: float, b: float, v: float) -> float:

    # UpperLowerFlag"l" gives lower strike level that gives max DdeltaDvol
    # UpperLowerFlag"l" gives upper strike level that gives max DdeltaDvol

    if is_lower:
        return S * exp(b * T - v * sqrt(T) * sqrt(4 + T * v * 2) / 2)
    else:
        return S * exp(b * T + v * sqrt(T) * sqrt(4 + T * v ** 2) / 2)

# What strike price that gives maximum gamma and vega


def max_gamma_vega_at_X(S: float, b: float, T: float, v: float) -> float:
    return S * exp((b + v * v / 2) * T)

# What asset price that gives maximum gamma


def max_gamma_at_S(x: float, b: float, T: float, v: float) -> float:
    return x * exp((-b - 3 * v * v / 2) * T)


# What asset price that gives maximum vega
def max_vega_at_S(X: float, b: float, T: float, v: float) -> float:
    return X * exp((-b + v * v / 2) * T)

# Delta for the generalized Black and Scholes formula


def in_the_money_probability(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        b: float,
        v: float
) -> float:
    d2 = (log(S / X) + (b - v ** 2 / 2) * T) / (v * sqrt(T))

    if is_call:
        return CND(d2)
    else:
        return CND(-d2)

# MirrorDeltaStrike, delta neutral straddle strike in the BSM formula


def delta_mirror_strike(S: float, T: float, b: float, v: float) -> float:
    return S * exp((b + v ** 2 / 2) * T)

# MirrorProbabilityStrike, probability neutral straddle strike in the BSM formula


def probability_mirror_strike(S: float, T: float, b: float, v: float) -> float:
    return S * exp((b - v ** 2 / 2) * T)

# MirrorDeltaStrike, general delta symmmetric strike in the BSM formula


def delta_mirror_call_put_strike(S: float, x: float, T: float, b: float, v: float) -> float:
    return S ** 2 / x * exp((2 * b + v ** 2) * T)

# Elasticity for the generalized Black and Scholes formula


def elasticity(
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float
) -> float:
    return (
        delta(is_call, S, X, T, r, b, v) * S /
        price(is_call, S, X, T, r, b, v)
    )

# Volatility estimate confidence interval


def confidence_interval_volatility(
        Alfa: float,
        n: int,
        VolatilityEstimate: float,
        is_lower: bool
) -> float:
    # UpperLower     ="L" gives the lower cofidence interval
    #                ="U" gives the upper cofidence interval
    # n: number of observations
    if is_lower:
        return VolatilityEstimate * sqrt((n - 1) / (CHIINV(Alfa / 2, n - 1)))
    else:
        return VolatilityEstimate * sqrt((n - 1) / (CHIINV(1 - Alfa / 2, n - 1)))


# Profit/Loss STD for the generalized Black and Scholes formula
def profit_loss_std(
        TypeFlag: Literal['a', 'p'],
        is_call: bool,
        S: float,
        X: float,
        T: float,
        r: float,
        b: float,
        v: float,
        NHedges: int
) -> float:

    if TypeFlag == "a":  # in dollars
        return (
            sqrt(pi / 4) *
            vega(S, X, T, r, b, v) *
            v / sqrt(NHedges)
        )
    elif TypeFlag == "p":  # in percent
        return (
            sqrt(pi / 4) *
            vega(S, X, T, r, b, v) *
            v / sqrt(NHedges) /
            price(is_call, S, X, T, r, b, v)
        )
    else:
        raise ValueError('invalid type flag')
