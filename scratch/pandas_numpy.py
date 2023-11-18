"""Example implementation using pandas and numpy for vector calculation"""

from typing import Callable, Union

from numpy import exp, log, sqrt
from numpy.typing import NDArray
import numpy as np
import pandas as pd
from scipy.stats import norm

BoolVector = Union["pd.Series[bool]", NDArray[np.bool_]]
FloatVector = Union["pd.Series[float]", NDArray[np.float_]]

cdf = norm.cdf
pdf = norm.pdf


def price(
        is_call: BoolVector,
        S: FloatVector,
        K: FloatVector,
        T: FloatVector,
        r: FloatVector,
        b: FloatVector,
        v: FloatVector
) -> FloatVector:
    d1 = (log(S / K) + T * (b + v ** 2 / 2)) / (v * sqrt(T))
    d2 = d1 - v * sqrt(T)

    return np.where(
        is_call,
        S * exp((b - r) * T) * cdf(d1) - K * exp(-r * T) * cdf(d2),
        K * exp(-r * T) * cdf(-d2) - S * exp((b - r) * T) * cdf(-d1)
    )


def price_df(df, is_call='is_call', S='S', K='K', T='T', r='r', b='b', v='v'):
    return price(
        df[is_call],
        df[S],
        df[K],
        df[T],
        df[r],
        df[b],
        df[v]
    )


def solve_ivol(
        p: FloatVector,
        price: Callable[[FloatVector], FloatVector],
        *,
        max_iterations: int = 20,
        epsilon=1e-8
) -> NDArray[np.float_]:
    epsilon_ = np.full(len(p), epsilon)

    v_lo: NDArray[np.float_] = np.full(len(p), 0.005)
    v_hi: NDArray[np.float_] = np.full(len(p), 4.0)
    p_lo = price(v_lo)
    p_hi = price(v_hi)

    n = 0
    v = v_lo + (p - p_lo) * (v_hi - v_lo) / (p_hi - p_lo)
    p1 = price(v)
    while np.any(np.abs(p - p1) > epsilon_) and n < max_iterations:
        n += 1

        v_lo = np.where(p1 < p, v, v_lo)
        v_hi = np.where(p1 < p, v_hi, v)

        p_lo = price(v_lo)
        p_hi = price(v_hi)
        v = v_lo + (p - p_lo) * (v_hi - v_lo) / (p_hi - p_lo)
        p1 = price(v)

    return v


def ivol(
        df: pd.DataFrame,
        is_call='is_call',
        S='S',
        K='K',
        T='T',
        r='r',
        b='b',
        p='p',
        *,
        max_iterations: int = 20,
        epsilon=1e-8
) -> FloatVector:
    return solve_ivol(
        df[p],
        lambda v: price(df[is_call], df[S], df[K], df[T], df[r], df[b], v),
        max_iterations=max_iterations,
        epsilon=epsilon
    )


price_data = pd.DataFrame([
    {'is_call': True, 'S': 110, 'K': 100, 'r': 0.1,
        'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': False, 'S': 110, 'K': 100,
        'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': True, 'S': 100, 'K': 100, 'r': 0.1,
        'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': False, 'S': 100, 'K': 100,
        'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': True, 'S': 100, 'K': 110, 'r': 0.1,
        'q': 0.08, 'T': 6/12, 'v': 0.125},
    {'is_call': False, 'S': 100, 'K': 110,
        'r': 0.1, 'q': 0.08, 'T': 6/12, 'v': 0.125},
])

price_data['b'] = price_data['r'] - price_data['q']
price_result = price_df(price_data)
print(price_result)

ivol_data = price_data.copy()
ivol_data['p'] = price_result
ivol_result = ivol(ivol_data)
print(ivol_result)
