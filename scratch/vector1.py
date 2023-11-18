"""Vector style"""

from typing import Union, Any, Callable

import numpy as np
import pandas as pd
from scipy.stats import norm

Vector = Union[pd.Series, np.ndarray]


def _price(is_call: Vector, S: Vector, K: Vector, T: Vector, r: Vector, b: Vector, v: Vector) -> np.ndarray:
    d1 = (np.log(S / K) + T * (b + v ** 2 / 2)) / (v * np.sqrt(T))
    d2 = d1 - v * np.sqrt(T)

    return np.where(
        is_call,
        S * np.exp((b - r) * T) * norm.cdf(d1) -
        K * np.exp(-r * T) * norm.cdf(d2),
        K * np.exp(-r * T) * norm.cdf(-d2) - S *
        np.exp((b - r) * T) * norm.cdf(-d1)
    )


def price(df, is_call='is_call', S='S', K='K', T='T', r='r', b='b', v='v'):
    return _price(
        df[is_call],
        df[S],
        df[K],
        df[T],
        df[r],
        df[b],
        df[v]
    )


def solve_ivol(
        p: Vector,
        price: Callable[[Vector], Vector],
        *,
        max_iterations: int = 20,
        epsilon=1e-8
) -> Vector:
    epsilon_ = np.full(len(p), epsilon)

    v_lo = np.full(len(p), 0.005)
    v_hi = np.full(len(p), 4.0)
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
) -> Vector:
    return solve_ivol(
        df[p],
        lambda v: _price(df[is_call], df[S], df[K], df[T], df[r], df[b], v),
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
price_result = price(price_data)
print(price_result)

ivol_data = price_data.copy()
ivol_data['p'] = price_result
ivol_result = ivol(ivol_data)
print(ivol_result)