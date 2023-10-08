"""Test utilities"""


def is_close_to(actual: float, expected: float, threshold: float) -> bool:
    return abs(actual - expected) < threshold
