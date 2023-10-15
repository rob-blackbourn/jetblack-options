"""Tests for distributions"""

from jetblack_options.distributions import CHIINV, CND

from .utils import is_close_to

def test_CHIINV():
    actual = CHIINV(0.050001, 10)
    expected = 18.306973
    assert is_close_to(actual, expected, 1e-6)

def test_cnd():
    actual = CND(1.23564285010596)
    expected = 0.89170432529243
    assert is_close_to(actual, expected, 1e-6)