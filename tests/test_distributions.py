"""Tests for distributions"""

from jetblack_options.distributions import CHIINV

from .utils import is_close_to

def test_CHIINV():
    actual = CHIINV(0.050001, 10)
    expected = 18.306973
    assert is_close_to(actual, expected, 1e-6)