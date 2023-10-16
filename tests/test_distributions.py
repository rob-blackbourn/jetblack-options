"""Tests for distributions"""

from statistics import NormalDist

from jetblack_options.distributions import CHIINV, CND, ND, CNDEV

from .utils import is_close_to

# def test_CHIINV():
#     actual = CHIINV(0.050001, 10)
#     expected = 18.306973
#     assert is_close_to(actual, expected, 1e-6)

def test_cnd():
    actual = CND(1.23564285010596)
    expected = 0.89170432529243
    assert is_close_to(actual, expected, 1e-6)


def test_nd():
    actual = ND(1.23564285010596)
    expected = 0.1859374114061063
    assert is_close_to(actual, expected, 1e-6)

def test_cndev():
    actual = CNDEV(0.7)
    expected = 0.5244005119066525
    assert is_close_to(actual, expected, 1e-6)


def test_cdf():
    nd = NormalDist()
    actual = nd.cdf(1.23564285010596)
    expected = 0.89170432529243
    assert is_close_to(actual, expected, 1e-6)


def test_pdf():
    nd = NormalDist()
    actual = nd.pdf(1.23564285010596)
    expected = 0.1859374114061063
    assert is_close_to(actual, expected, 1e-6)


def test_inv_cdf():
    nd = NormalDist()
    actual = nd.inv_cdf(0.7)
    expected = 0.5244005119066525
    assert is_close_to(actual, expected, 1e-6)
