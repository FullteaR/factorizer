import pytest
from factorizer import PplusOneFactorizer as Factorizer
from factorizer import TimeOutError


def test_ok_0():
    assert Factorizer().factorize(0) == (0, 0)


def test_ok_1():
    assert Factorizer().factorize(1) == (1, 1)


def test_ok_2():
    d, q = Factorizer().factorize(2)
    assert d * q == 2


def test_ok_3():
    d, q = Factorizer().factorize(3)
    assert d * q == 3


def test_ok_4():
    assert Factorizer().factorize(4) == (2, 2)


def test_ok_5():
    d, q = Factorizer().factorize(5)
    assert d * q == 5


def test_ok_even():
    d, q = Factorizer().factorize(6)
    assert {d, q} == {2, 3}


def test_ok_timeout_none():
    d, q = Factorizer(timeout=None).factorize(1081, M=10)
    assert {d, q} == {23, 47}


def test_ok_timeout_minus1():
    d, q = Factorizer(timeout=-1).factorize(1081, M=10)
    assert {d, q} == {23, 47}


def test_ok_timeout_5():
    d, q = Factorizer(timeout=5).factorize(1081, M=10)
    assert {d, q} == {23, 47}


def test_smooth_p_plus_1_24():
    # p=23, p+1=24=2^3*3; q=47, q+1=48=2^4*3 does NOT divide product(M=10)=2520
    d, q = Factorizer(A=3).factorize(1081, M=10)
    assert {d, q} == {23, 47}


def test_smooth_p_plus_1_14():
    # p=13, p+1=14=2*7; q=53, q+1=54=2*3^3 does NOT divide product(M=10)=2520
    d, q = Factorizer(A=3).factorize(689, M=10)
    assert {d, q} == {13, 53}


def test_smooth_p_plus_1_32():
    # p=31, p+1=32=2^5
    d, q = Factorizer(A=3).factorize(19003, M=300)
    assert {d, q} == {31, 613}


def test_smooth_p_plus_1_80():
    # p=79, p+1=80=2^4*5
    d, q = Factorizer(A=3).factorize(48427, M=300)
    assert {d, q} == {79, 613}


def test_smooth_p_plus_1_mersenne():
    # p=127, p+1=128=2^7
    d, q = Factorizer(A=3).factorize(16637, M=64)
    assert {d, q} == {127, 131}


def test_alternative_A():
    # A=3 with large M hits d=n for n=1135=5*227 (both factors' orders covered simultaneously)
    # A=5 succeeds because p=5 and q=227 orders do not coincide under this starting point
    d3, q3 = Factorizer(A=3).factorize(1135, M=100)
    assert {d3, q3} != {5, 227}
    d5, q5 = Factorizer(A=5).factorize(1135, M=10)
    assert {d5, q5} == {5, 227}


def test_timeout():
    with pytest.raises(TimeOutError):
        Factorizer(timeout=0.001).factorize(667, M=10**9)
