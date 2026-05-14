import pytest
from factorizer import ECMFactorizer as Factorizer
from factorizer import TimeOutError


def test_ok_simple_case_timeout_empty():
    divider = Factorizer()
    facts = divider.factorize(57)
    assert facts[0] * facts[1] == 57
    assert sorted(facts) == [3, 19]


def test_ok_simple_case_timeout_None():
    divider = Factorizer(timeout=None)
    facts = divider.factorize(57)
    assert facts[0] * facts[1] == 57
    assert sorted(facts) == [3, 19]


def test_ok_simple_case_timeout_minus1():
    divider = Factorizer(timeout=-1)
    facts = divider.factorize(57)
    assert facts[0] * facts[1] == 57
    assert sorted(facts) == [3, 19]


def test_ok_simple_case_timeout_5():
    divider = Factorizer(timeout=5)
    facts = divider.factorize(57)
    assert facts[0] * facts[1] == 57
    assert sorted(facts) == [3, 19]


def test_ok_0():
    divider = Factorizer()
    facts = divider.factorize(0)
    assert facts == (0, 0)


def test_ok_1():
    divider = Factorizer()
    facts = divider.factorize(1)
    assert facts == (1, 1)


def test_ok_2():
    divider = Factorizer()
    facts = divider.factorize(2)
    assert facts == (2, 1) or facts == (1, 2)


def test_ok_3():
    divider = Factorizer()
    facts = divider.factorize(3)
    assert facts == (1, 3) or facts == (3, 1)


def test_ok_4():
    divider = Factorizer()
    facts = divider.factorize(4)
    assert facts == (2, 2)


def test_ok_6():
    divider = Factorizer()
    facts = divider.factorize(6)
    assert facts == (2, 3) or facts == (3, 2)


def test_ok_even():
    divider = Factorizer(timeout=5)
    n = 2 * 3571
    facts = divider.factorize(n)
    assert facts == (2, 3571)


def test_ok_square():
    divider = Factorizer(timeout=10)
    n = 97 ** 2
    facts = divider.factorize(n)
    assert facts[0] * facts[1] == n
    assert sorted(facts) == [97, 97]


def test_ok_prime():
    divider = Factorizer(timeout=10)
    facts = divider.factorize(97, max_curves=50)
    assert facts == (1, 97) or facts == (97, 1)


def test_ok_medium():
    divider = Factorizer(timeout=30)
    n = 2885500349
    facts = divider.factorize(n, B1=10000, B2=100000, max_curves=200)
    assert facts[0] * facts[1] == n
    assert sorted(facts) == sorted([48611, 59359])


def test_ok_large_with_small_factor():
    divider = Factorizer(timeout=60)
    n = 18446744073709551617  # 2^64 + 1 = 274177 * 67280421310721
    facts = divider.factorize(n, B1=10000, B2=100000, max_curves=200)
    assert facts[0] * facts[1] == n
    assert 1 not in facts


def test_ok_custom_B1():
    divider = Factorizer(timeout=10)
    n = 221  # 13 * 17
    facts = divider.factorize(n, B1=100, B2=1000, max_curves=50)
    assert facts[0] * facts[1] == n
    assert sorted(facts) == [13, 17]


def test_ok_product_of_two_primes():
    divider = Factorizer(timeout=30)
    n = 1000000007 * 1000000009
    facts = divider.factorize(n, B1=50000, B2=500000, max_curves=300)
    assert facts[0] * facts[1] == n
    assert 1 not in facts


def test_ng_timeout():
    divider = Factorizer(timeout=1)
    n = 17094896531810236860130284769490321915294047711310368136394905170386978916759334950349117125605720936295486193376534502996558346954298821541237678202872064373159864453975455700275218336138295073109837853052911442982249175483147894454735060228013876333548456070708161754022653432247114865681934678336369669630310417775703125135381535339600363990582556226874678712573925886711666498382111760727570796847908646797570295191362260110871270300928038772025787352463090061558150045455638071691578922379413464004696514186854056461649591756647526422083918100655183395966280798720122926503397649310956601613684952853243575941193
    with pytest.raises(TimeOutError):
        divider.factorize(n, B1=100000, B2=1000000, max_curves=10000)
