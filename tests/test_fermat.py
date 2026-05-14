import pytest
from factorizer import FermatFactorizer as Factorizer
from factorizer import TimeOutError


def test_ok_simple_case_timeout_empty():
    divider = Factorizer()
    facts = divider.factorize(57)
    assert facts == (3, 19)


def test_ok_simple_case_timeout_None():
    divider = Factorizer(timeout=None)
    facts = divider.factorize(57)
    assert facts == (3, 19)


def test_ok_simple_case_timeout_minus1():
    divider = Factorizer(timeout=-1)
    facts = divider.factorize(57)
    assert facts == (3, 19)


def test_ok_simple_case_timeout_5():
    divider = Factorizer(timeout=5)
    facts = divider.factorize(57)
    assert facts == (3, 19)


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
    assert facts == (1, 2) or facts == (2, 1)


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
    divider = Factorizer(timeout=5)
    n = 97**2
    facts = divider.factorize(n)
    assert facts == (97, 97)


def test_ok_prime():
    divider = Factorizer(timeout=5)
    facts = divider.factorize(97)
    assert facts == (1, 97)


def test_ng_both_big_timeout():
    divider = Factorizer(timeout=1)
    n = 144483604528043653279487
    with pytest.raises(TimeOutError):
        facts = divider.factorize(n)


def test_ok_max_iter_sufficient():
    divider = Factorizer(max_iter=10000)
    facts = divider.factorize(57)
    assert facts == (3, 19)


def test_ok_max_iter_none():
    divider = Factorizer(max_iter=None)
    facts = divider.factorize(57)
    assert facts == (3, 19)


def test_ng_max_iter_zero():
    # max_iter=0: no steps taken, returns (1, n) for any non-trivial composite
    divider = Factorizer(max_iter=0)
    n = 3 * 9999991
    d, q = divider.factorize(n)
    assert d * q == n
    assert d == 1


def test_ng_max_iter_too_small():
    # factors are far apart, Fermat needs many steps
    divider = Factorizer(max_iter=1)
    n = 3 * 9999991
    d, q = divider.factorize(n)
    assert d * q == n
    assert d == 1


def test_ok_both_very_big():
    divider = Factorizer(timeout=5)
    n = 17094896531810236860130284769490321915294047711310368136394905170386978916759334950349117125605720936295486193376534502996558346954298821541237678202872064373159864453975455700275218336138295073109837853052911442982249175483147894454735060228013876333548456070708161754022653432247114865681934678336369669630310417775703125135381535339600363990582556226874678712573925886711666498382111760727570796847908646797570295191362260110871270300928038772025787352463090061558150045455638071691578922379413464004696514186854056461649591756647526422083918100655183395966280798720122926503397649310956601613684952853243575941193
    facts = divider.factorize(n)
    assert facts == (130747453251718202865367599610984256344049976791406031337536420100072758990497257029481575106895700732897562067500860599588866977471195605540443529953122402713181190128830424125700106569684526840827976181922955582670800429750099675311621080502576741880521124560628733981949754441379686189254362646041261993641, 130747453251718202865367599610984256344049976791406031337536420100072758990497257029481575106895700732897562067500860599588866977471195605540443529953122402713181190128830424125700106569684526840827976181922955582670800429750099675311621080502576741880521124560628733981949754441379686189254362646041262989473)
