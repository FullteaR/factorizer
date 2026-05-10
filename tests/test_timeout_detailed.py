import pytest
import time
from factorizer import (
    BruteForceFactorizer,
    FermatFactorizer,
    PollardsRhoFactorizer,
    TimeOutError,
)

# BruteForce needs well over 1 s for this: 144483604528043653279487 = 2147483647 × 67280421310721
HARD_N = 144483604528043653279487
EASY_N = 57  # = 3 × 19


# ---------------------------------------------------------------------------
# Float timeout
# Currently broken: range(float * 100) raises TypeError in Python 3.
# e.g. timeout=1.5 → range(150.0) → TypeError
# Affects all C++-backed factorizers (FactorDB passes timeout to requests directly).
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("Factorizer", [BruteForceFactorizer, FermatFactorizer, PollardsRhoFactorizer])
def test_float_timeout_completes(Factorizer):
    divider = Factorizer(timeout=5.0)
    d, q = divider.factorize(EASY_N)
    assert d * q == EASY_N


@pytest.mark.parametrize("Factorizer", [BruteForceFactorizer, FermatFactorizer, PollardsRhoFactorizer])
def test_float_timeout_fractional_completes(Factorizer):
    divider = Factorizer(timeout=1.5)
    d, q = divider.factorize(EASY_N)
    assert d * q == EASY_N


def test_float_timeout_fires():
    divider = BruteForceFactorizer(timeout=1.5)
    with pytest.raises(TimeOutError):
        divider.factorize(HARD_N)


# ---------------------------------------------------------------------------
# timeout=0 edge case
# range(0 * 100) = range(0) → empty loop → else clause → immediate TimeOutError
# ---------------------------------------------------------------------------

def test_timeout_zero_raises_immediately():
    divider = BruteForceFactorizer(timeout=0)
    with pytest.raises(TimeOutError):
        divider.factorize(EASY_N)


# ---------------------------------------------------------------------------
# Timing accuracy
# The current implementation polls at 10 ms intervals.
# TimeOutError should fire in [timeout, timeout + some_slack].
# ---------------------------------------------------------------------------

def test_timeout_fires_within_2x():
    timeout = 1
    divider = BruteForceFactorizer(timeout=timeout)
    start = time.monotonic()
    with pytest.raises(TimeOutError):
        divider.factorize(HARD_N)
    elapsed = time.monotonic() - start
    assert elapsed < timeout * 2, f"fired at {elapsed:.3f}s, expected <{timeout * 2}s"


def test_timeout_not_too_early():
    timeout = 1
    divider = BruteForceFactorizer(timeout=timeout)
    start = time.monotonic()
    with pytest.raises(TimeOutError):
        divider.factorize(HARD_N)
    elapsed = time.monotonic() - start
    # 10% margin for OS scheduling jitter
    assert elapsed >= timeout * 0.9, f"fired at {elapsed:.3f}s, expected >={timeout * 0.9}s"


def test_timeout_polling_granularity():
    """Actual wait should not exceed timeout + 100 ms (10 polling ticks)."""
    timeout = 1
    divider = BruteForceFactorizer(timeout=timeout)
    start = time.monotonic()
    with pytest.raises(TimeOutError):
        divider.factorize(HARD_N)
    elapsed = time.monotonic() - start
    assert elapsed < timeout + 0.1, f"polling granularity too coarse: {elapsed:.3f}s"


# ---------------------------------------------------------------------------
# Reuse after timeout
# After TimeOutError the daemon thread is still running in the background.
# The same instance must be safely reusable for a subsequent factorize() call.
# ---------------------------------------------------------------------------

def test_reuse_after_timeout_easy_n():
    divider = BruteForceFactorizer(timeout=1)
    with pytest.raises(TimeOutError):
        divider.factorize(HARD_N)
    facts = divider.factorize(EASY_N)
    assert facts == (3, 19)


def test_reuse_after_timeout_result_correct():
    """Result must not be contaminated by the stale daemon thread's eventual write."""
    divider = BruteForceFactorizer(timeout=1)
    with pytest.raises(TimeOutError):
        divider.factorize(HARD_N)
    # Give the daemon thread a moment to potentially write a stale entry
    time.sleep(0.05)
    d, q = divider.factorize(EASY_N)
    assert d * q == EASY_N


def test_reuse_multiple_times_after_timeout():
    divider = BruteForceFactorizer(timeout=1)
    for _ in range(3):
        with pytest.raises(TimeOutError):
            divider.factorize(HARD_N)
    facts = divider.factorize(EASY_N)
    assert facts == (3, 19)
