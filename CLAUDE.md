# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Build C++/Cython extension in-place (required before running tests or importing)
python setup.py build_ext --inplace

# Editable install for development
pip install -e .
```

Build requires: `libboost-dev` (for Boost.Multiprecision), `cython>=0.29`, `setuptools>=60.5`

## Testing

```bash
pytest tests/                        # Run all tests
pytest tests/test_rho.py             # Run a single test module
pytest tests/ -v                     # Verbose output
```

After any change to `.pyx`, `.pxd`, `.cpp`, or `.hpp` files, rebuild with `python setup.py build_ext --inplace` before running tests.

## Architecture

This project is a Python integer factorization library backed by C++ algorithms, bridged via Cython.

**Layer stack:**
1. **Python interface** — `src/factorizer.pyx` defines Python classes for each algorithm, handles threading/timeout logic, and converts Python `int` ↔ C++ `string`-represented big integers.
2. **Cython declarations** — `src/factorizer.pxd` declares the external C++ function signatures.
3. **C++ implementations** — `src/*_cpp.{hpp,cpp}` contain the actual algorithms; all use `boost::multiprecision::cpp_int` for arbitrary-precision arithmetic.
4. **Utilities** — `src/utils.{hpp,cpp}` provide shared math helpers: `ilog`, `jacobi_symbol`, `next_prime` (Miller-Rabin).

**Algorithms and their best use cases:**
- `BruteForceFactorizer` — trial division, best for small numbers or small factors
- `FermatFactorizer` — difference of squares, best when factors are close in size
- `PollardsRhoFactorizer` — general-purpose; configurable constant `c`
- `PminusOneFactorizer` — Pollard's p-1 with step 1 (bound M) and step 2 (bound L); effective when p-1 is smooth
- `RSAPrivateKeyFactorizer` — recovers primes from RSA (n, e, d)
- `FactorDBFactorizer` — pure Python; queries factordb.com via `requests`

**All factorizers share a `timeout` parameter** (seconds, float). Internally, `BaseClass` in `factorizer.pyx` runs the C++ function in a daemon thread and joins with the timeout.

**Big integer passing:** Python `int` values are converted to `str` before crossing the Cython boundary; C++ receives them as `std::string` and constructs `cpp_int` from them. Results are returned the same way.

## Key Files

| File | Purpose |
|------|---------|
| `src/factorizer.pyx` | All Python-visible classes and timeout logic |
| `src/factorizer.pxd` | Cython extern declarations for C++ functions |
| `src/utils.hpp` / `src/utils.cpp` | Shared number-theory helpers |
| `setup.py` | Builds all `.cpp` + `.pyx` into a single `factorizer` extension (C++14, `-O3`) |
| `tests/` | pytest suite — one file per factorizer |
