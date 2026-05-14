#include "utils.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
#include <boost/random.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>

using Bint = boost::multiprecision::cpp_int;

namespace ecm {

static Bint modn(Bint const& a, Bint const& n) {
    Bint r = a % n;
    if (r < 0) r += n;
    return r;
}

static Bint mulmod(Bint const& a, Bint const& b, Bint const& n) {
    return modn(a * b, n);
}

static Bint addmod(Bint const& a, Bint const& b, Bint const& n) {
    return modn(a + b, n);
}

static Bint submod(Bint const& a, Bint const& b, Bint const& n) {
    return modn(a - b, n);
}

struct Point {
    Bint X, Y, Z;
};

static Point infinity() {
    return {ZERO, ONE, ZERO};
}

static bool is_inf(Point const& P) {
    return P.Z == 0;
}

// Point doubling in Jacobian coordinates for y^2 = x^3 + ax + b
static Point ec_double(Point const& P, Bint const& a, Bint const& n) {
    if (is_inf(P) || P.Y == 0) return infinity();

    Bint Y2 = mulmod(P.Y, P.Y, n);
    Bint S = mulmod(Bint(4), mulmod(P.X, Y2, n), n);
    Bint Z2 = mulmod(P.Z, P.Z, n);
    Bint Z4 = mulmod(Z2, Z2, n);
    Bint M = addmod(mulmod(Bint(3), mulmod(P.X, P.X, n), n),
                    mulmod(a, Z4, n), n);
    Bint X3 = submod(mulmod(M, M, n), mulmod(Bint(2), S, n), n);
    Bint Y4 = mulmod(Y2, Y2, n);
    Bint Y3 = submod(mulmod(M, submod(S, X3, n), n),
                     mulmod(Bint(8), Y4, n), n);
    Bint Z3 = mulmod(Bint(2), mulmod(P.Y, P.Z, n), n);

    return {X3, Y3, Z3};
}

// Point addition in Jacobian coordinates
static Point ec_add(Point const& P, Point const& Q, Bint const& a, Bint const& n) {
    if (is_inf(P)) return Q;
    if (is_inf(Q)) return P;

    Bint Z1_2 = mulmod(P.Z, P.Z, n);
    Bint Z2_2 = mulmod(Q.Z, Q.Z, n);
    Bint U1 = mulmod(P.X, Z2_2, n);
    Bint U2 = mulmod(Q.X, Z1_2, n);
    Bint S1 = mulmod(P.Y, mulmod(Q.Z, Z2_2, n), n);
    Bint S2 = mulmod(Q.Y, mulmod(P.Z, Z1_2, n), n);

    Bint H = submod(U2, U1, n);
    Bint R = submod(S2, S1, n);

    if (H == 0) {
        if (R == 0) return ec_double(P, a, n);
        return infinity();
    }

    Bint H2 = mulmod(H, H, n);
    Bint H3 = mulmod(H, H2, n);
    Bint U1H2 = mulmod(U1, H2, n);
    Bint X3 = submod(submod(mulmod(R, R, n), H3, n),
                     mulmod(Bint(2), U1H2, n), n);
    Bint Y3 = submod(mulmod(R, submod(U1H2, X3, n), n),
                     mulmod(S1, H3, n), n);
    Bint Z3 = mulmod(H, mulmod(P.Z, Q.Z, n), n);

    return {X3, Y3, Z3};
}

// Scalar multiplication via double-and-add
static Point ec_mul(Point P, Bint k, Bint const& a, Bint const& n) {
    if (k == 0 || is_inf(P)) return infinity();
    if (k < 0) {
        k = -k;
        P.Y = modn(-P.Y, n);
    }

    Point result = infinity();

    while (k > 0) {
        if (k % 2 != 0) {
            result = ec_add(result, P, a, n);
        }
        P = ec_double(P, a, n);
        k >>= 1;
    }

    return result;
}

static std::vector<unsigned long> sieve_primes(unsigned long MAX) {
    if (MAX < 2) return {};
    std::vector<bool> v(MAX + 1, true);
    v[0] = v[1] = false;
    unsigned long sq = (unsigned long)std::sqrt((double)MAX);
    std::vector<unsigned long> primes;
    for (unsigned long i = 2; i <= MAX; i++) {
        if (!v[i]) continue;
        primes.push_back(i);
        if (i <= sq) {
            for (unsigned long j = i * i; j <= MAX; j += i) {
                v[j] = false;
            }
        }
    }
    return primes;
}

} // namespace ecm


std::string ECMFactorizer_cppfunc(std::string s, unsigned long B1, unsigned long B2, unsigned long max_curves) {
    using namespace ecm;
    Bint n(s);

    using engine_type = boost::random::independent_bits_engine<boost::mt19937, 1024, Bint>;
    engine_type gen;

    std::vector<unsigned long> primes = sieve_primes(B2);

    for (unsigned long curve = 0; curve < max_curves; curve++) {
        Bint x0 = modn(gen(), n);
        Bint y0 = modn(gen(), n);
        Bint a = modn(gen(), n);

        Point P = {x0, y0, ONE};

        // Stage 1: multiply P by each prime power up to B1
        bool overshot = false;
        bool found = false;
        for (size_t i = 0; i < primes.size() && primes[i] <= B1; i++) {
            unsigned long p = primes[i];
            unsigned long pp = p;
            while (pp <= B1 / p) pp *= p;

            P = ec_mul(P, Bint(pp), a, n);

            if (is_inf(P)) {
                overshot = true;
                break;
            }

            Bint g = gcd(P.Z, n);
            if (g > 1 && g < n) return g.str();
            if (g == n) { overshot = true; break; }
        }

        if (overshot) continue;

        Bint g = gcd(P.Z, n);
        if (g > 1 && g < n) return g.str();
        if (g == n) continue;

        // Stage 2: step through primes in (B1, B2] using precomputed gaps
        if (B2 > B1 && !is_inf(P)) {
            auto it = std::lower_bound(primes.begin(), primes.end(),
                                       (unsigned long)(B1 + 1));
            if (it != primes.end() && *it <= B2) {
                std::map<unsigned long, Point> gap_cache;

                Point Q = ec_mul(P, Bint(*it), a, n);
                if (is_inf(Q)) continue;
                Bint accum = Q.Z;

                unsigned long prev = *it;
                ++it;
                unsigned long count = 0;
                bool stage2_done = false;

                for (; it != primes.end() && *it <= B2; ++it) {
                    unsigned long gap = *it - prev;
                    if (gap_cache.find(gap) == gap_cache.end()) {
                        gap_cache[gap] = ec_mul(P, Bint(gap), a, n);
                    }
                    Q = ec_add(Q, gap_cache[gap], a, n);
                    if (is_inf(Q)) break;
                    accum = mulmod(accum, Q.Z, n);
                    prev = *it;
                    count++;

                    if (count % 100 == 0) {
                        g = gcd(accum, n);
                        if (g > 1 && g < n) return g.str();
                        if (g == n) { stage2_done = true; break; }
                        accum = ONE;
                    }
                }

                if (!stage2_done) {
                    g = gcd(accum, n);
                    if (g > 1 && g < n) return g.str();
                }
            }
        }
    }

    return "1";
}
