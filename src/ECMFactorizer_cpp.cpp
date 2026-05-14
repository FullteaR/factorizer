#include "utils.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
#include <boost/random.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <random>

using Bint = boost::multiprecision::cpp_int;

namespace ecm {

// Montgomery multiplication context
struct MontCtx {
    Bint n;
    Bint R_mask;  // R - 1 where R = 2^k
    Bint R2;      // R^2 mod n (for toMont)
    Bint n_prime; // -n^(-1) mod R
    unsigned k;   // bit width, R = 2^k

    Bint MONT_1, MONT_2, MONT_3, MONT_4, MONT_8;

    MontCtx(Bint const& n_) : n(n_) {
        k = (unsigned)msb(n) + 1;
        Bint R = Bint(1) << k;
        R_mask = R - 1;

        Bint R_mod_n = R % n;
        R2 = (R_mod_n * R_mod_n) % n;

        // n^(-1) mod R via iterative extended GCD
        Bint n_inv = mod_inverse(n, R);
        n_prime = R - n_inv;

        MONT_1 = toMont(ONE);
        MONT_2 = toMont(TWO);
        MONT_3 = toMont(Bint(3));
        MONT_4 = toMont(Bint(4));
        MONT_8 = toMont(Bint(8));
    }

    static Bint mod_inverse(Bint a, Bint const& m) {
        Bint m0 = m, x0 = 0, x1 = 1;
        if (m == 1) return ZERO;
        while (a > 1) {
            Bint q = a / m0;
            Bint t = m0;
            m0 = a % m0;
            a = t;
            t = x0;
            x0 = x1 - q * x0;
            x1 = t;
        }
        if (x1 < 0) x1 += m;
        return x1;
    }

    Bint REDC(Bint const& T) const {
        Bint m = ((T & R_mask) * n_prime) & R_mask;
        Bint t = (T + m * n) >> k;
        if (t >= n) t -= n;
        return t;
    }

    Bint toMont(Bint const& a) const {
        return REDC(a * R2);
    }

    Bint fromMont(Bint const& aR) const {
        return REDC(aR);
    }

    Bint mul(Bint const& aR, Bint const& bR) const {
        return REDC(aR * bR);
    }

    Bint add(Bint const& aR, Bint const& bR) const {
        Bint r = aR + bR;
        if (r >= n) r -= n;
        return r;
    }

    Bint sub(Bint const& aR, Bint const& bR) const {
        Bint r = aR - bR;
        if (r < 0) r += n;
        return r;
    }
};

struct Point {
    Bint X, Y, Z;
};

static Point infinity(MontCtx const& ctx) {
    return {ZERO, ctx.MONT_1, ZERO};
}

static bool is_inf(Point const& P) {
    return P.Z == 0;
}

// Jacobian doubling: y^2 = x^3 + ax + b
static Point ec_double(Point const& P, Bint const& a_mont, MontCtx const& ctx) {
    if (is_inf(P) || P.Y == 0) return infinity(ctx);

    Bint Y2 = ctx.mul(P.Y, P.Y);
    Bint S  = ctx.mul(ctx.MONT_4, ctx.mul(P.X, Y2));
    Bint Z2 = ctx.mul(P.Z, P.Z);
    Bint Z4 = ctx.mul(Z2, Z2);
    Bint M  = ctx.add(ctx.mul(ctx.MONT_3, ctx.mul(P.X, P.X)),
                      ctx.mul(a_mont, Z4));
    Bint X3 = ctx.sub(ctx.mul(M, M), ctx.mul(ctx.MONT_2, S));
    Bint Y4 = ctx.mul(Y2, Y2);
    Bint Y3 = ctx.sub(ctx.mul(M, ctx.sub(S, X3)),
                      ctx.mul(ctx.MONT_8, Y4));
    Bint Z3 = ctx.mul(ctx.MONT_2, ctx.mul(P.Y, P.Z));

    return {X3, Y3, Z3};
}

// Jacobian addition
static Point ec_add(Point const& P, Point const& Q, Bint const& a_mont, MontCtx const& ctx) {
    if (is_inf(P)) return Q;
    if (is_inf(Q)) return P;

    Bint Z1_2 = ctx.mul(P.Z, P.Z);
    Bint Z2_2 = ctx.mul(Q.Z, Q.Z);
    Bint U1   = ctx.mul(P.X, Z2_2);
    Bint U2   = ctx.mul(Q.X, Z1_2);
    Bint S1   = ctx.mul(P.Y, ctx.mul(Q.Z, Z2_2));
    Bint S2   = ctx.mul(Q.Y, ctx.mul(P.Z, Z1_2));

    Bint H = ctx.sub(U2, U1);
    Bint R = ctx.sub(S2, S1);

    if (H == 0) {
        if (R == 0) return ec_double(P, a_mont, ctx);
        return infinity(ctx);
    }

    Bint H2   = ctx.mul(H, H);
    Bint H3   = ctx.mul(H, H2);
    Bint U1H2 = ctx.mul(U1, H2);
    Bint X3   = ctx.sub(ctx.sub(ctx.mul(R, R), H3),
                        ctx.mul(ctx.MONT_2, U1H2));
    Bint Y3   = ctx.sub(ctx.mul(R, ctx.sub(U1H2, X3)),
                        ctx.mul(S1, H3));
    Bint Z3   = ctx.mul(H, ctx.mul(P.Z, Q.Z));

    return {X3, Y3, Z3};
}

// NAF (Non-Adjacent Form) representation
static std::vector<int8_t> compute_naf(Bint k) {
    std::vector<int8_t> naf;
    while (k > 0) {
        if (bit_test(k, 0)) {
            int rem = (int)(k % 4);
            int8_t ki = (int8_t)(2 - rem);
            naf.push_back(ki);
            k -= ki;
        } else {
            naf.push_back(0);
        }
        k >>= 1;
    }
    return naf;
}

// Scalar multiplication using NAF
static Point ec_mul(Point const& P, Bint const& k, Bint const& a_mont, MontCtx const& ctx) {
    if (k == 0 || is_inf(P)) return infinity(ctx);

    auto naf = compute_naf(k);

    Point negP = {P.X, ctx.sub(ZERO, P.Y), P.Z};

    Point result = infinity(ctx);
    for (int i = (int)naf.size() - 1; i >= 0; i--) {
        result = ec_double(result, a_mont, ctx);
        if (naf[i] == 1) {
            result = ec_add(result, P, a_mont, ctx);
        } else if (naf[i] == -1) {
            result = ec_add(result, negP, a_mont, ctx);
        }
    }
    return result;
}

static Bint modn(Bint const& a, Bint const& n) {
    Bint r = a % n;
    if (r < 0) r += n;
    return r;
}

static bool is_singular(Bint const& a, Bint const& b, Bint const& n) {
    Bint disc = modn(Bint(4) * a * a * a + Bint(27) * b * b, n);
    return gcd(disc, n) == n;
}

} // namespace ecm


std::string ECMFactorizer_cppfunc(std::string s, unsigned long B1, unsigned long B2, unsigned long max_curves) {
    using namespace ecm;
    Bint n(s);

    MontCtx ctx(n);

    std::random_device rd;
    boost::mt19937 base_gen(rd());
    using engine_type = boost::random::independent_bits_engine<boost::mt19937, 1024, Bint>;
    engine_type gen(base_gen);

    std::vector<unsigned long> primes = generate_primes(B2);

    for (unsigned long curve = 0; curve < max_curves; curve++) {
        Bint x0 = modn(gen(), n);
        Bint y0 = modn(gen(), n);
        Bint a  = modn(gen(), n);
        Bint b  = modn(y0 * y0 - x0 * x0 * x0 - a * x0, n);

        if (is_singular(a, b, n)) continue;

        Bint a_mont = ctx.toMont(a);
        Point P = {ctx.toMont(x0), ctx.toMont(y0), ctx.MONT_1};

        // Stage 1
        bool overshot = false;
        for (size_t i = 0; i < primes.size() && primes[i] <= B1; i++) {
            unsigned long p = primes[i];
            unsigned long pp = p;
            while (pp <= B1 / p) pp *= p;

            P = ec_mul(P, Bint(pp), a_mont, ctx);

            if (is_inf(P)) {
                overshot = true;
                break;
            }

            Bint z_real = ctx.fromMont(P.Z);
            Bint g = gcd(z_real, n);
            if (g > 1 && g < n) return g.str();
            if (g == n) { overshot = true; break; }
        }

        if (overshot) continue;

        {
            Bint z_real = ctx.fromMont(P.Z);
            Bint g = gcd(z_real, n);
            if (g > 1 && g < n) return g.str();
            if (g == n) continue;
        }

        // Stage 2
        if (B2 > B1 && !is_inf(P)) {
            auto it = std::lower_bound(primes.begin(), primes.end(),
                                       (unsigned long)(B1 + 1));
            if (it != primes.end() && *it <= B2) {
                std::map<unsigned long, Point> gap_cache;

                Point Q = ec_mul(P, Bint(*it), a_mont, ctx);
                if (is_inf(Q)) continue;
                Bint accum = ctx.fromMont(Q.Z);

                unsigned long prev = *it;
                ++it;
                unsigned long count = 0;
                bool stage2_done = false;
                Bint g;

                for (; it != primes.end() && *it <= B2; ++it) {
                    unsigned long gap = *it - prev;
                    if (gap_cache.find(gap) == gap_cache.end()) {
                        gap_cache[gap] = ec_mul(P, Bint(gap), a_mont, ctx);
                    }
                    Q = ec_add(Q, gap_cache[gap], a_mont, ctx);
                    if (is_inf(Q)) break;
                    accum = (accum * ctx.fromMont(Q.Z)) % n;
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
