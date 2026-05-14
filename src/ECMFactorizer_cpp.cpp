#include "utils.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
#include <boost/random.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

using Bint = boost::multiprecision::cpp_int;

namespace ecm {

// Montgomery multiplication context (REDC algorithm)
struct MontCtx {
    Bint n;
    Bint R_mask;
    Bint R2;
    Bint n_prime;
    unsigned k;
    Bint MONT_1;

    MontCtx(Bint const& n_) : n(n_) {
        k = (unsigned)msb(n) + 1;
        Bint R = Bint(1) << k;
        R_mask = R - 1;
        Bint R_mod_n = R % n;
        R2 = (R_mod_n * R_mod_n) % n;
        Bint n_inv = mod_inverse(n, R);
        n_prime = R - n_inv;
        MONT_1 = toMont(ONE);
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

    Bint toMont(Bint const& a) const { return REDC(a * R2); }
    Bint fromMont(Bint const& aR) const { return REDC(aR); }
    Bint mul(Bint const& aR, Bint const& bR) const { return REDC(aR * bR); }

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

// Point on Montgomery curve By^2 = x^3 + Ax^2 + x
// Projective (X:Z) where x = X/Z
struct MontPoint {
    Bint X, Z;
};

// xDBL: double a point; a24 = (A+2)/4 in Montgomery domain
static MontPoint xDBL(MontPoint const& P, Bint const& a24, MontCtx const& ctx) {
    Bint s  = ctx.add(P.X, P.Z);
    Bint d  = ctx.sub(P.X, P.Z);
    Bint ss = ctx.mul(s, s);
    Bint dd = ctx.mul(d, d);
    Bint e  = ctx.sub(ss, dd);
    return {ctx.mul(ss, dd), ctx.mul(e, ctx.add(dd, ctx.mul(a24, e)))};
}

// xADD: differential addition P + Q given P-Q
static MontPoint xADD(MontPoint const& P, MontPoint const& Q,
                       MontPoint const& PmQ, MontCtx const& ctx) {
    Bint u    = ctx.mul(ctx.sub(P.X, P.Z), ctx.add(Q.X, Q.Z));
    Bint v    = ctx.mul(ctx.add(P.X, P.Z), ctx.sub(Q.X, Q.Z));
    Bint sum  = ctx.add(u, v);
    Bint diff = ctx.sub(u, v);
    return {ctx.mul(PmQ.Z, ctx.mul(sum, sum)),
            ctx.mul(PmQ.X, ctx.mul(diff, diff))};
}

// Montgomery ladder: compute [k]P
static MontPoint mont_ladder(MontPoint const& P, Bint const& k,
                              Bint const& a24, MontCtx const& ctx) {
    if (k == 0) return {ZERO, ZERO};
    if (k == 1) return P;

    MontPoint R0 = P;
    MontPoint R1 = xDBL(P, a24, ctx);

    int nbits = (int)msb(k);
    for (int i = nbits - 1; i >= 0; i--) {
        if (bit_test(k, i)) {
            R0 = xADD(R0, R1, P, ctx);
            R1 = xDBL(R1, a24, ctx);
        } else {
            R1 = xADD(R0, R1, P, ctx);
            R0 = xDBL(R0, a24, ctx);
        }
    }
    return R0;
}

static Bint modn(Bint const& a, Bint const& n) {
    Bint r = a % n;
    if (r < 0) r += n;
    return r;
}

static unsigned long ul_gcd(unsigned long a, unsigned long b) {
    while (b) { unsigned long t = b; b = a % b; a = t; }
    return a;
}

} // namespace ecm


std::string ECMFactorizer_cppfunc(std::string s, unsigned long B1, unsigned long B2, unsigned long max_curves) {
    using namespace ecm;
    Bint n(s);

    // Python BaseClass handles n=0, n=1, even n.
    // Quick trial division for small odd factors
    for (unsigned long d = 3; d < 1000 && Bint(d) * d <= n; d += 2) {
        if (n % d == 0) return Bint(d).str();
    }
    if (n < 7) return "1";

    MontCtx ctx(n);

    std::random_device rd;
    boost::mt19937 base_gen(rd());
    using engine_type = boost::random::independent_bits_engine<boost::mt19937, 1024, Bint>;
    engine_type gen(base_gen);

    std::vector<unsigned long> primes = generate_primes(B1);

    // BSGS Stage 2 parameters: D = 2*3*5*7 = 210
    constexpr unsigned long D = 210;
    std::vector<unsigned long> baby_residues;
    for (unsigned long r = 1; r <= D / 2; r++) {
        if (ul_gcd(r, D) == 1)
            baby_residues.push_back(r);
    }

    for (unsigned long curve = 0; curve < max_curves; curve++) {
        // Suyama's parameterization: group order divisible by 12
        Bint sigma = modn(gen(), n - 7) + 6;

        Bint u  = modn(sigma * sigma - 5, n);
        Bint v  = modn(4 * sigma, n);
        Bint u3 = modn(u * u * u, n);
        Bint v3 = modn(v * v * v, n);

        // a24 = (A+2)/4 = (v-u)^3 (3u+v) / (16 u^3 v)
        Bint vmu  = modn(v - u, n);
        Bint vmu3 = modn(vmu * vmu * vmu, n);
        Bint num   = modn(vmu3 * modn(3 * u + v, n), n);
        Bint denom = modn(16 * u3 * v, n);

        Bint g = gcd(denom, n);
        if (g > 1) {
            if (g < n) return g.str();
            continue;
        }

        Bint a24_val = modn(num * MontCtx::mod_inverse(denom, n), n);
        Bint a24 = ctx.toMont(a24_val);
        MontPoint P = {ctx.toMont(u3), ctx.toMont(v3)};

        // ── Stage 1 ──
        bool overshot = false;
        Bint accum = ctx.MONT_1;
        unsigned stage1_count = 0;

        for (size_t i = 0; i < primes.size() && primes[i] <= B1; i++) {
            unsigned long p = primes[i];
            unsigned long pp = p;
            while (pp <= B1 / p) pp *= p;

            P = mont_ladder(P, Bint(pp), a24, ctx);
            accum = ctx.mul(accum, P.Z);
            stage1_count++;

            if (stage1_count % 64 == 0) {
                g = gcd(ctx.fromMont(accum), n);
                if (g > 1 && g < n) return g.str();
                if (g == n) { overshot = true; break; }
                accum = ctx.MONT_1;
            }
        }

        if (!overshot) {
            g = gcd(ctx.fromMont(accum), n);
            if (g > 1 && g < n) return g.str();
            if (g == n) overshot = true;
        }
        if (overshot) continue;

        // ── Stage 2: Baby-step Giant-step ──
        if (B2 <= B1 || P.Z == 0) continue;

        // Baby steps: compute [r]P for each r coprime to D, 1 <= r <= D/2
        std::vector<MontPoint> mult(D / 2 + 1);
        mult[1] = P;
        mult[2] = xDBL(P, a24, ctx);
        for (unsigned long i = 3; i <= D / 2; i++)
            mult[i] = xADD(mult[i - 1], P, mult[i - 2], ctx);

        std::vector<MontPoint> baby;
        baby.reserve(baby_residues.size());
        for (unsigned long r : baby_residues)
            baby.push_back(mult[r]);
        mult.clear();

        // Giant step increment S = [D]P
        MontPoint S = mont_ladder(P, Bint(D), a24, ctx);

        unsigned long m_start = (B1 / D) + 1;
        unsigned long m_end   = (B2 / D) + 1;

        // Seed two consecutive giant-step positions
        MontPoint Qprev = mont_ladder(P, Bint((unsigned long long)m_start * D), a24, ctx);
        MontPoint Qcurr;
        bool have_two = (m_start < m_end);
        if (have_two)
            Qcurr = mont_ladder(P, Bint((unsigned long long)(m_start + 1) * D), a24, ctx);

        Bint s2_accum = ctx.MONT_1;
        unsigned long s2_count = 0;
        bool s2_done = false;

        // Accumulate: product of (X_giant * Z_baby - X_baby * Z_giant)
        // A zero delta means x-coords match mod some factor → gcd will reveal it
        auto accumulate = [&](MontPoint const& Q) {
            for (size_t j = 0; j < baby.size(); j++) {
                Bint delta = ctx.sub(ctx.mul(Q.X, baby[j].Z),
                                     ctx.mul(baby[j].X, Q.Z));
                if (delta != 0)
                    s2_accum = ctx.mul(s2_accum, delta);
            }
            s2_count++;
            if (s2_count % 100 == 0) {
                g = gcd(ctx.fromMont(s2_accum), n);
                if (g > 1 && g < n) { s2_done = true; return; }
                if (g == n) { s2_done = true; g = n; return; }
                s2_accum = ctx.MONT_1;
            }
        };

        // Process m_start
        accumulate(Qprev);
        if (s2_done) {
            if (g > 1 && g < n) return g.str();
            continue;
        }

        // Process m_start+1 .. m_end via differential addition
        if (have_two) {
            accumulate(Qcurr);
            if (s2_done) {
                if (g > 1 && g < n) return g.str();
                continue;
            }

            for (unsigned long m = m_start + 2; m <= m_end && !s2_done; m++) {
                MontPoint Qnext = xADD(Qcurr, S, Qprev, ctx);
                Qprev = Qcurr;
                Qcurr = Qnext;
                accumulate(Qcurr);
            }

            if (s2_done) {
                if (g > 1 && g < n) return g.str();
                continue;
            }
        }

        // Final gcd for remaining accumulator
        g = gcd(ctx.fromMont(s2_accum), n);
        if (g > 1 && g < n) return g.str();
    }

    return "1";
}
