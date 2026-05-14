#include "utils.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
using Bint = boost::multiprecision::cpp_int;

// Compute V_k(P, 1) mod n using the binary ladder:
//   V_{2m}   = V_m^2 - 2
//   V_{2m+1} = V_m * V_{m+1} - P
static Bint lucas_V(Bint P, unsigned long k, Bint const& n) {
    if (k == 0) return Bint(2);
    P = ((P % n) + n) % n;
    if (k == 1) return P;

    int bits = 0;
    for (unsigned long tmp = k; tmp; tmp >>= 1) bits++;

    Bint Vm  = P;
    Bint Vm1 = ((P * P - 2) % n + n) % n;

    for (int i = bits - 2; i >= 0; i--) {
        if ((k >> i) & 1) {
            Bint new_Vm  = ((Vm * Vm1 - P) % n + n) % n;
            Bint new_Vm1 = ((Vm1 * Vm1 - 2) % n + n) % n;
            Vm  = new_Vm;
            Vm1 = new_Vm1;
        } else {
            Bint new_Vm1 = ((Vm * Vm1 - P) % n + n) % n;
            Bint new_Vm  = ((Vm * Vm  - 2) % n + n) % n;
            Vm  = new_Vm;
            Vm1 = new_Vm1;
        }
    }
    return Vm;
}

std::string PplusOneFactorizer_step1_cppfunc(std::string s, unsigned long M, unsigned long A) {
    Bint n(s);
    if (n % 2 == 0) return "2";
    Bint v(A);
    for (unsigned long p : generate_primes(M)) {
        unsigned long f  = get_f(p, M);
        unsigned long pe = 1;
        for (unsigned long i = 0; i < f; i++) pe *= p;
        v = lucas_V(v, pe, n);
    }
    Bint d = gcd(v - 2, n);
    return d.str();
}
