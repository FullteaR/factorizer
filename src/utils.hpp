#ifndef UTILS_HPP_FACTORIZER_4e17c669_cf51_4c34_8ac1_ae3306060740
#define UTILS_HPP_FACTORIZER_4e17c669_cf51_4c34_8ac1_ae3306060740

#include <boost/multiprecision/cpp_int.hpp>
#include <vector>
using Bint = boost::multiprecision::cpp_int;

Bint ilog(Bint n);
Bint ilog(Bint n, Bint const& a);

//returns (a/p)
//https://en.wikipedia.org/wiki/Jacobi_symbol
Bint jacobi_symbol(Bint const& a, Bint const& p);

//return first(smallest) prime p s.t. p>=n
Bint next_prime(Bint n);

unsigned long get_f(unsigned long q, unsigned long M);
std::vector<unsigned long> generate_primes(unsigned long MAX);

#define ZERO Bint("0")
#define ONE Bint("1")
#define TWO Bint("2")
#define MINUS_ONE Bint("-1")


#endif // UTILS_HPP_FACTORIZER_4e17c669_cf51_4c34_8ac1_ae3306060740