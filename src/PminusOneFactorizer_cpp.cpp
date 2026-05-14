#include "utils.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
using Bint = boost::multiprecision::cpp_int;

std::string PminusOneFactorizer_step1_cppfunc(std::string s, unsigned long M){
    Bint n(s);
    Bint A("2");
    if(n%A==0){
        return A.str();
    }
    for(unsigned long p : generate_primes(M)){
        unsigned long f = get_f(p, M);
        Bint Mp = pow(Bint(p), f);
        A = powm(A, Mp, n);
    }
    Bint retval = gcd(A-1, n);
    return retval.str();
}

std::string PminusOneFactorizer_step2_cppfunc(std::string s, unsigned long M, unsigned long L){
    Bint n(s);
    Bint A("2");
    if(n%A==0){
        return A.str();
    }
    for(unsigned long p : generate_primes(L)){
        if(p<M){
            unsigned long f = get_f(p, M);
            Bint Mp = pow(Bint(p), f);
            A = powm(A, Mp, n);
        }
        else{
            Bint AK = powm(A, p, n);
            Bint candidate = gcd(AK-1, n);
            if(1<candidate and candidate<n){
                return candidate.str();
            }
        }
    }
    return "1";
}