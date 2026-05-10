#include "utils.hpp"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/integer.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <cmath>
#include <cassert>

using Bint = boost::multiprecision::cpp_int;


Bint ilog(Bint n){
    Bint result = ZERO;
    while(n>0){
        n = n/2;
        result++;
    }
    return result;
}

Bint ilog(Bint n, Bint const& a){
    Bint result = ZERO;
    while(n>0){
        n = n/a;
        result++;
    }
    return result;
}

Bint _jacobi_symbol(Bint const& a, Bint const& p){
    if(a==1){ // (11) p43 [Wada 2001]
        return ONE;
    }
    else if(a==-1){ // (14) p43 [Wada 2001]
        if(a%4==1){
            return ONE;
        }
        else{
            return MINUS_ONE;
        }
    }
    else if(a==2){ // (15) p43 [Wada 2001]// (15) p43 [Wada 2001]
        int tmp = (int)(a%8);
        if(tmp == 1 or tmp == 7){
            return ONE;
        }
        else{
            return MINUS_ONE;
        }
    }

    if(a%2==0){ // (13) p43 [Wada 2001]
        return _jacobi_symbol(2, p) * _jacobi_symbol(a/2, p);
    }
    else{ // (16) p43 [Wada 2001]
        if(a%4==1 or p%4==1){
            return _jacobi_symbol(p%a, a);
        }
        else{
            return MINUS_ONE * _jacobi_symbol(p%a, a);
        }
    }
}


Bint jacobi_symbol(Bint const& a, Bint const& p){
    //returns (a/p)
    //https://en.wikipedia.org/wiki/Jacobi_symbol
    assert(p%2==1);
    if(a == 0){
        return ZERO;
    }
    assert(gcd(a,p)==1);
    return _jacobi_symbol(a%p, p);
}


//return first(smallest) prime p s.t. p>=n
Bint next_prime(Bint const& n){
    if(n<=2){
        return TWO;
    }
    Bint retval = n;
    if(retval%2==1){
        retval++;
    }
    while(true){
        if(boost::multiprecision::miller_rabin_test(retval, 25)){
            return retval;
        }
        retval+=2;
    }
}