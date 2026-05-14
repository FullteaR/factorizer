from libcpp.string cimport string

cdef extern from "limits.h":
    cdef int INT_MAX
    cdef unsigned int UINT_MAX
    cdef unsigned long ULONG_MAX

cdef extern from "BruteForceFactorizer_cpp.hpp":
    string BruteForceFactorizer_cppfunc(string s) nogil

cdef extern from "FermatFactorizer_cpp.hpp":
    string FermatFactorizer_cppfunc(string s) nogil

cdef extern from "PollardsRhoFactorizer_cpp.hpp":
    string PollardsRhoFactorizer_cppfunc(string s, long c) nogil

cdef extern from "PminusOneFactorizer_cpp.hpp":
    string PminusOneFactorizer_step1_cppfunc(string s, unsigned long M) nogil

cdef extern from "PminusOneFactorizer_cpp.hpp":
    string PminusOneFactorizer_step2_cppfunc(string s, unsigned long M, unsigned long L) nogil

cdef extern from "RSAPrivateKeyFactorizer_cpp.hpp":
    string RSAPrivateKeyFactorizer_cppfunc(string s, string d, string e) nogil

cdef extern from "ECMFactorizer_cpp.hpp":
    string ECMFactorizer_cppfunc(string s, unsigned long B1, unsigned long B2, unsigned long max_curves) nogil