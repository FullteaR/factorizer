from libcpp.string cimport string

cdef extern from "BruteForceFactorizer_cpp.hpp":
    string BruteForceFactorizer_cppfunc(string s)

cdef extern from "FermatFactorizer_cpp.hpp":
    string FermatFactorizer_cppfunc(string s)