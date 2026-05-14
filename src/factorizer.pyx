from libcpp.string cimport string
import requests
from requests.exceptions import Timeout as requestsTimeoutError
import threading


class TimeOutError(Exception):
    pass



class BaseClass:

    def __init__(self, timeout=None):
        if timeout is None or timeout<0:
            self.timeout = None
        else:
            self.timeout = timeout

        self.result = {}
        
    def factorize(self, n, *args, **kwargs):
        assert n>=0
        if n==0:
            return (0, 0)
        elif n==1:
            return (1, 1)
        elif n%2==0: # some algorithms can't get correct answer for 2(this is a prime) or even numbers.
            return (2, n//2)
        args = (str(n).encode(),)+args
        thread_factorize = threading.Thread(target=self.factorize_wrap, args=args, kwargs=kwargs, name="_factorize", daemon=True)
        thread_factorize.start()
        thread_factorize.join(timeout=self.timeout)
        if thread_factorize.is_alive():
            raise TimeOutError


        d = self.result[str(n).encode()]
        d = int(d.decode())
        assert d*(n//d) == n
        return (d, n//d)

    
    def factorize_wrap(self, n, *args, **kwargs):
        d = self._factorize(n, *args, **kwargs)
        self.result[n] = d

    
    def _factorize(self, string n, *args, **kwargs):
        return b"1"


class BruteForceFactorizer(BaseClass):

    def __init__(self, timeout=None):
        super().__init__(timeout)

    def _factorize(self, string n, *args, **kwargs):
        cdef:
            string d
        with nogil:
            d = BruteForceFactorizer_cppfunc(n)
        return d


class FermatFactorizer(BaseClass):

    def __init__(self, max_iter=None, timeout=None):
        super().__init__(timeout)
        self.max_iter = max_iter

    def _factorize(self, string n, *args, **kwargs):
        cdef:
            string d
            long long max_iter
        max_iter = self.max_iter if self.max_iter is not None else -1
        with nogil:
            d = FermatFactorizer_cppfunc(n, max_iter)
        return d


class PollardsRhoFactorizer(BaseClass):

    def __init__(self, c=1, timeout=None):
        super().__init__(timeout)
        self.c = c

    def _factorize(self, string n, *args, **kwargs):
        cdef:
            string d
            long c
        c = self.c
        with nogil:
            d = PollardsRhoFactorizer_cppfunc(n, c)
        return d


class PminusOneFactorizer(BaseClass):

    def __init__(self, step=1, timeout=None):
        assert int(step) in (1,2)
        super().__init__(timeout)
        self.step = int(step)
    
    def factorize(self, n, M=100000, L=1000000, *args, **kwargs):
        assert M<ULONG_MAX
        if self.step==2:
            assert L<ULONG_MAX
            assert M<L
        return super().factorize(n, M, L, *args, **kwargs)
    
    def _factorize(self, string n, unsigned long M, unsigned long L, *args, **kwargs):
        cdef:
            string d
        if self.step==1:
            while True:
                with nogil:
                    d = PminusOneFactorizer_step1_cppfunc(n, M)
                if d!=n: #If d!=1 and d!=n, then that is a factor. If d==1, maybe trying bigger M makes good result.
                    return d
                else: #d==n. In this case, mostly M is too big. try smaller M.
                    M = M//2

        else:
            with nogil:
                d = PminusOneFactorizer_step2_cppfunc(n, M, L)
        return d

class RSAPrivateKeyFactorizer(BaseClass):

    def __init__(self, timeout=None):
        super().__init__(timeout)

    def factorize(self, n, d, e=65537, *args, **kwargs):
        return super().factorize(n, str(d).encode(), str(e).encode(), *args, **kwargs)

    def _factorize(self, string n, string d, string e, *args, **kwargs):
        cdef:
            string p
        with nogil:
            p = RSAPrivateKeyFactorizer_cppfunc(n, d, e)
        return p


class PplusOneFactorizer(BaseClass):

    def __init__(self, A=3, timeout=None):
        super().__init__(timeout)
        self.A = A

    def factorize(self, n, M=100000, *args, **kwargs):
        assert M < ULONG_MAX
        return super().factorize(n, M, *args, **kwargs)

    def _factorize(self, string n, unsigned long M, *args, **kwargs):
        cdef:
            string d
            unsigned long A
        A = self.A
        with nogil:
            d = PplusOneFactorizer_step1_cppfunc(n, M, A)
        return d


class ECMFactorizer(BaseClass):

    def __init__(self, timeout=None):
        super().__init__(timeout)

    def factorize(self, n, B1=10000, B2=100000, max_curves=100, *args, **kwargs):
        assert B1 >= 2
        assert B2 >= B1
        assert max_curves >= 1
        return super().factorize(n, B1, B2, max_curves, *args, **kwargs)

    def _factorize(self, string n, unsigned long B1, unsigned long B2, unsigned long max_curves, *args, **kwargs):
        cdef:
            string d
        with nogil:
            d = ECMFactorizer_cppfunc(n, B1, B2, max_curves)
        return d


class FactorDBFactorizer(BaseClass):
    
    def __init__(self, timeout=None):
        super().__init__(timeout)
        self.ENDPOINT = "http://factordb.com/api"

    
    def _factorize(self, n):
        payload = {"query": n}
        
        try:
            r = requests.get(self.ENDPOINT, params=payload, timeout=self.timeout)
            return r.json()
        except requestsTimeoutError:
            raise TimeOutError
        except Exception as e:
            raise e

    def factorize(self, n, raw_result=False):
        if n==0:
            return (0, 0)
        elif n==1:
            return (1, 1)
        result = self._factorize(n)["factors"]
        if raw_result:
            return result
        
        if len(result)==1 and result[0][1]==1:
            return (1, n)
        else:
            d = int(result[0][0])
            assert d*(n//d) == n
            return (d, n//d)
