from numpy cimport import_array, ndarray
from numpy import empty

cimport simple

def run(ndarray[double, mode="c"] x not None,
        ndarray[double, mode="c"] y not None):
    """
    Calls c_run() from simple.f90
    """
    cdef int n = len(x)
    assert len(y) == n
    cdef ndarray[double, mode="c"] sol = empty(n)
    simple.c_run(&n, &x[0], &y[0], &sol[0])
    return sol
