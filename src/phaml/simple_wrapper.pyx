from numpy cimport import_array, ndarray
from numpy import empty

cimport simple

H_UNIFORM   = 1
H_ADAPTIVE  = 2
P_UNIFORM   = 3
P_ADAPTIVE  = 4
HP_ADAPTIVE = 5

HP_BIGGER_ERRIND  =  1
HP_APRIORI        =  2
HP_PRIOR2P_E      =  3
HP_PRIOR2P_H1     =  4
HP_T3S            =  5
HP_ALTERNATE      =  6
HP_TYPEPARAM      =  7
HP_COEF_DECAY     =  8
HP_COEF_ROOT      =  9
HP_SMOOTH_PRED    = 10
HP_NEXT3P         = 11
HP_REFSOLN_EDGE   = 12
HP_REFSOLN_ELEM   = 13
HP_NLP            = 14


HIERARCHICAL_COEFFICIENT = 1
TRUE_DIFF                = 2
LOCAL_PROBLEM_H          = 4
LOCAL_PROBLEM_P          = 5
INITIAL_CONDITION        = 6
EXPLICIT_ERRIND          = 7
EQUILIBRATED_RESIDUAL    = 8
REFSOLN_ERREST           = 9



cdef class Phaml(object):
    cdef int allocated #this is set to 0 by Cython

    def __cinit__(self, triangle_files, int problem_number=1):
        cdef char *s = triangle_files
        cdef int triangle_files_len = len(triangle_files)
        simple.c_phaml_init(s, &triangle_files_len, &problem_number)
        self.allocated = 1

    def __dealloc__(self):
        if self.allocated == 1:
            simple.c_phaml_del()

    def solve(self, params={}):
        cdef double term_energy_err = params.get("term_energy_err", 1e-5)
        cdef int max_eq = params.get("max_eq", 50000)
        cdef int verbose = params.get("verbose", 1)
        cdef int reftype = params.get("reftype", HP_ADAPTIVE)
        cdef int hp_strategy = params.get("hp_strategy", HP_PRIOR2P_H1)
        cdef int derefine = params.get("derefine", 1)
        cdef int degree = params.get("degree", 1)
        cdef int error_estimator = params.get("error_estimator",
                EXPLICIT_ERRIND)
        cdef double lambda0 = params.get("lamda0", 0)
        cdef int lambda_smallest = params.get("lambda_smallest", 1)
        cdef int num_eval = params.get("num_eval", 1)
        simple.c_phaml_solve(&term_energy_err, &max_eq, &verbose, &reftype,
                &hp_strategy, &derefine, &degree, &error_estimator,
                &lambda0, &lambda_smallest, &num_eval)

    def get_mesh(self):
        cdef int n, nelem
        simple.c_phaml_get_mesh_info(&n, &nelem)
        cdef ndarray[double, mode="c"] xvert = empty(n, dtype="double")
        cdef ndarray[double, mode="c"] yvert = empty(n, dtype="double")
        cdef ndarray[int, ndim=2, mode="fortran"] element_vertices = \
                empty((nelem, 3), dtype="int32", order="F")
        cdef ndarray[int, mode="c"] orders = \
                empty(nelem, dtype="int32")
        simple.c_phaml_get_mesh(&n, &xvert[0], &yvert[0],
            &nelem, &element_vertices[0, 0], &orders[0])
        return xvert, yvert, element_vertices, orders

    def get_solution_values(self,
            ndarray[double, mode="c"] x not None,
            ndarray[double, mode="c"] y not None):
        cdef int n = len(x)
        assert len(y) == n
        cdef ndarray[double, mode="c"] values = empty(n, dtype="double")
        simple.c_phaml_get_solution_values(&n, &x[0], &y[0], &values[0])
        return values
