from numpy cimport import_array, ndarray
from numpy import empty

cimport simple

cdef class Phaml(object):

    def __init__(self, triangle_files):
        cdef char *s = triangle_files
        cdef int triangle_files_len = len(triangle_files)
        simple.c_phaml_init(s, &triangle_files_len)

    def solve(self):
        simple.c_phaml_solve()

    def get_mesh(self):
        cdef int n, nelem
        simple.c_phaml_get_mesh_info(&n, &nelem)
        cdef ndarray[double, mode="c"] xvert = empty(n, dtype="double")
        cdef ndarray[double, mode="c"] yvert = empty(n, dtype="double")
        cdef ndarray[int, mode="c"] element_vertices = empty(nelem,
                dtype="int32")
        cdef ndarray[int, mode="c"] orders = empty(nelem,
                dtype="int32")
        simple.c_phaml_get_mesh(&n, &xvert[0], &yvert[0],
            &nelem, &element_vertices[0], &orders[0])
        return xvert, yvert, element_vertices, orders
