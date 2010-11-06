cdef extern void c_phaml_init(char *triangle_files, int *triangle_files_len,
        int *problem_number)
cdef extern void c_phaml_solve(double *term_energy_err, int *max_eq,
        int *verbose, int *reftype, int *hp_strategy, int *derefine,
        int *degree, int *error_estimator, double *lambda0,
        int *lambda_smallest, int *num_eval)
cdef extern void c_phaml_get_mesh_info(int *n, int *nelem)
cdef extern void c_phaml_get_mesh(int *n, double *xvert, double *yvert,
        int *nelem, int *element_vertices, int *orders)
cdef extern void c_phaml_get_solution_values(int *n, double *x, double *y,
        double *values)
cdef extern void c_phaml_del()
