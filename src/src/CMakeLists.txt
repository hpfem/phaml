set(SRC
    dbinom.f
    basis.f90             hypretype.f90     petsc_interf.dum.f90
    blopex_c.dum.f90      petsctype.dum.f90
    krylov.f90            phaml.f90
    cpusec_f95.f90        lapack_solve.f90      phamltype.f90
    eigen.f90             linsys.f90            pvmpack1.f90
    errest.f90            linsys_io.f90         pvmpack2.f90
    evaluate.f90          linsystype.f90        quadrules.f90
    global.f90            refine_adapt.f90      calgo582.f
    graphics.dum.f90      loadbal.f90           refine_elements.f90
    refine_top.f90        templates.f
    graphmain.f90         messpass.dummy.f90    refine_uniform.f90
    grid.f90
    grid_init.f90         slave.f90
    grid_io.f90           stopwatch.f90
    gridtype.f90          superlu_fwrap.f90
    grid_util.f90         mpi2_stringf.f90      superlu_interf.dum.f90
    hash1.f90             mpipack1.f90          superlu_interf.f90
    mpipack2.f90          superlutype.f90
    hashext1.f90          mpi_stringf.f90       sysdep.f90
    #hashext2.f90
    hbmg.f90
    hp_strategies.f90
    hypre_interf.dum.f90  openmp_stubs.f90      zoltanP_r_f_dum.f90
    zoltan_interf.dum.f90
    hypre_interf.f90      petsc_init.dum.f90
    mumps_interf.dum.f90  mumps_struc_dum.f90
    sort.f
    linsys_util.f90
    make_linsys.f90
    )

if (NOT WITH_ARPACK)
    set(SRC ${SRC} arpack.dum.f)
endif (NOT WITH_ARPACK)
# Always use dummy parpack for now:
set(SRC ${SRC} parpack.dum.f)

if (WITH_PHAML_BLAS_LAPACK)
    set(SRC ${SRC} blas.f lapack.f)
endif (WITH_PHAML_BLAS_LAPACK)

add_library(phaml ${SRC})
