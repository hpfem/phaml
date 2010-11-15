cmake_minimum_required(VERSION 2.6)

project(phaml)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # for gfortran, we need to use -fPIC
    set(CMAKE_Fortran_FLAGS "-fPIC")
endif ()

set(WITH_ARPACK yes)
set(WITH_PHAML_BLAS_LAPACK no)
set(WITH_PYTHON yes)

if (WITH_ARPACK)
    find_package(ARPACK REQUIRED)
endif (WITH_ARPACK)

add_subdirectory(src)
add_subdirectory(examples)
if (WITH_PYTHON)
    add_subdirectory(phaml)
endif (WITH_PYTHON)

message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message("WITH_ARPACK: ${WITH_ARPACK}")
message("WITH_PHAML_BLAS_LAPACK: ${WITH_PHAML_BLAS_LAPACK}")
message("WITH_PYTHON: ${WITH_PYTHON}")
