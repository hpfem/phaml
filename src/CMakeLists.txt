cmake_minimum_required(VERSION 2.6)

project(phaml)
enable_language(Fortran)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # for gfortran, we need to use -fPIC
    set(CMAKE_Fortran_FLAGS "-fPIC")
endif ()

add_subdirectory(src)

message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")
