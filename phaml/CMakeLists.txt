set(SRC
    iso_c_utilities.f90
    pde_pointers.f90
    example1.f90
    example2.f90
    example_eigen.f90
    simple.f90

    simple_wrapper.c
)

add_library(simple_wrapper SHARED ${SRC})
target_link_libraries(simple_wrapper phaml)

if (WITH_ARPACK)
    target_link_libraries(simple_wrapper arpack)
endif (WITH_ARPACK)
