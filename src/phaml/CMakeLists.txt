set(SRC
    iso_c_utilities.f90
    pde_pointers.f90
    example1.f90
    example2.f90
    example_eigen.f90
    simple.f90

    simple_wrapper.c
)

find_package(Python REQUIRED)
find_package(NumPy REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})
include_directories(${PROJECT_BINARY_DIR}/src)

add_python_library(simple_wrapper ${SRC})
target_link_libraries(simple_wrapper phaml)
set_target_properties(simple_wrapper PROPERTIES PREFIX "")

if (WITH_ARPACK)
    target_link_libraries(simple_wrapper arpack parpack)
endif (WITH_ARPACK)
