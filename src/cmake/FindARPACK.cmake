FIND_LIBRARY(ARPACK_LIBRARY libarpack.a libarpack.so)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ARPACK DEFAULT_MSG ARPACK_LIBRARY)