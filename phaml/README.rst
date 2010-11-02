Python Wrappers
===============

How it Works
------------

Wrapping Fortran subroutines and functions to Python is a straightforward
process. First we introduce a file::

    simple.f90 ......... Fortran functions exported using the use
                           iso_c_binding Fortran module --- this makes sure
                           that all variable types are the same as in C

and::

    rdirac.pxd   ......... Cython definitions for the exported Fortran functions

Now we write the actual Python wrappers using Cython::

    rdirac_wrapper.pyx ... The actual Python wrappers calling functions from
                           the module rdirac just like they were regular C
                           functions

Then we call Cython::

    $ cython rdirac_wrapper.pyx

which generates::

    rdirac_wrapper.c ..... Generated file (by Cython), we commit it in the
                           repository, so that cython is not needed to compile
                           the Python wrappers

Finally we import things into the module ``rdirac``::

    __init__.py .......... Import things from rdirac_wrapper.so

How to Compile the Wrappers
---------------------------

The files ``c_rdirac.f90`` and ``rdirac_wrapper.c`` have to be compiled and
linked together with ``librdirac.a`` into a ``rdirac_wrapper.so`` Python
extension module.

Cython is not needed to compile the wrappers.

How to Update the Wrappers
--------------------------

Just modify the files above, and rerun Cython (using the command above). You
need to have Cython installed on your system.
