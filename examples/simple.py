"""
Python version of the example "simple".
"""
import os

from numpy import array

from phaml import simple

x = array([0, 0.6, 1])
y = array([0, 0.5, 1])
#sol = array([0, 0.5, 1])
current_dir = os.path.dirname(os.path.abspath(__file__))
domain_file = os.path.join(current_dir, "domain")
u = simple.python_gets.run(x, y, triangle_files=domain_file)
print u
print type(u)
print u.dtype
print simple.python_gets.xvert[:100]
print simple.python_gets.yvert
print simple.python_gets.element_vertices
print simple.python_gets.element_order
print simple.python_gets.nvert
print simple.python_gets.nelem

print len(simple.python_gets.element_vertices)
print len(simple.python_gets.element_order)
