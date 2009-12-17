from numpy import array

import simple
x = array([0, 0.6, 1])
y = array([0, 0.5, 1])
#sol = array([0, 0.5, 1])
u = simple.python_gets.run(x, y)
print u
print type(u)
print u.dtype
print simple.python_gets.xvert
print simple.python_gets.yvert
print simple.python_gets.element_vertices[:5, :]
print simple.python_gets.element_order
print simple.python_gets.nvert
print simple.python_gets.nelem

print len(simple.python_gets.element_vertices)
print len(simple.python_gets.element_order)
