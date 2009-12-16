from numpy import array

import simple
x = array([0, 0.6, 1])
y = array([0, 0.5, 1])
#sol = array([0, 0.5, 1])
u = simple.run(x, y)
print u
print type(u)
print u.dtype
