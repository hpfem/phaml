"""
Python version of the example "simple".
"""
import os

from numpy import array

from phaml import Phaml

current_dir = os.path.dirname(os.path.abspath(__file__))
domain_file = os.path.join(current_dir, "domain")
p = Phaml(domain_file)
p.solve()
x, y, elems, orders = p.get_mesh()
for n, elem in enumerate(elems):
    print "Element #%d: %s, order=%d" % (n, elem, orders[n])
    for i in elem:
        print "   ", x[i], y[i]
values = p.get_solution_values(x, y)
print values[:10]
