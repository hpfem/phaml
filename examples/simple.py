"""
Python version of the example "simple".
"""
import os

from numpy import array

from phaml import Phaml
from phaml.plot import plot_mesh_mpl, plot_sln_mayavi, convert_mesh

problem_number = 1
params = {
        }

current_dir = os.path.dirname(os.path.abspath(__file__))
domain_file = os.path.join(current_dir, "domain")
p = Phaml(domain_file, problem_number)
p.solve(params)
mesh_data = p.get_mesh()
print "Saving the mesh..."
polygons, orders = convert_mesh(*mesh_data)

import matplotlib
matplotlib.use("Agg")
f = plot_mesh_mpl(polygons, orders)
f.savefig("mesh.png")

print "Saving the solution..."
x, y, mesh, _ = mesh_data
values = p.get_solution_values(x, y)

mesh = [elem-1 for elem in mesh]
f = plot_sln_mayavi(x, y, mesh, values)
f.savefig("sln.png")
print "Done"
