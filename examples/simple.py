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
print p.get_mesh()
