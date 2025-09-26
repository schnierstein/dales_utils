"""

This scripts plots a slice of thl from the rico_roel case found in the cases folder of the dales repository run with 2 MPI processes.

"""


import numpy as np
import matplotlib.pyplot as plt
from crosssection import crosssection

params = {
        "nop" : 2,
        "sim_dir" : "/home/niklas/tests/rico_roel/",
        "variable" : "thlxz",
        "timestep" : 0.5,
        "direction" : 'xz',
        "iexpnr" : 1,
        "averaged" : False
        }

out = crosssection(params)


data = out["data"]
x_domain = out["x_domain"]
z_domain = out["z_domain"]

fig, ax = plt.subplots(1,1)

im = ax.imshow(data,extent = [x_domain[0],x_domain[-1],z_domain[0],z_domain[-1]])

plt.show()

