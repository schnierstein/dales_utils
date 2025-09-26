"""

This scripts plots a slice of thl from the rico_roel case found in the cases folder of the dales repository run with 2 MPI processes.

"""


import numpy as np
import matplotlib.pyplot as plt
from field_dump import field_dump

params = {
        "nop" : 2,
        "sim_dir" : "/home/niklas/tests/rico_roel/",
        "variable" : "qt",
        "timestep" : 0.5,
        "iexpnr" : 1,
        "averaged" : False
        }

out = field_dump(params)

height = 30

data = out["data"][height,:,:]
x_domain = out["x_domain"]
y_domain = out["y_domain"]

fig, ax = plt.subplots(1,1)

im = ax.imshow(data,extent = [x_domain[0],x_domain[-1],y_domain[0],y_domain[-1]])

plt.show()

