"""

This scripts plots a basic profile of thl from the rico_roel case found in the cases folder of the dales repository run with 2 MPI processes.

"""


import numpy as np
import matplotlib.pyplot as plt
from profile import profile


params = {
        "nop" : 2,
        "sim_dir" : "/home/niklas/tests/rico_roel",
        "variable" : "thl",
        "timestep" : 0.5,
        "iexpnr" : 1,
        "averaged" : False
        }

out = profile(params)

max_height = 30

data = out["data"][:max_height]
z_domain = out["z_domain"][:max_height]

fig, ax = plt.subplots(1,1)

ax.plot(data,z_domain,label="Profile after {} seconds".format(out["time"]))

ax.set_ylabel("Height in [{}]".format(out["z_domain_unit"]))
ax.set_xlabel("{} in [{}]".format(params["variable"],out["unit"]))

ax.legend()

plt.show()

