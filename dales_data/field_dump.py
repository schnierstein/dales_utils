"""
Dales data - field_dump.py
"""

import netCDF4 as nc
import numpy as np


def field_dump(params):

    """
    Function to gather field dump data. Input and output are realised as dictionaries.

    Input
    ----------

    The only input of the function is a dictionary.


    params["sim_dir"] : string
                        location of the simulation directory

    params["iexpnr"] : int
                       experiment number as specified in the namoptions, i.e.
                       to achieve "001" put in 1

    params["variable"] : string
                         name of the variable to load

    params["nop"] : int
                    number of processes from MPI parallelization

    params["timestep"] : float
                         value between 0 and 1 representing the timestep to be loaded
                         0 represents first value
                         1 represents maximum value of the simulation

    params["averaged"] : bool
                         if True timestep1 has to be given, will average field
                         between timestep and timestep1

    params["timestep1"] : float
                          has to be given as second timestep when averaged is True
                          give in same way as timestep

    Output
    ------

    The only output is a dictionary (here "out")

    out["data"] : 3D numpy array containing the field data

    out["z_domain"] : 1D numpy array containing vertical levels

    out["y_domain"] : 1D numpy array containing horizontal domain (y)

    out["x_domain"] : 1D numpy array containing horizontal domain (x)

    out["unit"] : unit of data

    out["z_domain_unit"] : unit of z_domain 

    out["y_domain_unit"] : unit of y_domain 

    out["x_domain_unit"] : unit of x_domain 

    out["time"] : simulation time corresponding to timestep given

    out["time1"] : when averaged=True was given gives simulation time of timestep1

    """

    averaged = params["averaged"]
    iexpnr = '{0:03d}'.format(params["iexpnr"])
    timestep = params["timestep"]
    variable = params["variable"]
    nop = params["nop"]
    sim_dir = params["sim_dir"]

    for i in range(nop):

        ff = 'fielddump.000.{:03d}.'.format(i) + iexpnr + '.nc'

        with nc.Dataset(sim_dir+ff, "r", format="NETCDF4") as rootgrp:

            if (i==0):
                ts = int(timestep*(rootgrp["time"][:].shape[0]-1))

                if averaged:
                    ts1 = int(params["timestep1"]*(rootgrp["time"][:].shape[0]-1))
                    time = rootgrp["time"][:][ts1]
                else:
                    ts1 = ts+1

                fielddump_data = np.mean(rootgrp[variable][:][ts:ts1,:,:,:],axis=0)
                time = rootgrp["time"][:][ts]
                unit = rootgrp[variable].units

                x_var = rootgrp.variables[variable].dimensions[3]
                y_var = rootgrp.variables[variable].dimensions[2]
                z_var = rootgrp.variables[variable].dimensions[1]


                x_domain = rootgrp[x_var][:]
                y_domain = rootgrp[y_var][:]
                z_domain = rootgrp[z_var][:]

                x_domain_unit = rootgrp[x_var].units
                y_domain_unit = rootgrp[y_var].units
                z_domain_unit = rootgrp[z_var].units

            else:
                
                fielddump_data = \
                        np.hstack( \
                                (fielddump_data, \
                                np.mean(rootgrp[variable][:][ts:ts1,:,:,:],axis=0)))

                y_domain = np.append(y_domain,rootgrp[y_var][:])



    f_dump = {
            "data" : fielddump_data,
            "time" : time,
            "unit" : unit,
            "x_domain" : x_domain,
            "y_domain" : y_domain,
            "z_domain" : z_domain,
            "x_domain_unit" : x_domain_unit,
            "y_domain_unit" : y_domain_unit,
            "z_domain_unit" : z_domain_unit
            }
    if averaged :
        f_dump["time1"] = time1

    return f_dump

