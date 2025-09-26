"""
Dales data - profile.py

"""
import netCDF4 as nc
import numpy as np

def profile(params):

    """
    Function to gather profile data. Input and output are realised as dictionaries.

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
                         if True timestep1 has to be given, will average profile
                         between timestep and timestep1

    params["timestep1"] : float
                          has to be given as second timestep when averaged is True
                          give in same way as timestep

    Output
    ------

    The only output is a dictionary (here "out")

    out["data"] : 1D numpy array containing the profile data

    out["z_domain"] : 1D numpy array containing vertical levels (z)

    out["unit"] : unit of data

    out["z_domain_unit"] : unit of z_domain 

    out["time"] : simulation time corresponding to timestep given

    out["time1"] : when averaged=True was given gives simulation time of timestep1

    """
    averaged = params["averaged"]
    iexpnr = '{0:03d}'.format(params["iexpnr"])
    timestep = params["timestep"]
    variable = params["variable"]
    sd = params["sim_dir"]

    if sd[-1] != "/":
        sd = sd + "/"

    pf = 'profiles.' + iexpnr + '.nc'
    
    with nc.Dataset(sd+pf,"r",format="NETCDF4") as rootgrp:
        ts = int(timestep*(rootgrp["time"][:].shape[0]-1))
        
        if averaged:
            ts1 = int(params["timestep1"]*(rootgrp["time"][:].shape[0]-1))
            time1 = rootgrp["time"][:][ts1]
        else:
            ts1 = ts+1
        
        profile_data = np.mean(rootgrp[variable][:][ts:ts1,:],axis=0)
        time = rootgrp["time"][:][ts]
        unit = rootgrp[variable].units

        z_var = rootgrp.variables[variable].dimensions[1]
        z_domain = rootgrp[z_var][:]
        z_domain_unit = rootgrp[z_var].units


    profile  = { 
            "data" : profile_data,
            "time" : time,
            "z_domain" : z_domain,
            "z_domain_unit" : z_domain_unit,
            "unit" : unit}

    if averaged :
        profile["time1"] = time1
    
    return profile
