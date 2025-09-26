"""
Dales data - crosssection.py
"""

import netCDF4 as nc
import numpy as np

def crosssection(params):

    """
    Function to gather crosssection data. Input and output are realised as dictionaries.

    Input
    ----------

    The only input of the function is a dictionary.


    params["sim_dir"] : string
                        location of the simulation directory

    params["iexpnr"] : int
                       experiment number as specified in the namoptions, i.e.
                       to achieve "001" put in 1

    params["nop"] : int
                    number of processes from MPI parallelization

    params["variable"] : string
                         name of the variable to load

    params["direction"] : string
                          Decide which crosssection to load
                          Options: 'xy','xz','yz'

    params["timestep"] : float
                         value between 0 and 1 representing the timestep to be loaded
                         0 represents first value
                         1 represents maximum value of the simulation

    params["averaged"] : bool
                         if True timestep1 has to be given, will average crosssection
                         between timestep and timestep1

    params["timestep1"] : float
                          has to be given as second timestep when averaged is True
                          give in same way as timestep

    Output
    ------

    The only output is a dictionary (here "out")

    out["data"] : 2D numpy array containing the crosssection data

    out["z_domain"] : 1D numpy array containing vertical levels, if 'z' in direction

    out["y_domain"] : 1D numpy array containing horizontal domain (y), if 'y' in direction

    out["x_domain"] : 1D numpy array containing horizontal domain (x), if 'x' in direction

    out["unit"] : unit of data

    out["z_domain_unit"] : unit of z_domain, if 'z' in direction

    out["y_domain_unit"] : unit of y_domain, if 'y' in direction

    out["x_domain_unit"] : unit of x_domain, if 'x' in direction 

    out["time"] : simulation time corresponding to timestep given

    out["time1"] : when averaged=True was given gives simulation time of timestep1

    """

    averaged = params["averaged"]
    #time_dim = params["time_dim"]
    time_dim = False
    iexpnr = '{0:03d}'.format(params["iexpnr"])
    direction = params["direction"]
    timestep = params["timestep"]
    variable = params["variable"]
    sd = params["sim_dir"]
    nop = params["nop"]



    if (averaged and time_dim):
        print("ERROR: Can't have averaged and time_dim, untested behavior")
        return 1

    if direction == 'xz':

        cf = 'crossxz.x000y000.' + iexpnr + '.nc'

        with nc.Dataset(sd+cf,"r",format="NETCDF4") as rootgrp:
            ts = int(timestep*(rootgrp["time"][:].shape[0]-1))
            
            if averaged:
                ts1 = int(params["timestep1"]*(rootgrp["time"][:].shape[0]-1))
                time1 = rootgrp["time"][:][ts1]
            else:
                ts1 = ts+1
            
            if time_dim:
                cross_data = rootgrp[variable][:][ts:ts1,:,:]
            else:
                cross_data = np.mean(rootgrp[variable][:][ts:ts1,:,:],axis=0)
            
            time = rootgrp["time"][:][ts]
            unit = rootgrp[variable].units

            x_var = rootgrp.variables[variable].dimensions[2]
            z_var = rootgrp.variables[variable].dimensions[1]

            
            x_domain = rootgrp[x_var][:]
            z_domain = rootgrp[z_var][:]

            x_domain_unit = rootgrp[x_var].units
            z_domain_unit = rootgrp[z_var].units
            

    elif direction == 'xy':

        for i in range(nop):
            #TODO: check why 0002
            cf = 'crossxy.0002.x000y{:03d}.'.format(i) + iexpnr + '.nc'
            with nc.Dataset(sd+cf,"r",format="NETCDF4") as rootgrp:

                if (i==0):
                    ts = int(timestep*(rootgrp["time"][:].shape[0]-1))
                    if averaged:
                        ts1 = int(params["timestep1"]*(rootgrp["time"][:].shape[0]-1))
                        time1 = rootgrp["time"][:][ts1]
                    else:
                        ts1 = ts+1
                    cross_data = np.mean(rootgrp[variable][:][ts:ts1,:,:],axis=0)
                    time = rootgrp["time"][:][ts]
                    unit = rootgrp[variable].units

                    x_var = rootgrp.variables[variable].dimensions[2]
                    y_var = rootgrp.variables[variable].dimensions[1]

                    x_domain = rootgrp[x_var][:]
                    y_domain = rootgrp[y_var][:]

                    x_domain_unit = rootgrp[x_var].units
                    y_domain_unit = rootgrp[y_var].units

                else:
                    cross_data = \
                            np.vstack(
                                    (cross_data,np.mean(rootgrp[variable][:][ts:ts1,:,:],axis=0)))

                    y_domain = np.append(y_domain,rootgrp[y_var][:])


    elif direction == 'yz':

        for i in range(nop):
            cf = 'crossyz.x000y{:03d}.'.format(i) + iexpnr + '.nc'
            with nc.Dataset(sd+cf,"r",format="NETCDF4") as rootgrp:

                if (i==0):
                    ts = int(timestep*rootgrp["time"][:].shape[0])
                    if averaged:
                        ts1 = int(params["timestep1"]*(rootgrp["time"][:].shape[0]-1))
                        time1 = rootgrp["time"][:][ts1]
                    else:
                        ts1 = ts+1
                    cross_data = np.mean(rootgrp[variable][:][ts:ts1,:,:],axis=0)
                    time = rootgrp["time"][:][ts]
                    unit = rootgrp[variable].units

                    y_var = rootgrp.variables[variable].dimensions[2]
                    z_var = rootgrp.variables[variable].dimensions[1]

                    y_domain = rootgrp[y_var][:]
                    z_domain = rootgrp[z_var][:]

                    z_domain_unit = rootgrp[z_var].units
                    y_domain_unit = rootgrp[y_var].units

                else:
                    cross_data = \
                            np.hstack(
                                    (cross_data,np.mean(rootgrp[variable][:][ts:ts1,:,:],axis=0)))
                    y_domain = np.append(y_domain,rootgrp[y_var][:])


    else:
        pass
        #TODO handle error

    crosssection = { 
            "data" : cross_data,
            "time" : time,
            "unit" : unit}

    if direction == 'xy':
        crosssection["x_domain"] = x_domain
        crosssection["x_domain_unit"] = x_domain_unit
        crosssection["y_domain"] = y_domain
        crosssection["y_domain_unit"] = y_domain_unit

    elif direction == 'xz':
        crosssection["x_domain"] = x_domain
        crosssection["x_domain_unit"] = x_domain_unit
        crosssection["z_domain"] = z_domain
        crosssection["z_domain_unit"] = z_domain_unit

    elif direction == 'yz':
        crosssection["y_domain"] = y_domain
        crosssection["y_domain_unit"] = y_domain_unit
        crosssection["z_domain"] = z_domain
        crosssection["z_domain_unit"] = z_domain_unit
    else:
        pass

    if averaged:
        crosssection["time1"] = time1
    
    return crosssection
