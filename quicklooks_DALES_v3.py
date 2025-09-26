import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import sys
import os
import datetime
import matplotlib as mpl
import warnings

import matplotlib.patches as patches

warnings.filterwarnings("ignore", category=DeprecationWarning)

# plt.rcParams.update({"text.usetex" : True})
# mpl.rcParams.update({'text.latex.preamble': r"\usepackage{amsmath,siunitx,microtype}"})
mpl.rcParams.update({"font.size": 14})

import matplotlib.gridspec as gridspec

import argparse

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


def profile_time(
    variable,
    sdir,
    filename="profiles.001.nc",
    iexpnr=1,
    contour_variable="xxx",
    difference_variable="xxx",
    addition_variable="xxx",
    cmap="viridis",
    unit=" ",
    data_scale=1,
    addition_data_scale=-9999,
    title=" ",
    variable_name=" ",
    contour_variable_name=" ",
    data_scale_contour=1,
    diverging_cmap=False,
    contour_unit=" ",
    **parser_dic,
):

    if parser_dic["no_profiles"]:
        return

    data = nc.Dataset(os.path.join(sdir, filename))

    gridspec = {
        "hspace": 0.25,
        "wspace": 0.25,
        "left": 0.1,
        "right": 0.90,
        "top": 0.90,
        "bottom": 0.1,
        "width_ratios": [1, 1.4],
        "height_ratios": [1, 1],
    }

    # fig, axes = plt.subplots(2,1, figsize=(13.8,9.8))

    gs_kw = dict(width_ratios=[1.4, 1], height_ratios=[1, 1])

    fig, axd = plt.subplot_mosaic(
        [["upper left", "upper right"], ["lower left", "lower right"]],
        gridspec_kw=gridspec,
        figsize=(19.2, 10.8),
        constrained_layout=False,
    )
    axes = [axd["upper right"], axd["lower right"]]

    pax = [axd["upper left"], axd["lower left"]]

    d = data[variable][:] * data_scale

    time = data["time"][:] / 3600
    zm = data["zm"][:]
    zt = data["zt"][:]

    for axe in axes:

        axe.grid()
        axe.set_ylabel(r"Height [m]")
        axe.set_xlabel(r"Time [h]")

    Z, T = np.meshgrid(zm, time)

    d = np.rot90(d)
    Z = np.rot90(Z)
    T = np.rot90(T)

    variable_string = variable

    ind = len(zm[zm < parser_dic["zoom_height"]])
    ind = len(zm) - ind

    if variable == "qt":
        d[d < 0] = 0

    if difference_variable != "xxx":

        d = d - np.rot90(data[difference_variable][:]) * data_scale
        variable_string = variable_string + "-" + difference_variable

    if addition_variable != "xxx":

        if addition_data_scale == -9999:
            addition_data_scale = data_scale

        if addition_data_scale < 0:
            sign = "-"
        elif addition_data_scale > 0:
            sign = "+"

        d = d + np.rot90(data[addition_variable][:]) * addition_data_scale
        variable_string = variable_string + sign + addition_variable

    if contour_variable != "xxx":

        d2 = np.rot90(data[contour_variable][:]) * data_scale_contour
        variable_string = variable_string + "-contour-" + contour_variable

    max_val = np.amax(d)
    min_val = np.amin(d)

    if (abs(max_val + min_val) != (abs(max_val) + abs(min_val))) or diverging_cmap:

        vmax = np.amax([abs(max_val), abs(min_val)])
        vmin = -1 * np.amax([abs(max_val), abs(min_val)])

    else:

        vmax = max_val
        vmin = min_val

    pc0 = axes[0].pcolormesh(
        T, Z, d, shading="nearest", vmax=vmax, vmin=vmin, cmap=cmap
    )

    tmin = np.amin(time)
    tmax = np.amax(time)

    colors_tind = ["C0", "C1", "C2", "C3", "C4", "C5"]

    count_tind = 0

    xmax = -1e10
    xmin = 1e10

    for j in [0.1, 0.3, 0.5, 0.7, 0.9]:
        # for j in [0, 0.125, 0.25]:

        tind = tmin + (tmax - tmin) * j
        tind = np.argmin(np.abs(time - tind))
        label_tind = f"t={time[tind]:3.2f}"

        if np.amax(d[:, tind]) > xmax:
            xmax = np.amax(d[:, tind])

        if np.amin(d[:, tind]) < xmin:
            xmin = np.amin(d[:, tind])

        pax[0].plot(
            np.flip(d[:, tind]), zm[:], label=label_tind, color=colors_tind[count_tind]
        )

        count_tind += 1

    if xmax != xmin:
        xdist = np.abs(np.abs(xmin) - np.abs(xmax))
        pax[0].set_xlim([xmin - 0.05 * xdist, xmax + 0.05 * xdist])

    del xmax, xmin

    if parser_dic["mosaic_mode"] == True:

        if variable == "qt":
            rs_var = "qv"

            rs_dset = nc.Dataset(parser_dic["rs_file"], "r")
            rs_time = parser_dic["rs_time"]
            rs_data = rs_dset[rs_var][:] * 1000
            rs_height = rs_dset["GPSHgt"][:]

            tind = np.argmin(np.abs(time - 7))

            pax[0].plot(rs_data, rs_height, label=f"RS, {rs_time} UTC", color="black")
            pax[1].plot(
                rs_data[rs_height < parser_dic["zoom_height"]],
                rs_height[rs_height < parser_dic["zoom_height"]],
                label=f"RS, {rs_time} UTC",
                color="black",
            )

            pax[0].plot(np.flip(d[:, tind]), zm[:], label="7h sim time", color="C3")
            pax[1].plot(
                np.flip(d[ind:, tind]),
                zm[: len(zm) - ind],
                label="7h sim time",
                color="C3",
            )

        if variable == "ql":
            rs_var = "RH"

            rs_dset = nc.Dataset(parser_dic["rs_file"], "r")
            rs_time = parser_dic["rs_time"]
            rs_data = rs_dset[rs_var][:]
            rs_height = rs_dset["GPSHgt"][:]

            tind = np.argmin(np.abs(time - 7))

            cloud_points = np.array([])

            if rs_data[0] == 100:
                cloud_points = np.append(cloud_points, rs_height[0])

            cloud_threshold = parser_dic["RH_cloud_threshold"]

            for i in range(1, len(rs_data) - 1):

                if rs_data[i] >= cloud_threshold:

                    if (
                        rs_data[i - 1] >= cloud_threshold
                        and rs_data[i + 1] < cloud_threshold
                    ) or (
                        rs_data[i - 1] < cloud_threshold
                        and rs_data[i + 1] >= cloud_threshold
                    ):
                        cloud_points = np.append(cloud_points, rs_height[i])

            if rs_data[-1] == 100:
                cloud_points = np.append(cloud_points, rs_height[-1])

            cloud_points = cloud_points.reshape((len(cloud_points) // 2, 2))

            first_cloud = True

            for arr in cloud_points:

                cloud_poly = [
                    [0 - 1e3, arr[0]],
                    [0 - 1e3, arr[1]],
                    [1e3, arr[1]],
                    [1e3, arr[0]],
                ]

                if first_cloud:
                    rh_label = f"RS RH>{cloud_threshold:.1f}%"
                else:
                    rh_label = None

                poly0 = patches.Polygon(
                    cloud_poly,
                    alpha=0.3,
                    closed=True,
                    label=rh_label,
                )
                pax[0].add_patch(poly0)

                poly1 = patches.Polygon(
                    cloud_poly,
                    alpha=0.3,
                    closed=True,
                    label=rh_label,
                )
                pax[1].add_patch(poly1)

                first_cloud = False

        if variable == "thv":
            rs_var = "T"

            rs_dset = nc.Dataset(parser_dic["rs_file"], "r")
            rs_time = parser_dic["rs_time"]
            rs_data = rs_dset[rs_var][:]
            rs_height = rs_dset["GPSHgt"][:]
            rs_pressure = rs_dset["p"][:]
            rs_qv = rs_dset["qv"][:]

            rs_data = (
                rs_data * (rs_pressure[0] / rs_pressure) ** 0.2854 * (1 + 0.61 * rs_qv)
            )

            pax[0].plot(rs_data, rs_height, label=f"RS, {rs_time} UTC", color="black")
            pax[1].plot(
                rs_data[rs_height < parser_dic["zoom_height"]],
                rs_height[rs_height < parser_dic["zoom_height"]],
                label=f"RS, {rs_time} UTC",
                color="black",
            )

            tind = np.argmin(np.abs(time - 7))

            pax[0].plot(np.flip(d[:, tind]), zm[:], label="7h sim time", color="C3")
            pax[1].plot(
                np.flip(d[ind:, tind]),
                zm[: len(zm) - ind],
                label="7h sim time",
                color="C3",
            )

    pax[0].set_ylim([zm[0], zm[-1]])
    pax[0].grid()

    pax[0].legend()
    pax[0].set_ylabel("Height [m]")
    pax[0].set_xlabel(variable_name + " " + unit)

    divider0 = make_axes_locatable(axes[0])
    cax0 = divider0.append_axes("right", size="5%", pad=0.05)
    cbar0 = plt.colorbar(pc0, cax=cax0)

    cbar0.ax.get_yaxis().labelpad = 15
    cbar0.ax.set_ylabel(variable_name + " " + unit, rotation=270)

    if contour_variable != "xxx":

        minprof = np.amin(d2)
        maxprof = np.amax(d2)

        # if minprof == maxprof:
        #    minprof = minprof - 1
        #    maxprof = maxprof + 1

        levels = [
            minprof + abs(maxprof - minprof) * 0.05,
            # minprof+abs(maxprof-minprof)*0.25,
            (minprof + maxprof) / 2,
            # minprof+abs(maxprof-minprof)*0.75,
            minprof + abs(maxprof - minprof) * 0.95,
        ]

        cp0 = axes[0].contour(T, Z, d2, levels=levels, cmap="gist_rainbow")

        cax01 = divider0.append_axes("right", size="5%", pad=0.90)
        cbar01 = plt.colorbar(cp0, cax=cax01)

        cbar01.ax.get_yaxis().labelpad = 15
        cbar01.ax.set_ylabel(contour_variable_name + " " + contour_unit, rotation=270)

    max_val = np.amax(d[ind:, :])
    min_val = np.amin(d[ind:, :])

    # print(f"{variable}: max value {max_val} , min value {min_val}")

    if (abs(max_val + min_val) != (abs(max_val) + abs(min_val))) or diverging_cmap:

        vmax = np.amax([abs(max_val), abs(min_val)])
        vmin = -1 * np.amax([abs(max_val), abs(min_val)])

    else:

        vmax = max_val
        vmin = min_val

    pc1 = axes[1].pcolormesh(
        T[ind:, :],
        Z[ind:, :],
        d[ind:, :],
        shading="nearest",
        vmax=vmax,
        vmin=vmin,
        cmap=cmap,
    )

    tmin = np.amin(time)
    tmax = np.amax(time)

    colors_tind = ["C0", "C1", "C2", "C3", "C4", "C5"]

    count_tind = 0

    xmax = -1e10
    xmin = 1e10

    for j in [0.1, 0.3, 0.5, 0.7, 0.9]:
        # for j in [0, 0.125, 0.25]:

        tind = tmin + (tmax - tmin) * j
        tind = np.argmin(np.abs(time - tind))
        label_tind = f"t={time[tind]:3.2f}"

        if np.amax(d[ind:, tind]) > xmax:
            xmax = np.amax(d[ind:, tind])

        if np.amin(d[ind:, tind]) < xmin:
            xmin = np.amin(d[ind:, tind])

        pax[1].plot(
            np.flip(d[ind:, tind]),
            zm[: len(zm) - ind],
            label=label_tind,
            color=colors_tind[count_tind],
        )

        count_tind += 1

    if xmax != xmin:
        xdist = np.abs(np.abs(xmin) - np.abs(xmax))
        pax[1].set_xlim([xmin - 0.05 * xdist, xmax + 0.05 * xdist])

    del xmax, xmin

    pax[1].set_ylim([zm[0], zm[len(zm) - ind - 1]])

    pax[1].legend()
    pax[1].grid()

    pax[1].set_ylabel("Height [m]")
    pax[1].set_xlabel(variable_name + " " + unit)

    axes[0].set_ylim([np.amin(Z), np.amax(Z)])
    axes[1].set_ylim([np.amin(T), np.amax(Z[ind:, :])])

    axes[0].set_xlim([np.amin(T), np.amax(T)])
    axes[1].set_xlim([np.amin(T[ind:, :]), np.amax(T[ind:, :])])

    divider1 = make_axes_locatable(axes[1])
    cax1 = divider1.append_axes("right", size="5%", pad=0.05)
    cbar1 = plt.colorbar(pc1, cax=cax1)

    cbar1.ax.get_yaxis().labelpad = 15

    cbar1.ax.set_ylabel(variable_name + " " + unit, rotation=270)

    if contour_variable != "xxx":

        minprof = np.amin(d2[ind:, :])
        maxprof = np.amax(d2[ind:, :])

        levels = [
            minprof + abs(maxprof - minprof) * 0.05,
            # minprof+abs(maxprof-minprof)*0.25,
            (minprof + maxprof) / 2,
            # minprof+abs(maxprof-minprof)*0.75,
            minprof + abs(maxprof - minprof) * 0.95,
        ]

        cp1 = axes[1].contour(
            T[ind:, :], Z[ind:, :], d2[ind:, :], levels=levels, cmap="gist_rainbow"
        )

        cax11 = divider1.append_axes("right", size="5%", pad=0.90)
        cbar11 = plt.colorbar(cp1, cax=cax11)

        cbar11.ax.get_yaxis().labelpad = 15
        cbar11.ax.set_ylabel(contour_variable_name + " " + contour_unit, rotation=270)

    fig.suptitle(sdir.split("/")[-2] + "\n" + title)

    # plt.colorbar(pc0)
    try:
        os.mkdir(sim_dir + "quicklooks/profile_time/")
    except:
        print("Could not create profile_time sub-folder; might already exist")

    plt.savefig(
        sim_dir + f"quicklooks/profile_time/profile_time_single_{variable_string}.png",
        dpi=300,
    )
    plt.close(fig)


def timeseries(
    variables,
    sdir,
    tmfile,
    unit=" ",
    title=" ",
    labels=[" "],
    data_scales=[-1],
    ylabel=" ",
    **parser_dic,
):

    if parser_dic["no_timeseries"]:
        return

    fig, ax = plt.subplots(1, 1, figsize=(7.2, 7.2))

    if "nc" in tmfile:

        data = nc.Dataset(sdir + tmfile)

    elif "tmsurf" in tmfile:

        data = np.genfromtxt(sdir + tmfile, names=True)

    else:

        print(f"File {tmfile} not found")
        return

    time = data["time"][:] / 3600

    max_value = -1e30
    min_value = 1e30

    count = 0

    variable_string = ""

    for var in variables:

        variable_string += var

        if labels[0] == " ":
            l = ""
            for i in range(len(var.split("_"))):
                l += var.split("_")[i]
            label = l

        else:
            label = labels[count]

        d = data[var][:]

        if data_scales[0] != -1:
            d = d * data_scales[count]

        ax.plot(time, d, label=label)

        if np.amax(d) > max_value:
            max_value = np.amax(d)

        if np.amin(d) < min_value:
            min_value = np.amin(d)

        count += 1

    ax.set_xlim([time[0], time[-1]])
    y_max = max_value + abs(max_value) * 0.05
    y_min = min_value - abs(max_value) * 0.05

    ax.set_ylim([y_min, y_max])

    ax.set_xlabel("Time [h]")
    ax.set_ylabel(ylabel)

    # ax.plot([11,11],[0,15000], color="black", linestyle="dotted")

    ax.grid()
    ax.legend()

    fig.suptitle(title)

    try:
        os.mkdir(sim_dir + "quicklooks/tmser/")
    except:
        print("Could not create tmser sub-folder; might already exist")

    plt.savefig(sim_dir + f"quicklooks/tmser/tmser_{variable_string}.png", dpi=300)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Modify plotting script")

    required_args_group = parser.add_argument_group("Required arguments")
    mosaic_mode_group = parser.add_argument_group("MOSAiC mode arguments")

    required_args_group.add_argument(
        "-sd",
        "--sim_dir",
        type=str,
        required=True,
        help="File directory which contains simulation output",
    )

    parser.add_argument(
        "-nop",
        "--number_of_processes",
        type=int,
        required=False,
        help="Number of MPI processes used for the simulation",
    )

    parser.add_argument(
        "--mosaic_mode",
        action="store_true",
        required=False,
        default=False,
        help="Enables addional plotting capabilities with MOSAiC data, requires path to observational data",
    )

    parser.add_argument(
        "-zh",
        "--zoom_height",
        type=float,
        required=False,
        default=1500,
        help="Defines the maximum height shown for the zoomed-in plots",
    )

    parser.add_argument(
        "--no_profiles",
        required=False,
        default=False,
        action="store_true",
        help="Choose if profiles should not be plotted; False by default",
    )

    parser.add_argument(
        "--no_timeseries",
        required=False,
        default=False,
        action="store_true",
        help="Choose if time series should not be plotted; False by default",
    )

    mosaic_mode_group.add_argument(
        "--rs_file",
        type=str,
        required=False,
        help="MOSAiC radiosonde file, if MOSAiC mode is enabled",
    )

    mosaic_mode_group.add_argument(
        "--rs_time",
        type=int,
        required=False,
        help="In MOSAiC mode, defines time of radiosonde",
    )

    mosaic_mode_group.add_argument(
        "--RH_cloud_threshold",
        type=float,
        required=False,
        default=98.0,
        help="In MOSAiC mode, sets relative humidity threshold"
        + " recognized as cloud (default=98)",
    )

    P = parser.parse_args()

    if P.mosaic_mode and ((P.rs_file is None) or (P.rs_time is None)):
        parser.error("--mosaic_mode requires --rs_file and --rs_time to be set")

    sim_dir = P.sim_dir

    if sim_dir[-1] != "/":
        sim_dir = sim_dir + "/"

    try:
        os.mkdir(sim_dir + "quicklooks")
    except:
        print("Could not create quicklooks sub-folder; might already exist")

    variables = [
        "ql",
        "qt",
        "sv008",
        "sv010",
        "sv012",
        "lwu",
        "lwd",
        "swu",
        "swd",
        "thltend",
        "thl",
        "sv006",
        "u",
        "v",
        "w2r",
        "u2r",
        "v2r",
        "thv",
        "wqtt",
        "wthlt",
    ]

    titles = [
        "Liquid water",
        "Total water",
        "Ice water",
        "Snow water",
        "Graupel water",
        "LW up",
        "LW down",
        "SW up",
        "SW down",
        "Total radiative tendency",
        "Liquid water potential temperature",
        "CCN Number concentration",
        "West-East velocity",
        "South-North velocity",
        "w2r",
        "u2r",
        "v2r",
        "thv",
        "wqtt",
        "wthlt",
    ]

    units = [
        r"g/kg",
        r"g/kg",
        r"g/kg",
        r"g/kg",
        r"g/kg",
        "W/m2",
        "W/m2",
        "W/m2",
        "W/m2",
        "K/h",
        "K",
        "Unit",
        "m/s",
        "m/s",
        "m2/s2",
        "m2/s2",
        "m2/s2",
        "K",
        "g/kg m/s",
        "K m/s",
    ]

    data_scales = [
        1e3,
        1e3,
        1e3,
        1e3,
        1e3,
        1,
        1,
        1,
        1,
        3600,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1000,
        1,
    ]

    variable_names = [
        "Liquid water",
        "Total water",
        "Ice Water",
        "Snow Water",
        "Graupel Water",
        "LW up",
        "LW down",
        "SW up",
        "SW down",
        "Tot. rad. tend.",
        "Liquid water pot. temp.",
        "CCN Number concentration",
        "West-East velocity",
        "South-North velocity",
        "w2r",
        "u2r",
        "v2r",
        "Virtual potential temperature",
        "Total moisture flux",
        "Total Theta L flux",
    ]

    for i in range(len(variables)):

        profile_time(
            variables[i],
            sim_dir,
            data_scale=data_scales[i],
            unit=units[i],
            title=titles[i],
            variable_name=variable_names[i],
            **vars(P),
        )

    profile_time(
        "qt",
        sim_dir,
        addition_variable="ql",
        data_scale=1000,
        addition_data_scale=-1000,
        unit="g/kg",
        variable_name="Water vapor",
        title="Water vapor",
        **vars(P),
    )

    profile_time(
        "lwu",
        sim_dir,
        unit=r"W/m2",
        addition_variable="lwd",
        variable_name="LW net",
        title="Longwave net radiation",
        diverging_cmap=True,
        cmap="seismic",
        **vars(P),
    )

    profile_time(
        "swu",
        sim_dir,
        unit=r"W/m2",
        addition_variable="swd",
        variable_name="SW net",
        diverging_cmap=True,
        title="Shortwave net radiation",
        cmap="seismic",
        **vars(P),
    )

    profile_time(
        "qt",
        sim_dir,
        contour_variable="ql",
        data_scale=1000,
        data_scale_contour=1000,
        unit="g/kg",
        contour_unit="g/kg",
        variable_name="Total water",
        contour_variable_name="Liquid water",
        title="Total water with liquid water contour",
        **vars(P),
    )

    profile_time(
        "thltend",
        sim_dir,
        contour_variable="ql",
        data_scale=3600,
        data_scale_contour=1000,
        diverging_cmap=True,
        cmap="seismic",
        unit="K/h",
        contour_unit="g/kg",
        variable_name="Tot. rad. tend.",
        contour_variable_name="Liquid water",
        title="Total radiative tendency with liquid water contour",
        **vars(P),
    )

    profile_time(
        "qt",
        sim_dir,
        contour_variable="sv008",
        data_scale=1000,
        data_scale_contour=1000,
        unit="g/kg",
        contour_unit="g/kg",
        variable_name="Total water",
        contour_variable_name="Ice water",
        title="Total water with ice water contour",
        **vars(P),
    )

    profile_time(
        "ql",
        sim_dir,
        contour_variable="sv008",
        data_scale=1000,
        data_scale_contour=1000,
        unit="g/kg",
        contour_unit="g/kg",
        variable_name="Liquid water",
        contour_variable_name="Ice water",
        title="Liquid water with ice water contour",
        **vars(P),
    )

    profile_time(
        "lwu",
        sim_dir,
        addition_variable="lwd",
        contour_variable="ql",
        data_scale=1,
        data_scale_contour=1000,
        cmap="seismic",
        diverging_cmap=True,
        unit="W/m2",
        contour_unit="g/kg",
        variable_name="LW net",
        contour_variable_name="Liquid water",
        title="Net longwave radiation with liquid water contour",
        **vars(P),
    )

    profile_time(
        "lwu",
        sim_dir,
        addition_variable="lwd",
        contour_variable="sv008",
        data_scale=1,
        data_scale_contour=1000,
        cmap="seismic",
        diverging_cmap=True,
        unit="W/m2",
        contour_unit="g/kg",
        variable_name="LW net",
        contour_variable_name="Ice water",
        title="Net longwave radiation with ice water contour",
        **vars(P),
    )

    profile_time(
        "swu",
        sim_dir,
        addition_variable="swd",
        contour_variable="ql",
        data_scale=1,
        data_scale_contour=1000,
        cmap="seismic",
        diverging_cmap=True,
        unit="W/m2",
        contour_unit="g/kg",
        variable_name="SW net",
        contour_variable_name="Liquid water",
        title="Net shortwave radiation with liquid water contour",
        **vars(P),
    )

    timeseries(
        ["zi", "zb", "zc_av"],
        sim_dir,
        tmfile="tmser.001.nc",
        ylabel="Height",
        **vars(P),
    )

    timeseries(
        ["ctot_frac", "ci_frac", "cl_frac"],
        sim_dir,
        tmfile="tmser.001.nc",
        ylabel="Cloud fraction",
        **vars(P),
    )

    timeseries(
        ["wthvs", "wqls"],
        sim_dir,
        tmfile="tmsurf.001",
        ylabel="Surface flux in W/m2",
        data_scales=[1.216e3, 2.45e6 * 1.15],
        **vars(P),
    )

    timeseries(
        ["lwp_bar"],
        sim_dir,
        tmfile="tmser.001.nc",
        ylabel="Liquid path g/m2",
        data_scales=[1000],
        **vars(P),
    )
    timeseries(
        ["icwp_bar"],
        sim_dir,
        tmfile="tmser.001.nc",
        ylabel="Ice path g/m2",
        data_scales=[1000],
        **vars(P),
    )
    timeseries(
        ["giwp_bar"],
        sim_dir,
        tmfile="tmser.001.nc",
        ylabel="Graupel path g/m2",
        data_scales=[1000],
        **vars(P),
    )
    timeseries(
        ["siwp_bar"],
        sim_dir,
        tmfile="tmser.001.nc",
        ylabel="Snow path g/m2",
        data_scales=[1000],
        **vars(P),
    )
