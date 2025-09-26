import f90nml
import os


time_increment = 150  # in minutes, has to be synced with namoptions runtime

output_dicts = []
#    {
#        "namlist": "namfielddump",
#        "interval": [0, 16199],
#        "dtav": 1800,
#        "lfielddump": True,
#    },
#    {
#        "namlist": "namfielddump",
#        "interval": [16200, 26999],
#        "dtav": 450,
#        "lfielddump": True,
#    },
#    {
#        "namlist": "namfielddump",
#        "interval": [27000, 1e6],
#        "dtav": 1800,
#        "lfielddump": True,
#    },
#]


def update_startfile(startfile, increment):

    h = int(startfile[5:8])
    m = int(startfile[9:11])
    print(m)

    t_new = h * 60 + m + increment

    h_new = f"{t_new // 60:03d}"
    m_new = f"{t_new % 60:02d}"

    new_startfile = startfile[0:5] + h_new + "h" + m_new + "m" + startfile[12:]

    return [new_startfile, t_new * 60]


if os.path.exists("namoptions"):

    nml = f90nml.read("./namoptions")

    startfile = nml["run"]["startfile"]

    us = update_startfile(startfile, time_increment)

    nml["run"]["startfile"] = us[0]

    new_time = us[1]

    nml["run"]["lwarmstart"] = True

    nml["run"]["runtime"] = time_increment * 60

    for od in output_dicts:

        if (new_time >= od["interval"][0]) and (new_time <= od["interval"][1]):

            for k in od.keys():
                if k not in ["namlist", "interval"]:
                    nml[od["namlist"]][k] = od[k]

    os.remove("namoptions")

    f90nml.write(nml, "./namoptions")

else:
    print("ERROR: No namoptions")
