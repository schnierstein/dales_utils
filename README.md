# dales_utils


## Quicklook script

use as `python3 quicklooks_DALES_v3.py --sim_dir /PATH/TO/SIMULATION --zoom_height XXX`

Creates quicklooks folder in simulation directory with plots of various variables.

Zoom height is height of "interesting" features, shows up as extra subplot.

BUG: Crashed for some plots when there is no liquid water, still usable.


# dales_data

Script that help with reading in dales data, especially combining crosssections and fielddumps. They are ugly and slow, but they work (mostly).

Best look at the example and documention in code.


# dales run JUWELS

Three scripts working like this:

`runDALES.sh` is the script that is submitted to SLURM via `sbatch`.

If job has to be split, use `submission_script.sh` and sync with the `change_nml.py` script.
Creates jobs dependent on each other and `change_nml.py` changes the namelist between jobs to account for the right restarting file. Can change output options between jobs as well.
