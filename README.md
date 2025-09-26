# dales_utils


## Quicklook script

use as `python3 quicklooks_DALES_v3.py --sim_dir /PATH/TO/SIMULATION --zoom_height XXX`

Creates quicklooks folder in simulation directory with plots of various variables.

Zoom height is height of "interesting" features, shows up as extra subplot.

BUG: Crashed for some plots when there is no liquid water, still usable.


# dales_data

Script that help with reading in dales data, especially combining crosssections and fielddumps. They are ugly and slow, but they work (mostly).

Best look at the example and documention in code.

