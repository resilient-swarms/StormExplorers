import pathlib
import sys
module_dir = pathlib.Path(__file__).parent.resolve()
root_dir = module_dir.parent
model_dir = root_dir.joinpath("models")
asvlite_wrapper_dir = root_dir.joinpath("dependency", "ASVLite", "wrapper", "cython")
sys.path.insert(0, str(asvlite_wrapper_dir))

import glob
import pandas as pd

from storms_archive import Storms_archieve
from storm_tracking import Storm_tracking
from thrust_calibrator import Thrust_calibrator
from rudder_controller import Rudder_PID_controller
from asv import py_Asv, py_Asv_specification
from geometry import py_Coordinates_3D

# Wave glider specs
asv_spec = py_Asv_specification()
asv_spec.L_wl       = 2.1 # m
asv_spec.B_wl       = 0.6 # m
asv_spec.D          = 0.25 # m 
asv_spec.T          = 0.09 # m
asv_spec.max_speed  = 4.0 # m/s
asv_spec.disp       = 0.09 # m3
asv_spec.r_roll     = 0.2 # m
asv_spec.r_pitch    = 0.6 # m
asv_spec.r_yaw      = 0.6 # m
asv_spec.cog        = py_Coordinates_3D(1.05, 0.0, -3.0) # m

# Tune the rudder PID controller
rudder_controller   = Rudder_PID_controller(asv_spec) 

# Load the results of tuning the controller
rudder_tuning_results_dir = root_dir.joinpath("results", "rudder_controller", "tuning")
dfs = pd.DataFrame(columns="repeats P I D error_avg error_std".split())
i = 0
for file in glob.glob(str(rudder_tuning_results_dir)+"/*.txt"):
    df = pd.read_csv(file, header=0, names="P I D error_avg error_std".split(), sep=" ")
    for j in range(len(df)):
        dfs.loc[len(dfs)] = [i, df.loc[j, "P"], df.loc[j, "I"], df.loc[j, "D"], df.loc[j, "error_avg"], df.loc[j, "error_std"]]
    i += 1
# Find the PID for which the lowest error value was observed.
min_index = dfs["error_avg"].idxmin()
K = (dfs.loc[min_index, "P"], dfs.loc[min_index, "I"], dfs.loc[min_index, "D"])
print("Rudder PID controller tuning complete.")
print("P, I, D = ", K)
# Write the results to file
with open(str(model_dir.joinpath("rudder_PID")), "w+") as f:
    f.write("P I D\n")
    f.write("{} {} {}".format(K[0], K[1], K[2]))

# Calibrate wave glider thrust:
# - download wave and current data from Nov-2011 to Mar-2013 if it does not exist, 
# - generate training and validation data creating models for thrust tuning factor, 
# - generating models for thrust tuning factor. 
thrust_tuning = Thrust_calibrator(host_type="PC")

# Download storms archieve if it does not exist
storms_file = root_dir.joinpath("data", "storms")
storms_archieve = Storms_archieve(str(storms_file))
df = storms_archieve.get_storms()

# Run storm tracking experiments
storm_tracking = Storm_tracking("HPC")
storm_tracking.run_all()
storm_tracking.compute_performance()