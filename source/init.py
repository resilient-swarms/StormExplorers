import pathlib
import sys
module_dir = pathlib.Path(__file__).parent.resolve()
root_dir = module_dir.parent
model_dir = root_dir.joinpath("model")
asvlite_wrapper_dir = root_dir.joinpath("dependency", "ASVLite", "wrapper", "cython")
sys.path.insert(0, str(asvlite_wrapper_dir))

from storms_archive import Storms_archieve
from thrust_calibrator import Thrust_calibrator
from rudder_controller import Rudder_PID_controller
from asv import py_Asv, py_Asv_specification
from geometry import py_Coordinates_3D

# Tune the rudder PID controller
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
rudder_controller   = Rudder_PID_controller(asv_spec) # Initialise and tune the controller. 

# Runing this will:
# - download wave and current data from Nov-2011 to Mar-2013 if it does not exist, 
# - generate training and validation data creating models for thrust tuning factor, 
# - generating models for thrust tuning factor. 
thrust_tuning = Thrust_calibrator(host_type="PC")

# Download storms archieve if it does not exist
storms_file = root_dir.joinpath("data", "storms")
storms_archieve = Storms_archieve(str(storms_file))
df = storms_archieve.get_storms()
