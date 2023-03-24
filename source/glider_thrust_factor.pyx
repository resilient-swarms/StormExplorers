import pathlib
import sys
module_dir = pathlib.Path(__file__).parent.resolve()
root_dir = module_dir.parent
model_dir = root_dir.joinpath("model")
asvlite_wrapper_dir = root_dir.joinpath("dependency", "ASVLite", "wrapper", "cython")
sys.path.insert(0, str(asvlite_wrapper_dir))
import os
import math
cimport epsg
import pyproj
import random
from concurrent.futures import as_completed, ProcessPoolExecutor
from tqdm import tqdm # to show a progress bar
from cpython.datetime cimport timedelta
import pandas as pd
import numpy as np
from rudder_controller import Rudder_PID_controller
from netcdf cimport NetCDF_wave, NetCDF_current
cimport cds
from sea_surface cimport py_Sea_surface
from asv cimport py_Asv_specification, py_Asv
from geometry cimport py_Coordinates_3D 
from geometry import py_normalise_angle_PI, py_normalise_angle_2PI
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from joblib import load

_MPS_TO_KNOTS = 1.94384449 # 1 m/s = 1.94384449 knots

# Load the model for thrust tuning
_file_name = model_dir.joinpath("thrust_tuning_lin_reg_4.joblib")
_thrust_tuning_factor_model = load(str(_file_name)) 

def get_thrust_tuning_factor(wave_ht, zonal_current_velocity, meridional_current_velocity, vehicle_heading):
    current_speed = math.sqrt(zonal_current_velocity*zonal_current_velocity + meridional_current_velocity*meridional_current_velocity) * _MPS_TO_KNOTS
    theta = math.atan2(zonal_current_velocity, meridional_current_velocity)
    theta = py_normalise_angle_2PI(theta)
    relative_current_direction = abs(py_normalise_angle_PI(theta - vehicle_heading)) # (0, PI) radians
    return _thrust_tuning_factor_model.predict([[wave_ht, current_speed, relative_current_direction]])


# ===================================================================================

# Compute the thrust tuning factors for each of the sea states by comparing simulated 
# speed with the real-world speed.

# ===================================================================================

# These variables are declared globaly and not as a member of class __Thrust_tuning
# to reduce memory consumption. While multiprocessing using concurrent.futures.ProcessPoolExecutor,
# if this variable was in the class then the large numpy arrays in NetCDF will be copied several times
# thereby increasing the memory consumption. 
# Care should be taken to make sure this vairalbe is only written to while NOT multi-processing. 
_nc_data_wave = {} 
_nc_data_current = {}
class _Thrust_tuning:
    def __init__(self, host_type):
        '''
        Valid values for host_type are the strings - "PC", "HPC"
        '''
        self.host_type = host_type
        self.geodesic = pyproj.Geod(ellps='WGS84') # Initialise it once here. Used in __get_heading().
        self.proximity_margin = 20.0 # m
        # Files and directories
        self.__download_wave_data()
        self.__download_ocean_current_data()
        self.output_dir = root_dir.joinpath(*"results/glider_thrust/".split("/"))    
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        # Dataframe to hold both simulation and real-world data
        self.df = None
    
    def __load_log_files(self, log_files):
        for log_file in log_files:
            # Load log data to dataframe
            df = pd.read_csv(log_file)
            df["Timestamp(UTC)"] = pd.to_datetime(df["Timestamp(UTC)"]) # Convert column Timestamp to datetime format
            df["year"] = pd.DatetimeIndex(df["Timestamp(UTC)"]).year
            df["month"] = pd.DatetimeIndex(df["Timestamp(UTC)"]).month
            df["Latitude2"] = df[" Latitude"].shift(-1)
            df["Longitude2"] = df[" Longitude"].shift(-1)
            coordinates_PCS = df.apply(lambda row: epsg.GCS_to_PCS(row[" Longitude"], row[" Latitude"]), axis=1)
            df["p1.x"] = [coordinate[0] for coordinate in coordinates_PCS]
            df["p1.y"] = [coordinate[1] for coordinate in coordinates_PCS]
            # df["p1.x"], df["p1.y"] = df.apply(lambda row: epsg.GCS_to_PCS(row[" Longitude"], row[" Latitude"]), axis=1)
            df["p2.x"] = df["p1.x"].shift(-1)
            df["p2.y"] = df["p1.y"].shift(-1)
            df["delta_x"] = df.apply(lambda row: row["p2.x"] - row["p1.x"] if (row["p2.x"] * row["p1.x"] >= 0.0) else row["p2.x"] + row["p1.x"], axis=1)
            df["delta_y"] = df.apply(lambda row: row["p2.y"] - row["p1.y"] if (row["p2.y"] * row["p1.y"] >= 0.0) else row["p2.y"] + row["p1.y"], axis=1)
            df["distance(m)"] = ((df["delta_x"]**2 + df["delta_y"]**2)**0.5)
            df["delta_T(s)"] = (df["Timestamp(UTC)"].shift(-1) - df["Timestamp(UTC)"]).astype('timedelta64[s]')
            df["delta_T_simulated(s)"] = pd.NA # Will be set in rull_all()
            df["speed(knots)"] = _MPS_TO_KNOTS * df["distance(m)"]/ df["delta_T(s)"]
            df["speed_simulated(knots)"] = pd.NA # Will be set in rull_all()
            # Compute the initial heading for the wave glider
            df["heading(northing)"] = df.apply(lambda row: self.__get_heading(row[" Longitude"], row[" Latitude"], row["Longitude2"], row["Latitude2"]), axis=1)
            # Read the significant wave heights from the NetCDF files.
            for group_name, group in df.groupby(["year", "month"]):
                year, month = group_name
                file_name = "{}_{}.nc".format(year, str(month).zfill(2))
                self.__load_wave_nc_files([file_name])
                rows = group.index
                for i in rows:
                    df.loc[i, "simulated_wave_height(m)"] = self.__get_wave_data_at(df.loc[i, " Longitude"], df.loc[i, " Latitude"], df.loc[i, "Timestamp(UTC)"])[0]
            # Read the zonal and meridional velocities from the NetCDF files.
            for group_name, group in df.groupby(["year", "month"]):
                year, month = group_name
                self.__load_current_nc_files([(year, month)])
                rows = group.index
                for i in rows:
                    v_zonal, v_meridional = self.__get_ocean_current_at(df.loc[i, " Longitude"], df.loc[i, " Latitude"], df.loc[i, "Timestamp(UTC)"])
                    v = math.sqrt(v_zonal*v_zonal + v_meridional*v_meridional) # m/s
                    theta = math.atan2(v_zonal, v_meridional)
                    theta = py_normalise_angle_2PI(theta)
                    df.loc[i, "current(knots)"] = v * _MPS_TO_KNOTS # knots
                    df.loc[i, "current_direction(northing)"] = theta # radians
            # Relative current direction
            for i in range(len(df)):
                # range is (-PI, PI)
                df.loc[i, "relative_current_direction"] = py_normalise_angle_PI(df.loc[i, "current_direction(northing)"] - df.loc[i, "heading(northing)"])
            # Set the appropriate tuning factor for each row
            df["tuning_factor"] = 1.0 # Initialise 
            df["error_msg"] = pd.NA 
            # Clean the data
            df.drop(df.tail(1).index,inplace=True) # Drop the last row as it has incomplete data
            # df = df[df.month != 8] # Delete data for month of August as the data is incomplete.
            df = df[(df["speed(knots)"] > 0.0) & (df["speed(knots)"] < 10)] # Delete rows where 0m/s <= speed >= 10m/s. 
                                            # Wave glide will never achieve this speed and 
                                            # therefore this could be error in the data. 
            df = df[(df["delta_T(s)"] >= 1800) & (df["delta_T(s)"] <= 7200)] # Delete rows where time taken to waypoint is greater than 2hrs or less than 30 mins.
            df = df[df["distance(m)"] > self.proximity_margin] # Delete all rows where the distance to the waypoint  
                                                               # is less than the proximity margin. 
            df.reset_index(drop=True, inplace=True) # Index must be reset if rows dropped because run_all() uses iloc to iterate over rows.      
            if self.df is None:
                self.df = df
            else:
                self.df = pd.concat([self.df, df], axis=0, ignore_index=True)
            print("Loaded log file {}".format(log_file))
           
    def __download_wave_data(self):
        '''
        Downloads the wave data for the period Nov-2011 to Feb-2013 from CDS.
        '''
        nc_files = ["2011_11.nc",
                    "2011_12.nc",
                    "2012_01.nc",
                    "2012_02.nc",
                    "2012_03.nc",
                    "2012_04.nc",
                    "2012_05.nc",
                    "2012_06.nc",
                    "2012_07.nc",
                    "2012_08.nc",
                    "2012_09.nc",
                    "2012_10.nc",
                    "2012_11.nc",
                    "2012_12.nc",
                    "2013_01.nc",
                    "2013_02.nc",
                    "2013_03.nc"]
        for file in nc_files:
            cds_dir = root_dir.joinpath("data", "cds", "pacx", "waves")
            nc_file = cds_dir.joinpath(file)
            if not os.path.exists(str(nc_file)):
                # Make the cds dir if it does not exist.
                cds_dir.mkdir(parents=True, exist_ok=True)
                # nc file does not exist, download it.
                year, month = file[:-3].split("_")
                north = 40
                south = -40
                east = 180
                west = -180
                cds.get_wave_data(int(year), int(month), north, south, east, west, str(nc_file))
            else:
                print("Found {}".format(file))

    def __download_ocean_current_data(self):
        '''
        Downloads the ocean current data for the period Nov-2011 to Feb-2013 from CDS.
        '''
        tar_files = ["2011_11.tar.gz",
                    "2011_12.tar.gz",
                    "2012_01.tar.gz",
                    "2012_02.tar.gz",
                    "2012_03.tar.gz",
                    "2012_04.tar.gz",
                    "2012_05.tar.gz",
                    "2012_06.tar.gz",
                    "2012_07.tar.gz",
                    "2012_08.tar.gz",
                    "2012_09.tar.gz",
                    "2012_10.tar.gz",
                    "2012_11.tar.gz",
                    "2012_12.tar.gz",
                    "2013_01.tar.gz",
                    "2013_02.tar.gz",
                    "2013_03.tar.gz"]
        for file in tar_files:
            cds_dir = root_dir.joinpath("data", "cds", "pacx", "ocean_currents")
            tar_file = cds_dir.joinpath(file)
            if not os.path.exists(str(tar_file)):
                # Make the cds dir if it does not exist.
                cds_dir.mkdir(parents=True, exist_ok=True)
                # tar file does not exist, download it.
                year, month = file[:-7].split("_")
                cds.get_ocean_current_data(int(year), int(month), str(tar_file))
            else:
                print("Found {}".format(file))
    
    def __load_wave_nc_files(self, file_names):
        _nc_data_wave.clear()
        cds_dir = root_dir.joinpath("data", "cds", "pacx", "waves")
        for file_name in file_names:
            nc_file = cds_dir.joinpath(file_name)
            _nc_data_wave[file_name] = NetCDF_wave(str(nc_file))
            print("Loaded wave data file {}".format(file_name))
    
    def __load_current_nc_files(self, year_month_pairs):
        _nc_data_current.clear()
        cds_dir = root_dir.joinpath("data", "cds", "pacx", "ocean_currents")
        for year, month in year_month_pairs:
            days_loaded = []
            days_not_loaded = []
            for day in range(1, 32):
                file_name = "{}_{}_{}.nc".format(year, str(month).zfill(2), str(day).zfill(2))
                nc_file = cds_dir.joinpath(file_name)
                if os.path.exists(nc_file):
                    _nc_data_current[file_name] = NetCDF_current(str(nc_file))
                    days_loaded.append(day)
                else:
                    days_not_loaded.append(day)
            print("Loaded ocean current data files for the period {}-{}-[{} to {}]{}".format(year, str(month).zfill(2), days_loaded[0], days_loaded[-1], 
                  "; could not find files for the days {}".format(days_not_loaded) if len(days_not_loaded) != 0 else ""))

    def __get_wave_data_at(self, longitude, latitude, time):
        year = time.year
        month = time.month
        file_name = "{}_{}.nc".format(year, str(month).zfill(2))
        wave_hs, wave_dp = _nc_data_wave[file_name].get_wave_data_at(longitude, latitude, time)
        return (wave_hs, wave_dp)
    
    def __get_ocean_current_at(self, longitude, latitude, time):
        year = time.year
        month = time.month
        day = time.day
        file_name = "{}_{}_{}.nc".format(year, str(month).zfill(2), str(day).zfill(2))
        v_zonal, v_meridional = _nc_data_current[file_name].get_ocean_current_at(longitude, latitude, time)
        if v_zonal == None or v_meridional == None:
            print("Ocean current data not available for ({}, {}) at time {} in file {}".format(longitude, latitude, time, file_name))
            v_zonal, v_meridional = (0, 0)
        return (v_zonal, v_meridional)
    
    def __get_heading(self, long1, lat1, long2, lat2):
        fwd_azimuth, back_azimuth, distance = self.geodesic.inv(long1, lat1, long2, lat2)
        fwd_azimuth = fwd_azimuth if fwd_azimuth >= 0.0 else (360 + fwd_azimuth)
        fwd_azimuth = fwd_azimuth * math.pi/180
        return fwd_azimuth # radians

    def run(self, params):
        i, proximity_margin, update_tuning_for_sea_state, thrust_tuning_factor = params
        is_simulation_complete = False
        # Wave glider specs
        cdef py_Asv_specification asv_spec = py_Asv_specification(
                                                L_wl = 2.1,
                                                B_wl = 0.6,
                                                D = 0.25, 
                                                T = 0.09,
                                                max_speed = 4.0,
                                                disp = 0.09,
                                                r_roll = 0.2,
                                                r_pitch = 0.6,
                                                r_yaw = 0.6,
                                                cog = py_Coordinates_3D(1.05, 0.0, -3.0))
        # Start position and attitude
        cdef py_Coordinates_3D start_position_GCS = py_Coordinates_3D(self.df.loc[i, " Longitude"], self.df.loc[i, " Latitude"])
        cdef double x, y
        x, y = epsg.GCS_to_PCS(start_position_GCS.x, start_position_GCS.y)
        cdef py_Coordinates_3D start_position_PCS = py_Coordinates_3D(x, y)
        cdef py_Coordinates_3D start_attitude = py_Coordinates_3D(0.0, 0.0, math.pi/180.0 * self.df.loc[i, "heading(northing)"]) 
        # Time
        cdef double time_step_size = 40.0 # milliseconds
        simulation_start_time = self.df.loc[i, "Timestamp(UTC)"]
        # Sea surface
        cdef int wave_rand_seed = random.randint(1,100)
        cdef int count_wave_spectral_directions = 5
        cdef int count_wave_spectral_frequencies = 7
        cdef double wave_hs, wave_dp
        wave_hs, wave_dp        = self.__get_wave_data_at(start_position_GCS.x, start_position_GCS.y, simulation_start_time) 
        v_zonal, v_meridional   = self.__get_ocean_current_at(start_position_GCS.x, start_position_GCS.y, simulation_start_time)
        if wave_hs == None or wave_hs == 0.0 or wave_dp == None:
            return (is_simulation_complete, i, float("NaN"), "Wave data not available.")
        cdef py_Sea_surface sea_surface = py_Sea_surface(wave_hs, wave_dp, wave_rand_seed, count_wave_spectral_directions, count_wave_spectral_frequencies)
        # Initialise the wave glider
        cdef py_Asv asv = py_Asv(asv_spec, sea_surface, start_position_PCS, start_attitude)
        if update_tuning_for_sea_state == True:
            thrust_tuning_factor = get_thrust_tuning_factor(wave_hs, v_zonal, v_meridional, asv.py_get_attitude().z)
        asv.py_wg_set_thrust_tuning_factor(thrust_tuning_factor)
        # Initialise the rudder controller
        rudder_controller = Rudder_PID_controller(asv_spec, [1.25, 0.25, 1.75])
        cdef double rudder_angle = 0.0  
        # Waypoint
        cdef py_Coordinates_3D waypoint = py_Coordinates_3D(self.df.loc[i, "Longitude2"], self.df.loc[i, "Latitude2"])
        cdef double waypoint_x, waypoint_y
        waypoint_x, waypoint_y = epsg.GCS_to_PCS(waypoint.x, waypoint.y) 
        # Initialise time to start of simulation
        time = simulation_start_time
        max_run_time = 3 * self.df.loc[i, "delta_T(s)"]
        kill_time = time + timedelta(seconds=max_run_time)       
        # Simulate till last waypoint
        cdef py_Coordinates_3D position
        cdef double dx, dy, distance, wg_longitude, wg_latitude
        cdef double new_hs, new_dp, current_hs, current_dp
        # cdef double delta_time_sea_state_1, delta_time_sea_state_2
        cdef bint is_sea_state_same
        while time < kill_time:
            time = time + timedelta(seconds=time_step_size/1000.0)
            # Check if reached the waypoint  
            position = asv.py_get_position_cog() # Get the current position in PCS.
            dx = position.x - waypoint_x 
            dy = position.y - waypoint_y 
            distance = math.sqrt(dx*dx + dy*dy) # m
            if distance <= proximity_margin:
                # Reached waypoint
                break 
            else:
                # Get the sea state
                wg_longitude, wg_latitude = epsg.PCS_to_GCS(position.x, position.y)
                new_hs, new_dp = self.__get_wave_data_at(wg_longitude, wg_latitude, time) 
                if new_hs == None or new_hs == 0.0 or new_dp == None:
                    return (is_simulation_complete, i, float("NaN"), "Wave data not available.")
                v_zonal, v_meridional = self.__get_ocean_current_at(wg_longitude, wg_latitude, time)
                if v_zonal == None or v_zonal == float("NaN"):
                    v_zonal = 0.0
                if v_meridional == None or v_meridional == float("NaN"):
                    v_meridional = 0.0
                # Compare the sea state with the current sea state
                current_hs = sea_surface.py_get_significant_height() # m
                current_dp = sea_surface.py_get_predominant_heading() # radians
                is_sea_state_same = (new_hs == current_hs) and (new_dp == current_dp)
                # If the sea state has changed then, set the new sea state in the wave glider
                if not is_sea_state_same:
                    sea_surface = py_Sea_surface(new_hs, new_dp, wave_rand_seed, count_wave_spectral_directions, count_wave_spectral_frequencies)
                    asv.py_set_sea_state(sea_surface)            
                # Set rudder angle
                rudder_angle = rudder_controller.get_rudder_angle(asv, py_Coordinates_3D(waypoint_x, waypoint_y))
                # Update thrust tuning factor
                if update_tuning_for_sea_state == True:
                    thrust_tuning_factor = get_thrust_tuning_factor(new_hs, v_zonal, v_meridional, asv.py_get_attitude().z)
                    asv.py_wg_set_thrust_tuning_factor(thrust_tuning_factor)
                # Compute dynamics
                try:
                    asv.py_set_ocean_current(v_zonal, v_meridional)
                    asv.py_wg_compute_dynamics(rudder_angle, time_step_size)
                except Exception as e:
                    return (is_simulation_complete, i, float("NaN"), str(e))
        if time >= kill_time:
            # Simulation incomplete
            return (is_simulation_complete, i, float("NaN"), "Exceeded kill time.")
        else:
            is_simulation_complete = True
            position = asv.py_get_position_cog() # Get the current position in PCS.
            dx = position.x - start_position_PCS.x 
            dy = position.y - start_position_PCS.y 
            distance = math.sqrt(dx*dx + dy*dy) # m
            speed = _MPS_TO_KNOTS * distance/((time - simulation_start_time).total_seconds()) # knots
            return (is_simulation_complete, i, time, speed)
    
    def run_benjamin(self):
        self.df = None # Clear the dataframe so that it will contain only Benjamin's log and not logs from tuning.
        benjamin_log_file = root_dir.joinpath(*"data/pacx/logs/1.1/data/0-data/PacX_NODC_master_folder/benjamin/moseg_benjamin.txt".split("/"))
        self.__load_log_files([benjamin_log_file])
        # Multithread the runs but groupby month for the multithreading.
        results = []
        for group_name, group in self.df.groupby(["year", "month"]):
            #Ref: https://github.com/tqdm/tqdm/discussions/1121
            # Multiprocessing 
            executor = ProcessPoolExecutor()
            year, month = group_name
            rows = group.index
            group.drop(group.index, inplace=True) # Group is a copy of the rows from self.df. Delete the group to save memory.
            month_1 = month 
            month_2 = (month+1) % 12 if (month+1) > 12 else (month+1)
            year_1 = year 
            year_2 = year+1 if month_2<month_1 else year
            file_name_1 = "{}_{}.nc".format(year_1, str(month_1).zfill(2))
            file_name_2 = "{}_{}.nc".format(year_2, str(month_2).zfill(2))
            self.__load_wave_nc_files([file_name_1, file_name_2])
            self.__load_current_nc_files([(year_1, month_1), (year_2, month_2)])
            update_tuning_for_sea_state = True
            tuning_factor = None # Will be set in run() since update_tuning_for_sea_state = True
            params = [(i, self.proximity_margin, update_tuning_for_sea_state, tuning_factor) for i in rows]
            jobs = [executor.submit(self.run, param) for param in params] # Parameters to be passed to method run()
            if self.host_type == "HPC":
                jobs = as_completed(jobs)
            elif self.host_type == "PC":
                # Running on a personal computer, show tqdm progress bar
                jobs = tqdm(as_completed(jobs), total=len(jobs)) # Get a progress bar.
                jobs.set_description("{} {}".format(year, str(month).zfill(2)))
            else:
                raise ValueError("Incorrect value for the variable host_type.")
            for job in jobs:
                results.append(job.result())
        for result in results:
            is_simulation_complete = result[0]
            if is_simulation_complete:
                is_simulation_complete, i, end_time, simulated_speed = result
                self.df.loc[i, "delta_T_simulated(s)"] = (end_time - self.df.loc[i, "Timestamp(UTC)"]).total_seconds() if is_simulation_complete else end_time
                self.df.loc[i, "speed_simulated(knots)"] = simulated_speed
                self.df.loc[i, "error_msg"] = pd.NA
            else:
                self.df.loc[i, "delta_T_simulated(s)"] = pd.NA
                self.df.loc[i, "speed_simulated(knots)"] = pd.NA
                self.df.loc[i, "error_msg"] = result[2]
        # Write the dataframe to file
        benjamin_output_dir = self.output_dir.joinpath("benjamin")    
        benjamin_output_dir.mkdir(parents=True, exist_ok=True)
        benjamin_output_file = benjamin_output_dir.joinpath("simulation_data.csv")
        self.df.to_csv(str(benjamin_output_file))
        print("Benjamin simulation results written to {}".format(benjamin_output_file))

    def _tune_instance(self, params):
        i, tuning_factor, recursion_depth = params
        update_tuning_for_sea_state = False
        params = [i, self.proximity_margin, update_tuning_for_sea_state, tuning_factor] 
        actual_speed = self.df.loc[i, "speed(knots)"]
        result = None 
        try:
            result = self.run(params)
        except Exception as e:
            # Simulation failed
            is_simulation_complete = False
            time = float("NaN")
            simulated_speed = float("NaN")
            tuning_factor = float("NaN")
            error_msg = str(e)
            return (is_simulation_complete, i, time, simulated_speed, tuning_factor, error_msg)
        is_simulation_complete = result[0]
        if is_simulation_complete:
            # Simulation was successful.
            is_simulation_complete, i, time, simulated_speed = result
            speed_ratio = simulated_speed / actual_speed
            if recursion_depth >= 15:
                # Tried too many times and yet cannot find a tuning factor. Stop.
                error_msg = "Recursion limit reached. Final tuning_factor = {}.".format(tuning_factor)
                tuning_factor = float("NaN")
                return (is_simulation_complete, i, time, simulated_speed, tuning_factor, error_msg)
            else:
                if speed_ratio < 0.9 or speed_ratio > 1.1:
                    # Speed ration not within the limit. Require more tuning
                    tuning_factor = tuning_factor / speed_ratio
                    recursion_depth = recursion_depth + 1
                    is_simulation_complete, i, time, simulated_speed, tuning_factor, error_msg = self._tune_instance([i, tuning_factor, recursion_depth])
                    return (is_simulation_complete, i, time, simulated_speed, tuning_factor, error_msg)
                else:
                    return (is_simulation_complete, i, time, simulated_speed, tuning_factor, None)
        else:
            # Simulation failed. 
            is_simulation_complete, i, time, msg = result
            simulated_speed = float("NaN")
            if msg == "Exceeded kill time.":
                # Thrust too small. Scale up the thrust.
                if recursion_depth >= 15:
                    # Tried too many times and yet cannot find a tuning factor. Stop.
                    error_msg = "Recursion limit reached. Final tuning_factor = {}.".format(tuning_factor)
                    tuning_factor = float("NaN")
                    return (is_simulation_complete, i, time, simulated_speed, tuning_factor, error_msg)
                else:
                    tuning_factor = tuning_factor * 3
                    recursion_depth = recursion_depth + 1
                    is_simulation_complete, i, time, simulated_speed, tuning_factor, error_msg = self._tune_instance([i, tuning_factor, recursion_depth])
                    return (is_simulation_complete, i, time, simulated_speed, tuning_factor, error_msg)
            else:
                tuning_factor = float("NaN")
                error_msg = msg
                return (is_simulation_complete, i, time, simulated_speed, tuning_factor, error_msg)
    
    def _run_tuning(self, log_files, out_file_name):
        self.df = None # Clear the dataframe. 
        self.__load_log_files(log_files)
        for month_name, month_group in self.df.groupby(["year", "month"]):
            executor = ProcessPoolExecutor()
            year, month = month_name
            month_1 = month 
            month_2 = (month+1) % 12 if (month+1) > 12 else (month+1)
            year_1 = year 
            year_2 = year+1 if month_2<month_1 else year
            file_name_1 = "{}_{}.nc".format(year, str(month_1).zfill(2))
            file_name_2 = "{}_{}.nc".format(year+1 if month_2<month_1 else year, str(month_2).zfill(2))
            self.__load_wave_nc_files([file_name_1, file_name_2])
            self.__load_current_nc_files([(year_1, month_1), (year_2, month_2)])
            tuning_factor = 1.0
            recursion_depth = 0
            params = [(i, tuning_factor, recursion_depth) for i in month_group.index]
            jobs = [executor.submit(self._tune_instance, param) for param in params] # Parameters to be passed to method run()
            if self.host_type == "HPC":
                jobs = as_completed(jobs)
            elif self.host_type == "PC":
                # Running on a personal computer, show tqdm progress bar
                jobs = tqdm(as_completed(jobs), total=len(jobs)) # Get a progress bar.
                jobs.set_description("{} {}".format(year, str(month).zfill(2)))
            else:
                raise ValueError("Incorrect value for the variable host_type.")
            for job in jobs:
                is_simulation_complete, i, end_time, simulated_speed, tuning_factor, error_msg = job.result()
                self.df.loc[i, "delta_T_simulated(s)"] = (end_time - self.df.loc[i, "Timestamp(UTC)"]).total_seconds() if is_simulation_complete else end_time
                self.df.loc[i, "tuning_factor"] = tuning_factor   
                self.df.loc[i, "error_msg"] = error_msg    
                self.df.loc[i, "speed_simulated(knots)"] = simulated_speed  
            # Write tuning results to file
            tuning_output_dir = self.output_dir.joinpath("tuning")
            tuning_output_dir.mkdir(parents=True, exist_ok=True)
            tuning_output_file = tuning_output_dir.joinpath(out_file_name) 
            self.df.to_csv(str(tuning_output_file))

    def tune_wg_thrust(self):
        # Training data
        training_files = [root_dir.joinpath(*"data/pacx/logs/1.1/data/0-data/PacX_NODC_master_folder/fontaine_maru/moseg_fontaine.txt".split("/")),
                          root_dir.joinpath(*"data/pacx/logs/1.1/data/0-data/PacX_NODC_master_folder/papa_mau/moseg_papamau.txt".split("/")),
                          root_dir.joinpath(*"data/pacx/logs/1.1/data/0-data/PacX_NODC_master_folder/piccard_maru/moseg_piccard.txt".split("/"))]
        self._run_tuning(training_files, "tuning_factors_for_training.csv")
        # Validation data
        validation_file = [root_dir.joinpath(*"data/pacx/logs/1.1/data/0-data/PacX_NODC_master_folder/benjamin/moseg_benjamin.txt".split("/"))]
        self._run_tuning(validation_file, "tuning_factors_for_validation.csv")
