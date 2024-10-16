import pathlib
import sys
module_dir = pathlib.Path(__file__).parent.resolve()
root_dir = module_dir.parent
asvlite_wrapper_dir = root_dir.joinpath("dependency", "ASVLite", "wrapper", "cython")
sys.path.insert(0, str(asvlite_wrapper_dir))

import os
import math
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import multiprocessing as mp
from datetime import datetime, timedelta
from shapely.geometry import Polygon, LineString, Point
from tqdm import tqdm
from sklearn.neighbors import KNeighborsClassifier

import epsg
import cds
from storms_archive import Storms_archieve
from rudder_controller import Rudder_PID_controller
from netcdf import NetCDF_wave, NetCDF_current
from sea_surface import py_Sea_surface
from asv import py_Asv_specification, py_Asv
from geometry import py_Coordinates_3D, py_normalise_angle_PI, py_normalise_angle_2PI 
from thrust_calibrator import Thrust_calibrator
from storm import Storm

_MPS_TO_KNOTS = 1.94384449 # 1 m/s = 1.94384449 knots
_SIMULATION_TIME_STEP_SIZE = 40 # milli-sec

class Storm_tracking:
    def __init__(self, host_type) -> None:
        '''
        Valid values for host_type are the strings - "PC", "HPC"
        '''
        self.host_type = host_type
        # Map boundary
        self.map_boundary_south = 15
        self.map_boundary_north = 50
        self.map_boundary_west = -100
        self.map_boundary_east = -45
        # Deployment region
        self.deployment_boundary_south = 16.5
        self.deployment_boundary_north = 21.5 
        self.deployment_boundary_west  = -86.0
        self.deployment_boundary_east  = -81.0
        x0, y0 = epsg.GCS_to_PCS(self.deployment_boundary_west, self.deployment_boundary_south)
        x1, y1 = epsg.GCS_to_PCS(self.deployment_boundary_east, self.deployment_boundary_north)
        self.deployment_region = Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
        # Deployment configuration
        self.swarm_size = 9
        self.deployment_names = []
        self.deployments = [] # array of array of Points.
        diagonal_NW_SE_geom = LineString([(x0, y1), (x1, y0)]) 
        diagonal_SW_NE_geom = LineString([(x0, y0), (x1, y1)]) 
        diagonal_length = diagonal_NW_SE_geom.length   # m
        swarm_spacing = diagonal_length/(self.swarm_size+1) # m

        # Config 0 - baseline MAX
        points = ((-86.0, 20.0), (-86.0, 19.0), (-85.0, 19.0), (-86.0, 18.0), (-85.0, 18.0), (-84.0, 18.0), (-85.0, 17.0), (-84.0, 17.0), (-83.0, 17.0))
        max_deployment = []
        for point in points:
            long, lat = point
            x, y = epsg.GCS_to_PCS(long, lat)
            max_deployment.append(Point(x, y))
        self.deployment_names.append("MAX")
        self.deployments.append(max_deployment)
        # Config 1 - baseline MIN
        points = ((-85.5, 21.5), (-84.5, 21.5), (-83.5, 21.5), (-84.5, 20.5), (-83.5, 20.5), (-82.5, 20.5), (-83.5, 19.5), (-82.5, 19.5), (-81.5, 18.5))
        min_deployment = []
        for point in points:
            long, lat = point
            x, y = epsg.GCS_to_PCS(long, lat)
            min_deployment.append(Point(x, y))
        self.deployment_names.append("MIN")
        self.deployments.append(min_deployment)

        # Config 2 - NW-SE
        NW_SE_deployment = []
        for i in range(1, self.swarm_size+1):
            dist = i*swarm_spacing
            p = diagonal_NW_SE_geom.interpolate(dist, normalized=False)
            NW_SE_deployment.append(p)
        self.deployment_names.append("NW-SE")
        self.deployments.append(NW_SE_deployment)
        # Config 3 - SW-NE
        SW_NE_deployment = []
        for i in range(1, self.swarm_size+1):
            dist = i*swarm_spacing
            p = diagonal_SW_NE_geom.interpolate(dist, normalized=False)
            SW_NE_deployment.append(p)
        self.deployment_names.append("SW-NE")
        self.deployments.append(SW_NE_deployment)
        # Config 4 - SQU
        p0 = SW_NE_deployment[0]
        p1 = NW_SE_deployment[-1]
        p2 = SW_NE_deployment[-1]
        p3 = NW_SE_deployment[0]
        p01 = Point((p0.x+p1.x)/2, (p0.y+p1.y)/2)
        p03 = Point((p0.x+p3.x)/2, (p0.y+p3.y)/2)
        p23 = Point((p2.x+p3.x)/2, (p2.y+p3.y)/2)
        p12 = Point((p1.x+p2.x)/2, (p1.y+p2.y)/2)
        p_c = Point((p03.x+p12.x)/2, (p03.y+p12.y)/2)
        square_deployment= [p0, p01, p1, p03, p_c, p12, p3, p23, p2]
        self.deployment_names.append("SQU")
        self.deployments.append(square_deployment)
        
        # Set the list of storms that pass through the deployment region
        self._set_storms()
        # Messages
        self.warning_start_time_not_found = "Couldn't find the predicted storm track that intersected the deployment region. Therefore, setting the start time for simulation equal to the first recorded time for the storm.\n"
        self.warning_end_time_not_found = "Couldn't find the the time when the eye of the storm exited the deployment region. Therefore, setting the end time for simulation equal to the last recorded time for the storm.\n"
    
    def _set_storms(self):
        self.storms = pd.DataFrame()
        storms_archieve = Storms_archieve("../data/storms")
        df_storms = storms_archieve.get_storms()
        for i in range(len(df_storms)):
            path = df_storms.loc[i, "path(long,lat)"]
            path = [epsg.GCS_to_PCS(long, lat) for (long, lat) in path]
            path = LineString(path)
            # Check for intersection of path with deployment_region
            intersection = self.deployment_region.intersection(path)
            if not intersection.is_empty:
                if self.storms.empty:
                    self.storms = df_storms.loc[i].to_frame().T
                else:
                    self.storms = pd.concat([self.storms, df_storms.loc[i].to_frame().T], ignore_index=True)
       
    def _simulate_wg(self, start_position, start_time, end_time, storm, output_file, progress_bar, return_dict):
        is_simulation_complete = False
        f = open(output_file, "w")  
        f.write("time_stamp,longitude,latitude,z,wg_heading,hs(m),dist(km)\n")   
        wg_name = str(output_file).split("/")[-1]     
        # Wave glider specs
        asv_spec = py_Asv_specification(
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
        start_position_GCS = py_Coordinates_3D(*start_position) # NOTE: start_position is the tuple (longitude, latitude)
        x, y = epsg.GCS_to_PCS(start_position_GCS.x, start_position_GCS.y)
        start_position_PCS = py_Coordinates_3D(x, y)
        start_attitude = py_Coordinates_3D(0.0, 0.0, math.pi/2) # All wave gliders start with a heading towards East.
        # Sea surface
        wave_rand_seed = 1
        count_component_waves = 21
        wave_hs, wave_dp        = storm.get_wave_data_at(start_position_GCS.x, start_position_GCS.y, start_time) 
        v_zonal, v_meridional   = storm.get_ocean_current_at(start_position_GCS.x, start_position_GCS.y, start_time)
        if wave_hs == None or wave_hs == 0.0 or wave_dp == None:
            return_dict[wg_name] = (is_simulation_complete, "Wave data not available.")
        sea_surface = py_Sea_surface(wave_hs, wave_dp, wave_rand_seed, count_component_waves)
        # Initialise the wave glider
        asv = py_Asv(asv_spec, sea_surface, start_position_PCS, start_attitude)
        # Initialise the rudder controller
        rudder_controller = Rudder_PID_controller(asv_spec)
        rudder_angle = 0.0  
        # Initialise the thrust calibrator
        thrust_calibrator = Thrust_calibrator(self.host_type)
        # Initialise time to start of simulation
        time = start_time - timedelta(seconds=_SIMULATION_TIME_STEP_SIZE/1000.0)  
        # Simulate till last waypoint
        while time <= end_time: 
            time = time + timedelta(seconds=_SIMULATION_TIME_STEP_SIZE/1000.0)
            # Update waypoint 
            position = asv.py_get_position_cog() # Get the current position in PCS.
            wg_longitude, wg_latitude = epsg.PCS_to_GCS(position.x, position.y)
            waypoint = storm.get_nearest_point_on_predicted_track(wg_longitude, wg_latitude, time)
            waypoint = py_Coordinates_3D(*waypoint) 
            waypoint_x, waypoint_y = epsg.GCS_to_PCS(waypoint.x, waypoint.y)  
            # Get the sea state
            new_hs, new_dp = storm.get_wave_data_at(wg_longitude, wg_latitude, time) 
            if new_hs == None or new_hs == 0.0 or new_dp == None:
                return_dict[wg_name] = (is_simulation_complete, "Wave data not available.")
            v_zonal, v_meridional = storm.get_ocean_current_at(wg_longitude, wg_latitude, time)
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
                sea_surface = py_Sea_surface(new_hs, new_dp, wave_rand_seed, count_component_waves)
                asv.py_set_sea_state(sea_surface)            
            # Set rudder angle
            rudder_angle = rudder_controller.get_rudder_angle(asv, py_Coordinates_3D(waypoint_x, waypoint_y))
            # Update thrust tuning factor
            thrust_tuning_factor = thrust_calibrator.get_thrust_tuning_factor(new_hs, v_zonal, v_meridional, asv.py_get_attitude().z)
            asv.py_wg_set_thrust_tuning_factor(thrust_tuning_factor)
            # Compute dynamics
            try:
                asv.py_set_ocean_current(v_zonal, v_meridional)
                asv.py_wg_compute_dynamics(rudder_angle, _SIMULATION_TIME_STEP_SIZE)
            except Exception as e:
                return_dict[wg_name] = (is_simulation_complete, str(e))
            # Get the distance to the centre of the storm
            storm_longitude, storm_latitude = storm.get_eye_location(time)
            storm_x, storm_y = epsg.GCS_to_PCS(storm_longitude, storm_latitude)
            position = asv.py_get_position_cog() # Get the current position in PCS.
            wg_longitude, wg_latitude = epsg.PCS_to_GCS(position.x, position.y)
            attitude = asv.py_get_attitude()
            distance_to_storm = Point(storm_x, storm_y).distance(Point(position.x, position.y)) / 1000.0 # Km
            # Save simulation data
            if time.second == 0 and time.microsecond == 0:
                f.write("{date},{x},{y},{z},{heading},{hs},{dist}\n".format( 
                                                        date=time.strftime("%Y-%m-%d %H:%M:%S.%f"),
                                                        x=wg_longitude, 
                                                        y=wg_latitude, 
                                                        z=position.z, 
                                                        heading=attitude.z, 
                                                        hs=current_hs,
                                                        dist=distance_to_storm))
            # Update progress bar
            if progress_bar is not None:
                progress_bar.update(1)
        is_simulation_complete = True
        msg = None
        if progress_bar is not None:
            progress_bar.close()
        f.close()
        return_dict[wg_name] = (is_simulation_complete, msg) 
    
    def _simulate_swarm_in_storm(self, storm_index):
        # Storm details
        storm_name = self.storms.loc[storm_index, "storm_name"]
        storm_year = self.storms.loc[storm_index, "year"]
        # Output directory for the storm
        storm_dir = root_dir.joinpath(*"results/north_atlantic_deployments/{}_{}".format(int(storm_year), storm_name).split("/"))
        storm_dir.mkdir(parents=True, exist_ok=True) # Make the dir if it does not exist.
        log_file_name = storm_dir.joinpath("log.txt")
        log_file = open(log_file_name, "w")
        # Set the start time of simulation.
        # From the list of predicted tracks find the first predicted track that intersects the deployment region.
        start_time_index = None
        simulation_start_time = None 
        predicted_paths = self.storms.loc[storm_index, "predicted_paths(long,lat)"]
        for i in range(len(predicted_paths)):
            predicted_path = LineString([epsg.GCS_to_PCS(long, lat) for (long, lat) in predicted_paths[i]])
            intersection = predicted_path.intersection(self.deployment_region)
            if not intersection.is_empty:
                start_time_index = i
                break 
        if start_time_index == None:
            start_time_index = 0
            log_file.write("[warning] {}".format(self.warning_start_time_not_found))
        simulation_start_time = self.storms.loc[storm_index, "time_stamps"][start_time_index]
        # Set the end time of simulation.
        # Stop the simulation when the eye of the storm is outside of the deployment region after entering it.
        end_time_index = None
        simulation_end_time = None
        storm_path = self.storms.loc[storm_index, "path(long,lat)"]
        has_entered_the_region = False 
        for i in range(len(storm_path)):
            long, lat = storm_path[i]
            point = Point(epsg.GCS_to_PCS(long, lat))
            if self.deployment_region.contains(point):
                has_entered_the_region = True 
            if has_entered_the_region:
                if not self.deployment_region.contains(point):
                    end_time_index = i
                    break
        if end_time_index == None:
            end_time_index = -1
            log_file.write("[warning] {}".format(self.warning_end_time_not_found))
        simulation_end_time = self.storms.loc[storm_index, "time_stamps"][end_time_index]
        # Init the storm
        global storm
        init_ocean_currents = True
        storm = Storm(storm_year, 
                    storm_name, 
                    self.map_boundary_north, 
                    self.map_boundary_south,
                    self.map_boundary_east,
                    self.map_boundary_west,
                    init_ocean_currents)
        # Plot the storm track
        self._plot_storm_track(storm_index, storm_dir)
        # Simulate deployments
        for i in range(len(self.deployments)):
            print("\n{} deployment config {}.".format(storm_name, i))
            deployment_dir = storm_dir.joinpath("config_{}".format(i))
            deployment_dir.mkdir(parents=True, exist_ok=True) # Make the dir if it does not exist.
            processes = [] 
            manager = mp.Manager()
            return_dict = manager.dict()
            progress_bars = []
            # Create the processes
            for j in range(self.swarm_size):
                count_simulation_steps = ((simulation_end_time - simulation_start_time).total_seconds()*1000)/_SIMULATION_TIME_STEP_SIZE
                progress_bar = None 
                if self.host_type == "PC":
                    progress_bar = tqdm(leave=False, total=count_simulation_steps, ncols=120)
                    progress_bar.set_description("wg{}".format(j))
                    progress_bars.append(progress_bar)
                elif self.host_type == "HPC":
                    # Don't set the tqdm bars if running on a cluster.
                    progress_bar = None 
                else:
                    raise ValueError("Incorrect value for the variable host_type.")
                # Output file for wave glider
                wg_output_file = deployment_dir.joinpath("wg{}".format(j))
                # Start position 
                deployment_point = self.deployments[i][j]
                longitude, latitude = epsg.PCS_to_GCS(deployment_point.x, deployment_point.y)
                wg_start_position = (longitude, latitude)
                args = [wg_start_position, 
                        simulation_start_time, 
                        simulation_end_time, 
                        storm,
                        wg_output_file,
                        progress_bar,
                        return_dict]
                process = mp.Process(target=self._simulate_wg, args=args)
                processes.append(process)
            # Start the simulations
            for process in processes:
                process.start()
            # Close processes
            for process in processes:
                process.join()
            # Log if simulation faced errors
            if(len(return_dict) != 0):
                for key, value in return_dict.items():
                    is_simulation_complete, message = value 
                    if not is_simulation_complete:
                        log_file.write("[error] Deployment config {}, {}: {}.\n".format(i, key, value))
            # Clear
            progress_bars.clear()
            processes.clear()
            return_dict.clear()
            # Plot the deployment result
            self._plot_deployment(storm_index, start_time_index, end_time_index, deployment_dir)
            # Compute deployment performance
            self._compute_deployment_performance(storm_index, start_time_index, end_time_index, deployment_dir)
        # Close the log file
        log_file.close()
    
    def _compute_deployment_performance(self, storm_index, start_time_index, end_time_index, deployment_dir):
        df = pd.DataFrame()
        for i in range(self.swarm_size):
            wg_output_file = deployment_dir.joinpath("wg{}".format(i))
            if os.path.exists(wg_output_file):
                df_wg = pd.read_csv(wg_output_file)
                df_wg.set_index("time_stamp")
                if i == 0:
                    df["time_stamp"] = df_wg["time_stamp"]
                    df.set_index("time_stamp")
                df["wg{}".format(i)] = df_wg["dist(km)"]
        df["min_dist(km)"] = df[["wg{}".format(i) for i in range(self.swarm_size)]].min(axis=1)
        output_file = deployment_dir.joinpath("performance.csv")
        df.to_csv(str(output_file))
    
    def _plot_storm_track(self, storm_index, storm_dir):
        # Plot the grids
        ax = plt.axes(projection=ccrs.PlateCarree())
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        # Plot map
        ax.coastlines(color='grey', linewidth=0.5)
        # Plot the boundaries of the deployment region
        longitudes = []
        latitudes = []
        edges = self.deployment_region.boundary
        for point in edges.coords:
            longitude, latitude = epsg.PCS_to_GCS(*point)
            longitudes.append(longitude)
            latitudes.append(latitude)
        plt.plot(longitudes, latitudes, color='red', linestyle='-', linewidth=1, transform=ccrs.Geodetic())
        # Plot storm track
        positions = self.storms.loc[storm_index, "path(long,lat)"]
        longitudes = [position[0] for position in positions]
        latitudes  = [position[1] for position in positions]
        plt.plot(longitudes, latitudes, color='black', linewidth=0.5, transform=ccrs.Geodetic(), label="Storm track")
        plot_path_file = storm_dir.joinpath("tc_path.png")
        plt.savefig(str(plot_path_file), dpi=300)
        plt.close() 
        plt.cla() # clear axis 
        plt.clf() # clear figure

    def _plot_deployment(self, storm_index, start_time_index, end_time_index, deployment_dir):
        # Plot the grids
        ax = plt.axes(projection=ccrs.PlateCarree())
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        # Plot map
        ax.coastlines(color='grey', linewidth=0.5)
        # Plot the boundaries of the deployment region
        longitudes = []
        latitudes = []
        edges = self.deployment_region.boundary
        for point in edges.coords:
            longitude, latitude = epsg.PCS_to_GCS(*point)
            longitudes.append(longitude)
            latitudes.append(latitude)
        plt.plot(longitudes, latitudes, color='red', linestyle='-', linewidth=1, transform=ccrs.Geodetic())
        # Plot storm track
        positions = self.storms.loc[storm_index, "path(long,lat)"][:end_time_index+1]
        longitudes = [position[0] for position in positions]
        latitudes  = [position[1] for position in positions]
        plt.plot(longitudes, latitudes, color='black', linewidth=0.5, transform=ccrs.Geodetic(), label="Storm track")
        # Create time markers
        marker_times = self.storms.loc[storm_index, "time_stamps"][:end_time_index+1:3]
        # Set time markers on storm track
        positions = self.storms.loc[storm_index, "path(long,lat)"][:end_time_index+1:3]
        longitudes = [position[0] for position in positions]
        latitudes  = [position[1] for position in positions]
        # Plot wave glider path
        for i in range(self.swarm_size):
            wg_output_file = deployment_dir.joinpath("wg{}".format(i))
            if os.path.exists(wg_output_file):
                df = pd.read_csv(wg_output_file) 
                if not df.empty:
                    df["time_stamp"] = pd.to_datetime(df["time_stamp"]) # Convert column Timestamp to datetime format.
                    # Ploting all the points from the file can be too much,therefor plot the position 
                    # of the vehicle for each min instead of each simulation time step
                    # simulation_time_step_size = _SIMULATION_TIME_STEP_SIZE/1000.0 # sec
                    # n = int(1/simulation_time_step_size * 60) 
                    n = 1
                    times = df["time_stamp"][::n]
                    longitudes = df["longitude"][::n]
                    latitudes = df["latitude"][::n]
                    wave_glider_path_plot = plt.plot(longitudes, latitudes, linestyle='-', linewidth=0.5, transform=ccrs.Geodetic())
                    # Add marker for the wave glider ids
                    if(len(longitudes) != 0 and len(latitudes) != 0):
                        plt.text(longitudes[0]-0.1, latitudes[0], "wg{}".format(i), fontsize=4, ha='right', va='top')
                    # Create markers for wg path
                    color = wave_glider_path_plot[0].get_color()
                    indices = []
                    for time in marker_times:
                        index = df['time_stamp'].sub(time).abs().idxmin()
                        indices.append(index)
                    for j in range(len(indices)):
                        if indices[j] != 0:
                            plt.text(df["longitude"][indices[j]], df["latitude"][indices[j]], str(j), fontsize=4, color=color)
        plot_path_file = deployment_dir.joinpath("wg_paths.png")
        plt.savefig(str(plot_path_file), dpi=300)
        plt.close() 
        plt.cla() # clear axis 
        plt.clf() # clear figure
    
    def run_all(self):
        count_storms_executed_in_parallel = math.ceil(mp.cpu_count() / self.swarm_size)
        i = 0
        while i < len(self.storms):
            processes = [] 
            for j in range(count_storms_executed_in_parallel):
                process = mp.Process(target=self._simulate_swarm_in_storm, args=[i])
                i += 1
                processes.append(process)
            # Start simulation of each storm
            for process in processes:
                process.start()
            # Close the processes
            for process in processes:
                process.join()

    def compute_performance(self):
        df = pd.DataFrame(columns=["storm_name", "year", "config", "min_dist(km)", "normalised_min_dist", "avg_dist(km)", "normalised_avg_dist"])
        results_dir = root_dir.joinpath("results", "north_atlantic_deployments")
        summary_plot = results_dir.joinpath("summary.png")
        # Compute min distance and avg distance
        for storm_index in range(len(self.storms)):
            storm_name = self.storms.loc[storm_index, "storm_name"]
            storm_year = self.storms.loc[storm_index, "year"]
            storm_dir = results_dir.joinpath("{}_{}".format(int(storm_year), storm_name))
            for i in range(len(self.deployments)):
                config_dir = storm_dir.joinpath("config_{}".format(i))
                file_name = config_dir.joinpath("performance.csv")
                peformance_df = None
                try:
                    peformance_df = pd.read_csv(str(file_name))
                    min_dist = peformance_df["min_dist(km)"].min()
                    avg_dist = peformance_df["min_dist(km)"].mean()
                    # Update df
                    new_row = pd.DataFrame({"storm_name":storm_name, "year":storm_year, "config":i, "min_dist(km)":min_dist, "avg_dist(km)":avg_dist}, index=[0]) 
                    df = pd.concat([df, new_row], ignore_index=True)
                except FileNotFoundError as e:
                    pass
        # Compute normalised min and avg distance.
        for storm_name, group in df.groupby(["storm_name", "year"]):
            # Normalised min_dist
            avg = group["min_dist(km)"].mean()
            std = group["min_dist(km)"].std()
            max = group["min_dist(km)"].max()
            for i in group.index:
                normalised_min_dist = 0.0
                if std != 0:
                    normalised_min_dist = (df.loc[i, "min_dist(km)"] - avg)/max 
                df.loc[i, "normalised_min_dist"] = normalised_min_dist
            # Normalised avg_dist
            avg = group["avg_dist(km)"].mean()
            std = group["avg_dist(km)"].std()
            max = group["avg_dist(km)"].max()
            for i in group.index:
                normalised_avg_dist = (df.loc[i, "avg_dist(km)"] - avg)/max 
                df.loc[i, "normalised_avg_dist"] = normalised_avg_dist
        # Cluster the normalised results for plotting.
        X = []
        y = []
        for config, group in df.groupby("config"):
            for i in group.index:
                xx = [group.loc[i, "normalised_min_dist"], group.loc[i, "normalised_avg_dist"]]
                yy = group.loc[i, "config"]
                X.append(xx)
                y.append(yy)
        X = np.array(X)
        y = np.array(y)
        knn = KNeighborsClassifier(n_neighbors=75, weights="distance") # Note: experiment with different values for n_neighbors
        knn.fit(X, y)
        resolution = 1000
        mins = X.min(axis=0) - 0.1
        maxs = X.max(axis=0) + 0.1
        xx, yy = np.meshgrid(np.linspace(mins[0], maxs[0], resolution),
                         np.linspace(mins[1], maxs[1], resolution))
        Z = knn.predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)
        plt.contourf(Z, extent=(mins[0], maxs[0], mins[1], maxs[1]), cmap="Pastel2", alpha=0.5)
        plt.contour(Z, extent=(mins[0], maxs[0], mins[1], maxs[1]), linewidths=0.25, colors='grey')
        # Scatter plot
        colors = ["green", "red", "black", "blue", "violet"]
        labels = self.deployment_names
        legend_elements = []
        for color, label in zip(colors, labels):
            legend_elements.append(Patch(color=color, label=label))
        for config, group in df.groupby("config"):
            config = int(config)
            plt.plot(group["normalised_min_dist"].to_list(), group["normalised_avg_dist"].to_list(), 'k.', color=colors[config], markersize=5)
        # plt.axhline(0, linewidth=0.5, color="grey")
        # plt.axvline(0, linewidth=0.5, color="grey")
        plt.xlabel("Normalised minimum distance", fontsize=8)
        plt.ylabel("Normalised average distance", fontsize=8)
        plt.xticks(fontsize=4) # set the new x-axis labels
        plt.yticks(fontsize=4) # set the new y-axis label size
        plt.legend(handles=legend_elements, prop={'size': 6})
        plt.savefig(str(summary_plot), bbox_inches='tight', dpi=600)
        # plt.show()
        plt.close()
        plt.cla() # clear axis 
        plt.clf() # clear figure    

if __name__ == '__main__':   
    storm_tracking = Storm_tracking("HPC")
    # storm_tracking.run_all()
    storm_tracking.compute_performance()