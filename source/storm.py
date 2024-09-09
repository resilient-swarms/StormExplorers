import pathlib
module_dir = pathlib.Path(__file__).parent.resolve()
root_dir = module_dir.parent

import os
import epsg
import pathlib
from shapely.geometry import LineString, Point
from shapely.ops import nearest_points
from netcdf import NetCDF_wave, NetCDF_current, NetCDF_wind, NetCDF_precipitation
import cds
from storms_archive import Storms_archieve
from datetime import datetime

class Storm:
    '''
    Class to simulate waves in the storm and the track of the storm. 
    :param year: the year of the storm.
    :param storm_name: name of the storm.
    :param north: north limit of the simulated sea surface.
    :param south: south limit of the simulated sea surface.
    :param east: east limit of the simulated sea surface.
    :param west: west limit of the simulated sea surface.
    '''
    def __init__(self, year, storm_name, north, south, east, west, init_ocean_current=False):
        '''
        Init the storm. Check if the storms database is available and creates it if necessary.
        From the storms database extract actual and predicted paths of for the storm. 
        Also, checks if the wave data is available and downloads if necessary. 
        '''
        self.map_boundary_north = north
        self.map_boundary_south = south
        self.map_boundary_east = east
        self.map_boundary_west = west
        # Check if the storms database exist and create it if necessary.
        db_file = root_dir.joinpath("data", "storms")
        storms = Storms_archieve(db_file) # db_file will be created if not existing.
        # Find the storm in the db
        self.storms_df = storms.get_storms()
        storm_index = self.storms_df.loc[(self.storms_df["year"]==year) & (self.storms_df["storm_name"]==storm_name)].index[0]
        # Download and load the wave and current data
        self.nc_data_wave = {}
        self.nc_data_wind = {}
        self.nc_data_precipitation = {}
        self.nc_data_current = {}
        self._download_and_load_wave_data(storm_index)
        self._download_and_load_wind_data(storm_index)
        self._download_and_load_precipitation_data(storm_index)
        if init_ocean_current:
            self._download_and_load_ocean_current_data(storm_index)
        # Copy storm tracks
        self._time_stamps = self.storms_df.loc[storm_index, "time_stamps"]
        self._actual_track_GCS = self.storms_df.loc[storm_index, "path(long,lat)"]
        self._actual_track_PCS = []
        for point in self._actual_track_GCS:
            point_x, point_y = epsg.GCS_to_PCS(*point)
            self._actual_track_PCS.append((point_x, point_y))
        self._predicted_tracks_GCS = self.storms_df.loc[storm_index, "predicted_paths(long,lat)"]
        self._predicted_tracks_PCS = []
        for path in self._predicted_tracks_GCS:
            path_PCS = []
            for point in path:
                point_x, point_y = epsg.GCS_to_PCS(*point)
                path_PCS.append((point_x, point_y))
            self._predicted_tracks_PCS.append(path_PCS)

    def _download_and_load_wave_data(self, storm_index):
        nc_files = []
        time_stamps = self.storms_df.loc[storm_index, "time_stamps"]
        for j in (0,-1):
            year  = time_stamps[j].year
            month = time_stamps[j].month
            nc_files.append("{}_{}.nc".format(year, str(month).zfill(2)))
        nc_files= set(nc_files) # Make a unique list of the required nc files
        # Check if the file exist, else download it.
        wave_dir = root_dir.joinpath("data", "cds", "north_atlantic", "waves")
        wave_dir.mkdir(parents=True, exist_ok=True) # Make the dir if it does not exist.
        for file in nc_files:
            nc_file = wave_dir.joinpath(file)
            if not os.path.exists(str(nc_file)):
                print(nc_file)
                # nc file does not exist, download it.
                year, month = file[:-3].split("_")
                cds.get_wave_data(int(year), int(month), 
                                  self.map_boundary_north, 
                                  self.map_boundary_south, 
                                  self.map_boundary_east, 
                                  self.map_boundary_west, 
                                  str(nc_file))
            else:
                print("Found {}".format(file))
            # Load the file
            self.nc_data_wave[file] = NetCDF_wave(str(nc_file))
            print("Loaded wave data file {}".format(file))

    def _download_and_load_wind_data(self, storm_index):
        nc_files = []
        time_stamps = self.storms_df.loc[storm_index, "time_stamps"]
        for j in (0,-1):
            year  = time_stamps[j].year
            month = time_stamps[j].month
            nc_files.append("{}_{}.nc".format(year, str(month).zfill(2)))
        nc_files= set(nc_files) # Make a unique list of the required nc files
        # Check if the file exist, else download it.
        wave_dir = root_dir.joinpath("data", "cds", "north_atlantic", "winds")
        wave_dir.mkdir(parents=True, exist_ok=True) # Make the dir if it does not exist.
        for file in nc_files:
            nc_file = wave_dir.joinpath(file)
            if not os.path.exists(str(nc_file)):
                print(nc_file)
                # nc file does not exist, download it.
                year, month = file[:-3].split("_")
                cds.get_wind_data(int(year), int(month),
                                  self.map_boundary_north,
                                  self.map_boundary_south,
                                  self.map_boundary_east,
                                  self.map_boundary_west,
                                  str(nc_file))
            else:
                print("Found {}".format(file))
            # Load the file
            self.nc_data_wind[file] = NetCDF_wind(str(nc_file))
            print("Loaded wind data file {}".format(file))

    def _download_and_load_precipitation_data(self, storm_index):
        nc_files = []
        time_stamps = self.storms_df.loc[storm_index, "time_stamps"]
        for j in (0,-1):
            year  = time_stamps[j].year
            month = time_stamps[j].month
            nc_files.append("{}_{}.nc".format(year, str(month).zfill(2)))
        nc_files= set(nc_files) # Make a unique list of the required nc files
        # Check if the file exist, else download it.
        wave_dir = root_dir.joinpath("data", "cds", "north_atlantic", "precipitation")
        wave_dir.mkdir(parents=True, exist_ok=True) # Make the dir if it does not exist.
        for file in nc_files:
            nc_file = wave_dir.joinpath(file)
            if not os.path.exists(str(nc_file)):
                print(nc_file)
                # nc file does not exist, download it.
                year, month = file[:-3].split("_")
                cds.get_precipitation_data(int(year), int(month),
                                  self.map_boundary_north,
                                  self.map_boundary_south,
                                  self.map_boundary_east,
                                  self.map_boundary_west,
                                  str(nc_file))
            else:
                print("Found {}".format(file))
            # Load the file
            self.nc_data_precipitation[file] = NetCDF_precipitation(str(nc_file))
            print("Loaded precipitation data file {}".format(file))

    def _download_and_load_ocean_current_data(self, storm_index):
        zip_files = []
        time_stamps = self.storms_df.loc[storm_index, "time_stamps"]
        for j in (0,-1):
            year  = time_stamps[j].year
            month = time_stamps[j].month
            zip_files.append("{}_{}.zip".format(year, str(month).zfill(2)))
        zip_files = set(zip_files) # Make a unique list of the required nc files
        # Check if the file exist, else download it.
        ocean_current_dir = root_dir.joinpath("data", "cds", "north_atlantic", "ocean_currents")
        ocean_current_dir.mkdir(parents=True, exist_ok=True) # Make the dir if it does not exist.
        for file in zip_files:
            zip_file = ocean_current_dir.joinpath(file)
            year, month = file[:-4].split("_")
            if not os.path.exists(str(zip_file)):
                # tar file does not exist, download it.
                cds.get_ocean_current_data(int(year), int(month), str(zip_file))
            else:
                print("Found {}".format(file))
            # Load the file
            days_loaded = []
            days_not_loaded = []
            for day in range(1, 32):
                file_name = "{}_{}_{}.nc".format(year, str(month).zfill(2), str(day).zfill(2))
                nc_file = ocean_current_dir.joinpath(file_name)
                if os.path.exists(nc_file):
                    self.nc_data_current[file_name] = NetCDF_current(str(nc_file))
                    days_loaded.append(day)
                else:
                    days_not_loaded.append(day)
            print("Loaded ocean current data files for the period {}-{}-[{} to {}]{}".format(year, str(month).zfill(2), days_loaded[0], days_loaded[-1], 
                  "; could not find files for the days {}".format(days_not_loaded) if len(days_not_loaded) != 0 else ""))
    
    def get_wave_data_at(self, longitude, latitude, time):
        '''
        Returns the significant wave height and mean wave direction for a given location and time.
        '''
        year = time.year
        month = time.month
        file_name = "{}_{}.nc".format(year, str(month).zfill(2))
        wave_hs, wave_dp = self.nc_data_wave[file_name].get_wave_data_at(longitude, latitude, time)
        return (wave_hs, wave_dp)

    def get_wind_speed_at(self, longitude, latitude, time):
        '''
        Returns the eastward and northward wind velocities for a given location and time.
        '''
        year = time.year
        month = time.month
        file_name = "{}_{}.nc".format(year, str(month).zfill(2))
        u, v = self.nc_data_wind[file_name].get_wind_speed_at(longitude, latitude, time)
        return (u, v)

    def get_precipitation_at(self, longitude, latitude, time):
        '''
        Returns the total precipitation for a given location and time.
        '''
        year = time.year
        month = time.month
        file_name = "{}_{}.nc".format(year, str(month).zfill(2))
        tp = self.nc_data_precipitation[file_name].get_precipitation_at(longitude, latitude, time)
        return tp

    def get_ocean_current_at(self, longitude, latitude, time):
        year = time.year
        month = time.month
        day = time.day
        file_name = "{}_{}_{}.nc".format(year, str(month).zfill(2), str(day).zfill(2))
        v_zonal, v_meridional = self.nc_data_current[file_name].get_ocean_current_at(longitude, latitude, time)
        if v_zonal == None or v_meridional == None:
            print("Ocean current data not available for ({}, {}) at time {} in file {}".format(longitude, latitude, time, file_name))
            v_zonal, v_meridional = (0, 0)
        return (v_zonal, v_meridional)
    
    def get_eye_location(self, date_time):
        '''
        Returns the longitude and latitude of the eye of the storm for a given datetime instance. 
        '''
        time_diffs = [abs((time - date_time).total_seconds()) for time in self._time_stamps]
        min_diff = min(time_diffs)
        for i in range(len(time_diffs)):
            if time_diffs[i] == min_diff:
                return self._actual_track_GCS[i]
            else:
                continue 

    def get_nearest_point_on_predicted_track(self, longitude, latitude, date_time):
        '''
        Returns the longitude and latitude of the nearest point from the predicted path of the storm for the given datetime instance.  
        '''
        time_diffs = [abs((time - date_time).total_seconds()) for time in self._time_stamps]
        min_diff = min(time_diffs)
        for i in range(len(time_diffs)):
            if time_diffs[i] == min_diff:
                index = i 
                break
            else:
                continue 
        predicted_track_PCS = self._predicted_tracks_PCS[index]
        # Find the nearest point on the predicted track
        location_longitude, location_latitude = epsg.GCS_to_PCS(longitude, latitude)
        location = Point(location_longitude, location_latitude)
        path = LineString(predicted_track_PCS)
        nearest_point = nearest_points(path, location)[0]
        nearest_point_longitude, nearest_point_latitude = epsg.PCS_to_GCS(nearest_point.x, nearest_point.y)
        return (nearest_point_longitude, nearest_point_latitude)

    def get_storm_track_bounding_box(self):
        '''
        Returns the north, south, east, west limits of the storm track. 
        '''
        longitudes = []
        latitudes = []
        for point in self._actual_track_GCS:
            longitude, latitude = point
            longitudes.append(longitude)
            latitudes.append(latitude)
        north = max(latitudes)
        south = min(latitudes)
        east = max(longitudes)
        west = min(longitudes)
        return (north, south, east, west)