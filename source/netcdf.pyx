import netCDF4 as nc 
from datetime import datetime, timedelta
import numpy
import math

cdef class NetCDF_wave:
    '''
    Class to extract wave data from an .nc file.
    '''
    def __cinit__(self, file_path) -> None:
        data = nc.Dataset(file_path)
        self._epoch_time = datetime(1900, 1, 1)
        # Longitude range
        self._first_longitude = data.variables["longitude"][0].item()
        self._last_longitude = data.variables["longitude"][-1].item()
        self._delta_longitude = data.variables["longitude"][1].item() - data.variables["longitude"][0].item()
        # Latitude range
        self._first_latitude = data.variables["latitude"][0].item()
        self._last_latitude = data.variables["latitude"][-1].item()
        self._delta_latitude = data.variables["latitude"][1].item() - data.variables["latitude"][0].item()
        # Time range
        self._first_time = data.variables["time"][0].item()
        self._last_time = data.variables["time"][-1].item()
        self._delta_time = data.variables["time"][1].item() - data.variables["time"][0].item()
        # Numpy array of significant wave hts and wave headings. 
        # These arrays consume significant memory, but can significantly improve runtime speed when 
        # fetching the wave data for a location and time.  
        self._hs = data.variables["swh"][:].data # numpy array of dtype=float32 of shape (len_time, len_latitude, len_longitude)
        self._dp = data.variables["mwd"][:].data # numpy array of dtype=float32 of shape (len_time, len_latitude, len_longitude)

    cdef tuple __get_wave_data_at(self, double longitude, double latitude, int time_diff_from_epoch_in_hrs):
        # Check if latitude, longitude and time are within range of the data in the netcdf file.
        if not ((self._first_latitude <= latitude <= self._last_latitude) or (self._first_latitude >= latitude >= self._last_latitude)):
            raise ValueError("Latitude {} is out of range [{}, {}]".format(latitude, self._first_latitude, self._last_latitude))
        if not ((self._first_longitude <= longitude <= self._last_longitude) or (self._first_longitude >= longitude >= self._last_longitude)):
            raise ValueError("Longitude {} is out of range [{}, {}]".format(longitude, self._first_longitude, self._last_longitude))
        cdef datetime time 
        cdef datetime time_0 
        cdef datetime time_n
        if not ((self._first_time <= time_diff_from_epoch_in_hrs <= self._last_time) or (self._first_time >= time_diff_from_epoch_in_hrs >= self._last_time)):
            time   = self._epoch_time + timedelta(hours=int(time_diff_from_epoch_in_hrs))
            time_0 = self._epoch_time + timedelta(hours=int(self.times[0]))
            time_n = self._epoch_time + timedelta(hours=int(self.times[-1]))
            raise ValueError("Time {} is out of range [{}, {}]".format(time, time_0, time_n))
        # Index 
        cdef int index_time = int((time_diff_from_epoch_in_hrs - self._first_time) // self._delta_time)
        cdef int index_latitude = round((latitude - self._first_latitude) / self._delta_latitude)
        cdef int index_longitude = round((longitude - self._first_longitude) / self._delta_longitude)
        cdef double hs = self._hs[index_time, index_latitude, index_longitude]
        cdef double dp = self._dp[index_time, index_latitude, index_longitude]
        return (hs,dp)

    def get_wave_data_at(self, longitude, latitude, time) -> tuple:
        '''
        Returns the significant wave height and mean wave direction for a given location and time.
        '''
        # Get diff from epoch time in hrs
        time_diff_from_epoch = time - self._epoch_time # time difference w.r.t epoch time
        time_diff_from_epoch_in_sec = time_diff_from_epoch.total_seconds()
        time_diff_from_epoch_in_hrs = divmod(time_diff_from_epoch_in_sec, 3600)[0]
        # NetCDF files have longitude in the range [-180.0, 179.5]. Therefore if the given longitude > 179.5, then round to nearest 0.5
        if longitude >= 179.75 and longitude <= 180.0:
            longitude = -180.0
        elif longitude > 179.5 and longitude < 179.75:
            longitude = 179.5
        wave_hs, wave_dp = self.__get_wave_data_at(longitude, latitude, time_diff_from_epoch_in_hrs)
        # If the location was close to land then wave_hs could be a negative value; correct it by shifting the position slightly to the west. 
        if wave_hs <= 0.0:
            longitude_altered = max(longitude - 0.5, -180.0)
            wave_hs, wave_dp = self.__get_wave_data_at(longitude_altered, latitude, time_diff_from_epoch_in_hrs)
        # Check if the wave data is still -ve
        if wave_hs <= 0.0:
            longitude_altered = min(longitude + 0.5, 180.0)
            wave_hs, wave_dp = self.__get_wave_data_at(longitude_altered, latitude, time_diff_from_epoch_in_hrs)
        # Check if the wave data is still -ve
        if wave_hs <= 0.0:
            latitude_altered = max(latitude - 0.5, -90.0)
            wave_hs, wave_dp = self.__get_wave_data_at(longitude, latitude_altered, time_diff_from_epoch_in_hrs)
        # Check if the wave data is still -ve
        if wave_hs <= 0.0:
            latitude_altered = min(latitude + 0.5, 90.0)
            wave_hs, wave_dp = self.__get_wave_data_at(longitude, latitude_altered, time_diff_from_epoch_in_hrs)
        if wave_hs <= 0.0:
            return (None, None) 
        else:
            wave_dp = wave_dp * math.pi/180 
            return (wave_hs, wave_dp) # m, radians


cdef class NetCDF_wind:
    '''
    Class to extract wind data from an .nc file.
    '''
    def __cinit__(self, file_path) -> None:
        data = nc.Dataset(file_path)
        self._epoch_time = datetime(1900, 1, 1)
        # Longitude range
        self._first_longitude = data.variables["longitude"][0].item()
        self._last_longitude = data.variables["longitude"][-1].item()
        self._delta_longitude = data.variables["longitude"][1].item() - data.variables["longitude"][0].item()
        # Latitude range
        self._first_latitude = data.variables["latitude"][0].item()
        self._last_latitude = data.variables["latitude"][-1].item()
        self._delta_latitude = data.variables["latitude"][1].item() - data.variables["latitude"][0].item()
        # Time range
        self._first_time = data.variables["time"][0].item()
        self._last_time = data.variables["time"][-1].item()
        self._delta_time = data.variables["time"][1].item() - data.variables["time"][0].item()
        # Numpy array of eastward and northward wind velocities.
        self._u = data.variables["u10"][:].data # numpy array of dtype=float32 of shape (len_time, len_latitude, len_longitude)
        self._v = data.variables["v10"][:].data # numpy array of dtype=float32 of shape (len_time, len_latitude, len_longitude)

    cdef tuple __get_wind_speed_at(self, double longitude, double latitude, int time_diff_from_epoch_in_hrs):
        # Check if latitude, longitude and time are within range of the data in the netcdf file.
        if not ((self._first_latitude <= latitude <= self._last_latitude) or (self._first_latitude >= latitude >= self._last_latitude)):
            raise ValueError("Latitude {} is out of range [{}, {}]".format(latitude, self._first_latitude, self._last_latitude))
        if not ((self._first_longitude <= longitude <= self._last_longitude) or (self._first_longitude >= longitude >= self._last_longitude)):
            raise ValueError("Longitude {} is out of range [{}, {}]".format(longitude, self._first_longitude, self._last_longitude))
        cdef datetime time
        cdef datetime time_0
        cdef datetime time_n
        if not ((self._first_time <= time_diff_from_epoch_in_hrs <= self._last_time) or (self._first_time >= time_diff_from_epoch_in_hrs >= self._last_time)):
            time   = self._epoch_time + timedelta(hours=int(time_diff_from_epoch_in_hrs))
            time_0 = self._epoch_time + timedelta(hours=int(self.times[0]))
            time_n = self._epoch_time + timedelta(hours=int(self.times[-1]))
            raise ValueError("Time {} is out of range [{}, {}]".format(time, time_0, time_n))
        # Index
        cdef int index_time = int((time_diff_from_epoch_in_hrs - self._first_time) // self._delta_time)
        cdef int index_latitude = round((latitude - self._first_latitude) / self._delta_latitude)
        cdef int index_longitude = round((longitude - self._first_longitude) / self._delta_longitude)
        cdef double u = self._u[index_time, index_latitude, index_longitude]
        cdef double v = self._v[index_time, index_latitude, index_longitude]
        return (u, v)

    def get_wind_speed_at(self, longitude, latitude, time) -> tuple:
        '''
        Returns the eastward and northward component of wind velocities for a given location and time.
        '''
        # Get diff from epoch time in hrs
        time_diff_from_epoch = time - self._epoch_time # time difference w.r.t epoch time
        time_diff_from_epoch_in_sec = time_diff_from_epoch.total_seconds()
        time_diff_from_epoch_in_hrs = divmod(time_diff_from_epoch_in_sec, 3600)[0]
        # NetCDF files have longitude in the range [-180.0, 179.5]. Therefore if the given longitude > 179.5, then round to nearest 0.5
        if longitude >= 179.75 and longitude <= 180.0:
            longitude = -180.0
        elif longitude > 179.5 and longitude < 179.75:
            longitude = 179.5
        u, v = self.__get_wind_data_at(longitude, latitude, time_diff_from_epoch_in_hrs)
        return (u, v) # m/s


cdef class NetCDF_precipitation:
    '''
    Class to extract precipitation data from an .nc file.
    '''
    def __cinit__(self, file_path) -> None:
        data = nc.Dataset(file_path)
        self._epoch_time = datetime(1900, 1, 1)
        # Longitude range
        self._first_longitude = data.variables["longitude"][0].item()
        self._last_longitude = data.variables["longitude"][-1].item()
        self._delta_longitude = data.variables["longitude"][1].item() - data.variables["longitude"][0].item()
        # Latitude range
        self._first_latitude = data.variables["latitude"][0].item()
        self._last_latitude = data.variables["latitude"][-1].item()
        self._delta_latitude = data.variables["latitude"][1].item() - data.variables["latitude"][0].item()
        # Time range
        self._first_time = data.variables["time"][0].item()
        self._last_time = data.variables["time"][-1].item()
        self._delta_time = data.variables["time"][1].item() - data.variables["time"][0].item()
        # Numpy array of total precipitation.
        self._tp = data.variables["tp"][:].data # numpy array of dtype=float32 of shape (len_time, len_latitude, len_longitude)

    cdef double __get_precipitation_at(self, double longitude, double latitude, int time_diff_from_epoch_in_hrs):
        # Check if latitude, longitude and time are within range of the data in the netcdf file.
        if not ((self._first_latitude <= latitude <= self._last_latitude) or (self._first_latitude >= latitude >= self._last_latitude)):
            raise ValueError("Latitude {} is out of range [{}, {}]".format(latitude, self._first_latitude, self._last_latitude))
        if not ((self._first_longitude <= longitude <= self._last_longitude) or (self._first_longitude >= longitude >= self._last_longitude)):
            raise ValueError("Longitude {} is out of range [{}, {}]".format(longitude, self._first_longitude, self._last_longitude))
        cdef datetime time
        cdef datetime time_0
        cdef datetime time_n
        if not ((self._first_time <= time_diff_from_epoch_in_hrs <= self._last_time) or (self._first_time >= time_diff_from_epoch_in_hrs >= self._last_time)):
            time   = self._epoch_time + timedelta(hours=int(time_diff_from_epoch_in_hrs))
            time_0 = self._epoch_time + timedelta(hours=int(self.times[0]))
            time_n = self._epoch_time + timedelta(hours=int(self.times[-1]))
            raise ValueError("Time {} is out of range [{}, {}]".format(time, time_0, time_n))
        # Index
        cdef int index_time = int((time_diff_from_epoch_in_hrs - self._first_time) // self._delta_time)
        cdef int index_latitude = round((latitude - self._first_latitude) / self._delta_latitude)
        cdef int index_longitude = round((longitude - self._first_longitude) / self._delta_longitude)
        cdef double tp = self._tp[index_time, index_latitude, index_longitude]
        return tp

    def get_precipitation_at(self, longitude, latitude, time) -> tuple:
        '''
        Returns the total precipitation for a given location and time.
        '''
        # Get diff from epoch time in hrs
        time_diff_from_epoch = time - self._epoch_time # time difference w.r.t epoch time
        time_diff_from_epoch_in_sec = time_diff_from_epoch.total_seconds()
        time_diff_from_epoch_in_hrs = divmod(time_diff_from_epoch_in_sec, 3600)[0]
        # NetCDF files have longitude in the range [-180.0, 179.5]. Therefore if the given longitude > 179.5, then round to nearest 0.5
        if longitude >= 179.75 and longitude <= 180.0:
            longitude = -180.0
        elif longitude > 179.5 and longitude < 179.75:
            longitude = 179.5
        tp = self.__get_precipitation_at(longitude, latitude, time_diff_from_epoch_in_hrs)
        return tp # m


cdef class NetCDF_current:
    '''
    Class to extract zonal and meridional velocities from a .nc file.
    '''
    def __cinit__(self, file_path) -> None:
        data = nc.Dataset(file_path)
        self._epoch_time = datetime(1950, 1, 1)  
        # Longitude range
        self._first_longitude = data.variables["longitude"][0].item()
        self._last_longitude = data.variables["longitude"][-1].item()
        self._delta_longitude = data.variables["longitude"][1].item() - data.variables["longitude"][0].item()
        # Latitude range
        self._first_latitude = data.variables["latitude"][0].item()
        self._last_latitude = data.variables["latitude"][-1].item()
        self._delta_latitude = data.variables["latitude"][1].item() - data.variables["latitude"][0].item()
        # Time range
        self._time = data.variables["time"][0].item()
        # Ocean current velocities
        self._zonal_velocities = data.variables["ugos"][:].data # numpy array of dtype=float32 of shape (len_time, len_latitude, len_longitude)
        self._meridional_velocities = data.variables["vgos"][:].data # numpy array of dtype=float32 of shape (len_time, len_latitude, len_longitude)
        
    cdef tuple __get_ocean_current_at(self, double longitude, double latitude, int time_diff_from_epoch_in_days):
        # # Check if latitude, longitude and time are within range of the data in the netcdf file.
        # if not ((self._first_latitude <= latitude <= self._last_latitude) or (self._first_latitude >= latitude >= self._last_latitude)):
        #     raise ValueError("Latitude {} is out of range [{}, {}]".format(latitude, self._first_latitude, self._last_latitude))
        # if not ((self._first_longitude <= longitude <= self._last_longitude) or (self._first_longitude >= longitude >= self._last_longitude)):
        #     raise ValueError("Longitude {} is out of range [{}, {}]".format(longitude, self._first_longitude, self._last_longitude))
        cdef datetime time 
        if (self._time != time_diff_from_epoch_in_days):
            time  = self._epoch_time + timedelta(days=int(time_diff_from_epoch_in_days))
            time2 = self._epoch_time + timedelta(days=int(self._time))
            raise ValueError("Time {} out of range. Time on file {}.".format(time, time2))
        # Index 
        cdef int index_time = 0
        cdef int index_latitude = round((latitude - self._first_latitude) / self._delta_latitude)
        cdef int index_longitude = round((longitude - self._first_longitude) / self._delta_longitude)
        cdef double zonal_velocity = self._zonal_velocities[index_time, index_latitude, index_longitude]
        cdef double meridional_velocity = self._meridional_velocities[index_time, index_latitude, index_longitude]
        return (zonal_velocity, meridional_velocity)

    def get_ocean_current_at(self, longitude, latitude, time) -> tuple:
        '''
        Returns the zonal and meridional velocities for a given location and depth.
        '''
        # NetCDF files have longitude in the range [-180.0, 179.5]. Therefore if the given longitude > 179.5, then round to nearest 0.5
        # Get diff from epoch time in hrs
        time_diff_from_epoch = time - self._epoch_time # time difference w.r.t epoch time
        time_diff_from_epoch_in_sec = time_diff_from_epoch.total_seconds()
        time_diff_from_epoch_in_days = divmod(time_diff_from_epoch_in_sec, 86400)[0]
        if longitude > 179.875:
            longitude = 179.875
        elif longitude < -179.875:
            longitude = -179.875
        zonal, meridional = self.__get_ocean_current_at(longitude, latitude, time_diff_from_epoch_in_days)

        # If the location was close to land then wave_hs could be a negative value; correct it by shifting the position slightly to the west. 
        if abs(zonal) >= 100.0 or abs(meridional) >= 100:
            longitude_altered = max(longitude - 0.5, -180.0)
            zonal, meridional = self.__get_ocean_current_at(longitude_altered, latitude, time_diff_from_epoch_in_days)
        # Check if the wave data is still -ve
        if abs(zonal) >= 100.0 or abs(meridional) >= 100:
            longitude_altered = min(longitude + 0.5, 180.0)
            zonal, meridional = self.__get_ocean_current_at(longitude_altered, latitude, time_diff_from_epoch_in_days)
        # Check if the wave data is still -ve
        if abs(zonal) >= 100.0 or abs(meridional) >= 100:
            latitude_altered = max(latitude - 0.5, -90.0)
            zonal, meridional = self.__get_ocean_current_at(longitude, latitude_altered, time_diff_from_epoch_in_days)
        # Check if the wave data is still -ve
        if abs(zonal) >= 100.0 or abs(meridional) >= 100:
            latitude_altered = min(latitude + 0.5, 90.0)
            zonal, meridional = self.__get_ocean_current_at(longitude, latitude_altered, time_diff_from_epoch_in_days)
        if abs(zonal) >= 100.0 or abs(meridional) >= 100:
            return (None, None) 
        else:
            return (zonal, meridional) # m/s, m/s