cimport numpy as np
from cpython.datetime cimport datetime


cdef class NetCDF_wave:
    cdef datetime _epoch_time
    # Longitude range
    cdef double _first_longitude
    cdef double _last_longitude
    cdef double _delta_longitude
    # Latitude range
    cdef double _first_latitude 
    cdef double _last_latitude 
    cdef double _delta_latitude 
    # Time range
    cdef double _first_time 
    cdef double _last_time 
    cdef double _delta_time 
    # Numpy array of significant wave hts and wave headings. 
    cdef np.ndarray _hs 
    cdef np.ndarray _dp 
    # Methods
    cdef tuple __get_wave_data_at(self, double longitude, double latitude, int time)


cdef class NetCDF_wind:
    cdef datetime _epoch_time
    # Longitude range
    cdef double _first_longitude
    cdef double _last_longitude
    cdef double _delta_longitude
    # Latitude range
    cdef double _first_latitude
    cdef double _last_latitude
    cdef double _delta_latitude
    # Time range
    cdef double _first_time
    cdef double _last_time
    cdef double _delta_time
    # Numpy array of wind velocities.
    cdef np.ndarray _u # Eastward component of wind velocity (m/s)
    cdef np.ndarray _v # Northward component of wind velocity (m/s)
    # Methods
    cdef tuple __get_wind_speed_at(self, double longitude, double latitude, int time)


cdef class NetCDF_precipitation:
    cdef datetime _epoch_time
    # Longitude range
    cdef double _first_longitude
    cdef double _last_longitude
    cdef double _delta_longitude
    # Latitude range
    cdef double _first_latitude
    cdef double _last_latitude
    cdef double _delta_latitude
    # Time range
    cdef double _first_time
    cdef double _last_time
    cdef double _delta_time
    # Numpy array for total precipitation.
    cdef np.ndarray _tp # Total precipitation (m)
    # Methods
    cdef double __get_precipitation_at(self, double longitude, double latitude, int time)


cdef class NetCDF_current:
    cdef datetime _epoch_time
    # Longitude range
    cdef double _first_longitude
    cdef double _last_longitude
    cdef double _delta_longitude
    # Latitude range
    cdef double _first_latitude 
    cdef double _last_latitude 
    cdef double _delta_latitude 
    # Time range
    cdef double _time 
    # Numpy array of zonal and meridional velocities. 
    cdef np.ndarray _zonal_velocities 
    cdef np.ndarray _meridional_velocities
    # Methods
    cdef tuple __get_ocean_current_at(self, double longitude, double latitude, int time)