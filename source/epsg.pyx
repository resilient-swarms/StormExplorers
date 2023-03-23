from pyproj import Proj, Transformer

# Methods for tranforming GCS to PCS and vice versa. 

# Reference
# ---------
# [1] https://epsg.io/
# [2] https://epsg.io/4326
# [3] https://epsg.io/3857

transformer_epsg_GCS_to_PCS = Transformer.from_crs("epsg:4326", "epsg:3857", always_xy = True)
transformer_epsg_PCS_to_GCS = Transformer.from_crs("epsg:3857", "epsg:4326", always_xy = True)

cpdef tuple GCS_to_PCS(double longitude, double latitude):
    '''
    Function to convert Geographic Coordinate System (GCS) to Projected Coordinate System (PCS).
    Returns (longitude, latitude) in meters.
    '''
    cdef double x 
    cdef double y
    x, y = transformer_epsg_GCS_to_PCS.transform(longitude, latitude, errcheck = True) # x is longitude in meter and y is latitude in meter. 
    return (x, y)

cpdef tuple PCS_to_GCS(double x, double y):
    '''
    Function to convert Geographic Coordinate System (GCS) to Projected Coordinate System (PCS).
    x is the longitudes in meters and y is latitude in meters. 
    Returns (longitude, latitude).
    '''
    cdef double longitude
    cdef double latititude
    longitude, latitude = transformer_epsg_PCS_to_GCS.transform(x, y, errcheck = True)  
    return (longitude, latitude)