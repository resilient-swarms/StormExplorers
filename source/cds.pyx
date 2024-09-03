import cdsapi
import pathlib
import tarfile
import os

cpdef void _get_reanalysis_era5_single_levels_data(int year, int month, double north, double south, double east, double west, list variables, str file_path):
    '''
    Retrieve data from Climate Data Store (CDS): https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview
    '''
    client = cdsapi.Client()
    client.retrieve('reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'variable': variables,
            'year': year,
            'month': month,
            'day': list(range(32)),
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': [north, west, south, east],
            "format": "netcdf"
        }, file_path)


cpdef void get_wave_data(int year, int month, double north, double south, double east, double west, str file_path):
    '''
    Retrieve wave data from Climate Data Store (CDS): https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview
    '''
    cdef list variables = ['mean_wave_direction', 'significant_height_of_combined_wind_waves_and_swell']
    _get_reanalysis_era5_single_levels_data(year, month, north, south, east, west, variables, file_path)


cpdef void get_precipitation_data(int year, int month, double north, double south, double east, double west, str file_path):
    '''
    Retrieve total precipitation from Climate Data Store (CDS): https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview
    '''
    cdef list variables = ['total_precipitation',]
    _get_reanalysis_era5_single_levels_data(year, month, north, south, east, west, variables, file_path)


cpdef void get_wind_data(int year, int month, double north, double south, double east, double west, str file_path):
    '''
    Retrieve wind speed data from Climate Data Store (CDS): https://cds-beta.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview
    '''
    cdef list variables = ['10m_u_component_of_wind', '10m_v_component_of_wind',]
    _get_reanalysis_era5_single_levels_data(year, month, north, south, east, west, variables, file_path)

cpdef void get_ocean_current_data(int year, int month, str file_path):
    '''
    Retrieve ocean current data from Climate Data Store (CDS): https://cds-beta.climate.copernicus.eu/datasets/satellite-sea-level-global?tab=overview
    '''
    client = cdsapi.Client()
    days = []
    for day in range(1,32):
        days.append(str(day).zfill(2))
    client.retrieve('satellite-sea-level-global',
        {
            'year': year,
            'month': str(month).zfill(2),
            'day': days,
            'version': 'vDT2021',
            'variable': 'daily',
            'format': 'tgz',
        }, file_path)
    # Open the downloaded tar file
    tar_file = tarfile.open(file_path)
    # Extract files
    archieve_dir = pathlib.Path(file_path).parent.resolve()
    tar_file.extractall(archieve_dir)
    tar_file.close
    # Rename extracted files    dt_global_twosat_phy_l4_20190303_vDT2021
    for day in range(32):
        if os.path.exists(archieve_dir.joinpath("dt_global_twosat_phy_l4_{}{}{}_vDT2021.nc".format(year, str(month).zfill(2), str(day).zfill(2)))):
            os.rename(archieve_dir.joinpath("dt_global_twosat_phy_l4_{}{}{}_vDT2021.nc".format(year, str(month).zfill(2), str(day).zfill(2))), 
                    archieve_dir.joinpath(file_path[:-7]+"_{}.nc".format(str(day).zfill(2))))