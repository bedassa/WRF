
import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from wrf import getvar, to_np, get_cartopy, latlon_coords
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset, num2date
import pandas as pd

def extract_avg_t2m_in_box(wrf_file, T2, XLAT, XLONG, lat1, lat2, lon1, lon2, Times, start_time, end_time):
    """
    Extract time series of average 2-meter temperature over a lat-lon box.
    """
    import xarray as xr
    import numpy as np
    import pandas as pd

    # Open dataset
    ds = xr.open_dataset(wrf_file, decode_times=False, engine="netcdf4")

    # Standardize time variable
    if "XTIME" in ds:
        ds = ds.rename({"XTIME": "Times"})
        times = pd.to_datetime(ds[Times].values, unit='m', origin=pd.Timestamp('2022-03-01 00:00:00'))
        ds[Times] = times.round('T')
    elif "valid_time" in ds:
        ds = ds.rename({"valid_time": "Times"})
        times = pd.to_datetime(ds['Times'].values, unit='s', origin='unix')
        ds[Times] = times.round('T')

    # Convert T2 to Celsius
    T2m_full = ds[T2] - 273.15
    T2m_full = T2m_full.assign_coords(Times=ds[Times])

    # Select time range
    T2m_range = T2m_full.sel(Times=slice(start_time, end_time))

    # Get lat/lon
    lats = ds[XLAT].squeeze()
    lons = ds[XLONG].squeeze()

    # Build mask for box selection
    mask = (lats >= lat1) & (lats <= lat2) & (lons >= lon1) & (lons <= lon2)

    # Apply mask and average over lat/lon dimensions
    T2m_box = T2m_range.where(mask)
    T2m_avg = T2m_box.mean(dim=[dim for dim in T2m_box.dims if dim != "Times"], skipna=True)

    # Build time axis
    time_axis = T2m_avg["Times"].values

    return {
        'lat_range': (lat1, lat2),
        'lon_range': (lon1, lon2),
        'T2m_avg_C': T2m_avg.values,
        'time_axis': time_axis
    }

wrf_d01 = r'C:\Bedassa\SoilTemp_Processed\WRF\wrfout_d01_2022-03-01_00_00_00_T2SMOISTSLB.NC'
wrf_d02 = r'C:\Bedassa\SoilTemp_Processed\WRF\wrfout_d02_2022-03-01_00_00_00_T2SMOISTSLB.NC'
wrf_d03 = r'C:\Bedassa\SoilTemp_Processed\WRF\wrfout_d03_2022-03-01_00_00_00_T2SMOISTSLB.NC'
wrf_d04 = r'C:\Bedassa\SoilTemp_Processed\WRF\wrfout_d04_2022-03-01_00_00_00_T2SMOISTSLB.NC'

ERA5_9km = r'C:\\Bedassa\\SoilTemp_Processed\\WRF\\ERA5_2022_03_01_05.nc'
ERA5_31km = r'C:\\Bedassa\\SoilTemp_Processed\\WRF\\ERA5_31km_2022_03_01_05.nc'
Times ='Times'
XLAT = 'XLAT'
XLONG = 'XLONG'
T2 = 'T2'
# (T2m_avg.lat >= 59.36366272) & (T2m_avg.lat <= 61.50312805) &
#     (T2m_avg.lon >= 4.61058521 ) & (T2m_avg.lon <= 11.48741436),
   # D01
lat1 = 59.36366272
lat2 = 61.50312805
lon1 =4.61058521 
lon2 = 11.48741436
#d02 
# subset = T2m_avg.where(
#     (T2m_avg.lat >=60.14520264) & (T2m_avg.lat <= 60.90861511) &
#     (T2m_avg.lon >= 6.90286112 ) & (T2m_avg.lon <= 9.72412395),
#     drop=True
# )
lat1 = 60.14520264
lat2 = 60.90861511
lon1 = 6.90286112
lon2 = 9.72412395

# d04
#D04
# subset = T2m_avg.where(
#     (T2m_avg.lat >=60.38815308) & (T2m_avg.lat <= 60.68687439) &
#     (T2m_avg.lon >= 7.73356724 ) & (T2m_avg.lon <= 8.67006969),
#     drop=True
# )
lat1 = 60.38815308
lat2 = 60.68687439
lon1 = 7.73356724
lon2 = 8.67006969


start_time = '2022-03-02 00:00:00'
end_time = '2022-03-03 21:00:00'

T_1 = extract_avg_t2m_in_box(wrf_d01, T2, XLAT, XLONG, lat1, lat2, lon1, lon2, Times, start_time, end_time)
T_2 = extract_avg_t2m_in_box(wrf_d02, T2, XLAT, XLONG, lat1, lat2, lon1, lon2, Times, start_time, end_time)
T_3 = extract_avg_t2m_in_box(wrf_d03, T2, XLAT, XLONG, lat1, lat2, lon1, lon2, Times, start_time, end_time)
T_4 = extract_avg_t2m_in_box(wrf_d04, T2, XLAT, XLONG, lat1, lat2, lon1, lon2, Times, start_time, end_time)

def extract_avg_t2m_in_box_ERA5(wrf_file, T2, XLAT, XLONG, lat1, lat2, lon1, lon2, start_time, end_time):
    """
    Extract time series of average 2-meter temperature over a lat-lon box from ERA5 data.
    """
    import xarray as xr
    import numpy as np
    import pandas as pd

    # Open dataset
    ds = xr.open_dataset(wrf_file, decode_times=False, engine="netcdf4")

    # Standardize time variable
    if "XTIME" in ds:
        ds = ds.rename({"XTIME": "Times"})
        times = pd.to_datetime(ds["Times"].values, unit='m', origin=pd.Timestamp('2022-03-01 00:00:00'))
        ds["Times"] = times.round('T')
    elif "valid_time" in ds:
        ds = ds.rename({"valid_time": "Times"})
        times = pd.to_datetime(ds["Times"].values, unit='s', origin='unix')
        ds["Times"] = times.round('T')

    # Convert T2 to Celsius
    T2m_full = ds[T2] - 273.15
    T2m_full = T2m_full.assign_coords(Times=ds["Times"])

    # Select time range
    T2m_range = T2m_full.sel(Times=slice(start_time, end_time))

    # Extract lat/lon
    lats = ds[XLAT].values
    lons = ds[XLONG].values

    # If lat/lon are 1D, expand to 2D mesh
    if lats.ndim == 1 and lons.ndim == 1:
        lat_grid, lon_grid = np.meshgrid(lats, lons, indexing='ij')
    else:
        lat_grid, lon_grid = lats, lons

    # Create mask for the bounding box
    mask = (lat_grid >= lat1) & (lat_grid <= lat2) & (lon_grid >= lon1) & (lon_grid <= lon2)

    # Mask temperature data and average over spatial dimensions
    masked_T2m = T2m_range.where(mask)
    T2m_avg_time_series = masked_T2m.mean(dim=[dim for dim in masked_T2m.dims if dim != "Times"], skipna=True)

    # Extract time axis
    time_axis = T2m_avg_time_series["Times"].values

    return {
        'lat_range': (lat1, lat2),
        'lon_range': (lon1, lon2),
        'T2m_avg_C': T2m_avg_time_series.values,
        'time_axis': time_axis
    }


T2 = 't2m'
XLAT = 'latitude'
XLONG = 'longitude'

ERA5_9 = extract_avg_t2m_in_box_ERA5(ERA5_9km,T2,XLAT,XLONG, lat1, lat2, lon1, lon2, start_time, end_time)


ERA5_31 = extract_avg_t2m_in_box_ERA5(ERA5_31km,T2,XLAT,XLONG, lat1, lat2, lon1, lon2, start_time, end_time)



# Define expected time axis
time_axis = pd.date_range(start='2022-03-02 00:00:00', end='2022-03-03 21:00:00', freq='H')

# Trim each temperature series to match time axis length
t2_d01 = T_1['T2m_avg_C'][:len(time_axis)]
t2_d02 = T_2['T2m_avg_C'][:len(time_axis)]
t2_d03 = T_3['T2m_avg_C'][:len(time_axis)]
t2_d04 = T_4['T2m_avg_C'][:len(time_axis)]
t2_era5 = ERA5_9['T2m_avg_C'][:len(time_axis)]
t2_era5_31 = ERA5_31['T2m_avg_C'][:len(time_axis)]


# Plotting


plt.figure(figsize=(14, 6))
plt.plot(time_axis, t2_era5, label='ERA5 9km', linestyle='--', color='blue')
plt.plot(time_axis, t2_era5_31, label='ERA5 31km', linestyle='--', color='orange') 
plt.plot(time_axis, t2_d01, label='WRF 9km')
plt.plot(time_axis, t2_d02, label='WRF 3km')
plt.plot(time_axis, t2_d03, label='WRF 1km')
plt.plot(time_axis, t2_d04, label='WRF 330m')
# Labels and legend
plt.title("2m temperature (°C) time series mean over domain1")
plt.xlabel("Time")
plt.ylabel("Temperature (°C)")
plt.legend(title="", loc='best', ncol=6)
plt.grid(True)
plt.tight_layout()
plt.show()



