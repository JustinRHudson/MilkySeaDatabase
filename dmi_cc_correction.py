'''
    This scripts detrends the HadISST SST data so it can be used to produce a
    climate change corrected Dipole Mode Index for the Indian Oceal Dipole.
    The resulting data is saved to ./DATA/ as  CC_CORRECTED_DMI.txt

    This script assumes that HadISST_sst.nc has already been downloaded from:
    https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
'''

###### IMPORTS GO HERE ######
import numpy as np
from netCDF4 import Dataset
import datetime as dt
import numpy.ma as ma
import xarray as xr

###### GLOBALS GO HERE ######
IOD_WEST_BOX = [50,-10,70,10] #left edge, bottom edge, right edge, top edge
IOD_EAST_BOX = [90,-10,110,0]

###### FUNCTIONS GO HERE ######
def load_HadISST_data() -> tuple[np.ndarray]:
    '''
        Loads the HadISST data into numpy arrays.

        Inputs:
            None
        
        Returns:
            h_dates (np.ndarray): The dates coresponding to the HadISST SSTs
            h_lon (np.ndarray): The longitudes corresponding to the HadISST SSTs
            h_lat (np.ndarray): The latitudes corresponding to the HadISST SSTs
            h_sst (np.ndarray): The HadISST SSTs
    '''
    H_FILE = Dataset('./DATA/HadISST_sst.nc')
    h_time = H_FILE.variables['time'][:] #days since 1870,1,1
    h_dates = np.array([dt.datetime(1870,1,1)+dt.timedelta(days = int(t)) for t in h_time])
    h_lat = H_FILE.variables['latitude'][:]
    h_lon = H_FILE.variables['longitude'][:]
    h_sst = ma.getdata(H_FILE.variables['sst'][:]) #degrees C, fill -1e+30 (time,lat,lon)
    h_sst[np.where(h_sst <= -999)] = np.nan

    return h_dates,h_lat,h_lon,h_sst

# Let's detrend at each gridpoint
# Step 1: Pick a grid point and check it isn't land
# Step 2: Fit a line to the time series there
# Step 3: Apply correction so the slope is now zero
# Step 4: Go to next grid point and repeat
def check_if_grid_point_is_ocean(grid_point_ts:np.ndarray) -> bool:
    ''' Checks if a grid point is ocean based on the time series. Because
    I am looking at an SST dataset I can just see the ratio of NaNs to non-NaNs

    Inputs:
    grid_point_ts (np.ndarray): The Time Series of the grid point

    Returns:
    is_ocean (bool): A Boolean that is true if the grid cell is ocean
    '''

    is_ocean = True
    nan_list = np.isnan(grid_point_ts)

    if np.sum(nan_list) >=  1:
        is_ocean = False
    
    return is_ocean
        

def line_fit_to_grid_point(grid_point_ts:np.ndarray) -> tuple[float,float]:
    ''' Fits a line to a grid point time series and returns the slope and intercept
    of the line fit to the time series.
    
    Inputs: 
    grid_point_ts (np.ndarray): The Time series at the grid point
    
    Outputs: 
    slope (float): The slope of the line fit to the time series
    y-intercept (float): The y-intercept of the line fit to the time series
    '''

    fake_xs = np.arange(0,len(grid_point_ts),1)

    line_fit = np.polyfit(fake_xs,grid_point_ts,1)

    slope = line_fit[0]
    y_intercept = line_fit[1]

    return slope,y_intercept


def slope_correction(grid_point_ts:np.ndarray,slope:float,y_intercept:float) -> np.ndarray:
    ''' Corrects the slope of the time series so that the climate change signal has been removed

    Inputs:
    grid_point_ts (np.ndarray); The Time Series at the grid point
    slope (float): The slope of the line fit to grid_point_ts
    y_intercept(float): The y-intercept of the line fit to grid_point_ts

    Outputs:
    corrected_ts (np.ndarray): The slope corrected time series
    '''

    correct_xs = np.arange(0,len(grid_point_ts),1)

    correction_line = slope*correct_xs+y_intercept

    fit_dif = (grid_point_ts - correction_line)

    corrected_ts = fit_dif + y_intercept

    return corrected_ts

def correction_loop(sst_dataset:np.ndarray) -> np.ndarray:
    ''' Loops over the SST dataset and applies the time series correction at each
    grid point to remove the climate change signal
    

    Input:
    sst_dataset (np.ndarray): The sst dataset to have the climate change correction applied

    Output:
    corrected_dataset (np.ndarray): The climate change corrected dataset
    '''

    corrected_dataset = np.empty(sst_dataset.shape)
    corrected_dataset[:,:,:] = np.nan

    for i in range(sst_dataset.shape[1]):
        for j in range(sst_dataset.shape[2]):
            if check_if_grid_point_is_ocean(sst_dataset[:,i,j]):
                point_slope,point_int = line_fit_to_grid_point(sst_dataset[:,i,j])
                corrected_ts = slope_correction(sst_dataset[:,i,j],point_slope,point_int)
                corrected_dataset[:,i,j] = np.copy(corrected_ts)
    
    return corrected_dataset

def date_grouper(sst_dates:np.ndarray) -> np.ndarray:
    '''Gets the indices of dates for each calendar day and returns them
    so that the sst dataset can be indexed on those dates to get monthly means
    across the entire monthly sst dataset.
    
    Inputs:
    sst_dates (np.ndarray): The dates corresponding to the time axis of the sst dataset.
    
    Outputs:
    date_inds (np.ndarray): A (12) length numpy array where each element is a numpy array
    of varying length with the indices for month N of 12.'''

    jan_ray = []
    feb_ray = []
    mar_ray = []
    apr_ray = []
    may_ray = []
    jun_ray = []
    jul_ray = []
    aug_ray = []
    sep_ray = []
    oct_ray = []
    nov_ray = []
    dec_ray = []

    for i in range(len(sst_dates)):
        if sst_dates[i].month == 1:
            jan_ray.append(i)
        elif sst_dates[i].month == 2:
            feb_ray.append(i)
        elif sst_dates[i].month == 3:
            mar_ray.append(i)
        elif sst_dates[i].month == 4:
            apr_ray.append(i)
        elif sst_dates[i].month == 5:
            may_ray.append(i)
        elif sst_dates[i].month == 6:
            jun_ray.append(i)
        elif sst_dates[i].month == 7:
            jul_ray.append(i)
        elif sst_dates[i].month == 8:
            aug_ray.append(i)
        elif sst_dates[i].month == 9:
            sep_ray.append(i)
        elif sst_dates[i].month == 10:
            oct_ray.append(i)
        elif sst_dates[i].month == 11:
            nov_ray.append(i)
        elif sst_dates[i].month == 12:
            dec_ray.append(i)

    date_inds = np.array([np.array(jan_ray),np.array(feb_ray),np.array(mar_ray),np.array(apr_ray),
                          np.array(may_ray),np.array(jun_ray),np.array(jul_ray),np.array(aug_ray),
                          np.array(sep_ray),np.array(oct_ray),np.array(nov_ray),np.array(dec_ray)])
    
    return date_inds

def get_monthly_means(monthly_date_inds:np.ndarray,sst_dataset:np.ndarray) -> np.ndarray:
    '''Gets the mean of each month across the entire SST dataset and returns a
    12 x lat x lon numpy array containing the monthly means so anomalies can be
    calculated easily.
    
    Inputs:
    monthly_date_inds (np.ndarray): The indices for each month
    sst_dataset (np.ndarray): The sst dataset I am trying to get the anomalies for
    
    Outputs:
    monthly_means (np.ndarray): 12 x lat x lon numpy array containing the monthly means
    for the anomalies to be calculated relative to.'''

    monthly_means = np.empty((12,sst_dataset.shape[1],sst_dataset.shape[2]))

    for i in range(len(monthly_date_inds)):
        monthly_means[i,:,:] = np.nanmean(sst_dataset[monthly_date_inds[i]],axis = (0))
    
    return monthly_means

def get_anomalous_ssts(monthly_date_inds:np.ndarray,monthly_means:np.ndarray,sst_dataset:np.ndarray) -> np.ndarray:
    '''Calculates the anomalous SSTs.
    
    Inputs:
    monthly_date_inds (np.ndarray): The indices for each month
    sst_dataset (np.ndarray): The sst dataset I am trying to get the anomalies for
    monthly_means (np.ndarray): 12 x lat x lon numpy array containing the monthly means
    for the anomalies to be calculated relative to.
    
    Outputs:
    anomalous_sst (np.ndarray): The anomalous ssts I am trying to calculate.'''

    anomalous_sst = np.copy(sst_dataset)

    for i in range(len(monthly_date_inds)):
        anomalous_sst[monthly_date_inds[i]] = sst_dataset[monthly_date_inds[i]] - monthly_means[i]
    
    return anomalous_sst

def box_data(box:list,data_lats:np.ndarray,data_lons:np.ndarray,dataset:np.ndarray):
    '''Given a Box and the lats/lons for a dataset return the data within the box
    as time,lat,lon dataset.
    
    Inputs:
    box (list): A box detailing the left, bottom, right, and top edge of the region fo interest
    data_lats (np.ndarray): The latitudes for the dataset
    data_lons (np.ndarray): The longitudes for the dataset
    dataset (np.ndarray): The dataset I want to subset
    
    Outputs:
    box_time_series (np.ndarray): A 1D spatial mean time series of the box
    '''

    lat_inds = np.intersect1d(np.where(data_lats <= box[3]),np.where(data_lats >= box[1]))
    lon_inds = np.intersect1d(np.where(data_lons <= box[2]),np.where(data_lons >= box[0]))

    box_data = np.array(dataset[:,lat_inds[0]:lat_inds[-1]+1,lon_inds[0]:lon_inds[-1]+1])
    box_lats = np.array(data_lats[lat_inds[0]:lat_inds[-1]+1])
    box_lons = np.array(data_lons[lon_inds[0]:lon_inds[-1]+1])

    return box_lats,box_lons,box_data

def great_circle_distance(lat1:float,lon1:float,lat2:float,lon2:float) -> float:
        lat1,lat2 = np.deg2rad(lat1),np.deg2rad(lat2)
        delta_lon = np.deg2rad(np.abs(lon1-lon2))
        delta_sigma = np.arccos(np.sin(lat1)*np.sin(lat2) + np.cos(lat1)*np.cos(lat2)*np.cos(delta_lon))
        d = 6371. * delta_sigma #multiply by Radius of Earth in kilometers

        return d

#get the values of dx and dy between grid points on my center grid
def calc_dx_dy(latitudes:np.ndarray,longitudes:np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    #use great circle distance to get dx and dy on my center grids
    dx_array = np.zeros((len(longitudes),len(latitudes)))
    dy_array = np.zeros((len(longitudes),len(latitudes)))
    dlat = np.abs(latitudes[0]-latitudes[1])
    #dx is constant along a latitude and dy only varies with latitude
    for i in range(len(latitudes)):
        dx_array[:,i] += great_circle_distance(latitudes[i],longitudes[0],latitudes[i],longitudes[1])
        dy_array[:,i] += great_circle_distance(latitudes[i],longitudes[0],latitudes[i]+dlat,longitudes[0])
    
    #swap the axes of the dx and dy arrays
    dx_array = np.swapaxes(dx_array,0,1)
    dy_array = np.swapaxes(dy_array,0,1)
    
    return dx_array,dy_array


def area_average(data_lats:np.ndarray,data_lons:np.ndarray,dataset:np.ndarray) -> np.ndarray:
    '''Given a time,lat,lon dataset calculate the area average and return a 1D
    time series of this dataset.
    
    Inputs:
    data_lats (np.ndarray): The latitudes corresponding to the dataset
    data_lons (np.ndarray): The longitudes corresponding to the dataset
    dataset (np.ndarray): The dataset I am trying to take the area average of
    
    Outputs:
    area_ave_data (np.ndarray): The area averaged time series'''

    # get the DX and DY arrays to calculate area with
    dx_array,dy_array = calc_dx_dy(data_lats,data_lons)

    #assume close enough to rectangles locally that I can just multipy dx and dy
    grid_areas = dx_array*dy_array #units are km**2

    #mark which ones are not nans
    non_nan_inds = np.zeros((dataset.shape[1],dataset.shape[2]))
    for i in range(len(data_lats)):
        for j in range(len(data_lons)):
            if check_if_grid_point_is_ocean(dataset[:,i,j]):
                non_nan_inds[i,j] = 1
    
    #get the total area if the non-nan areas
    non_nan_area = 0
    for i in range(non_nan_inds.shape[0]):
        for j in range(non_nan_inds.shape[1]):
            if non_nan_inds[i,j]:
                non_nan_area += grid_areas[i,j]
    
    #now convert grid area to pct area of the non-nan area
    pct_area = grid_areas / non_nan_area

    #for each time step I need to multiply each non-nan area by their
    #pct area and take the sum
    area_ave_data = np.empty(dataset.shape[0])
    for i in range(dataset.shape[0]):
        area_ave_data[i] = np.nansum(dataset[i,:,:]*pct_area)
    
    return area_ave_data

def write_data(DMI_Timeseries:np.ndarray,DMI_Dates:np.ndarray) -> None:
    '''
        Writes the data to the ./DATA/CC_CORRECTED_DMI.txt file

        Parameters:
            DMI_Timeseries (np.ndarray): The calculated DMI timeseries
            DMI_Dates (np.ndarray): The dates corresponding to the DMI data
        
        Returns:
            None
    '''
    #set the reference date for the data
    ref_date = dt.datetime(1870,1,1)
    #open the file and write to it
    f = open('./DATA/CC_CORRECTED_DMI.txt','w')
    f.write('#This is a climate change corrected Dipole Mode Index (DMI) for the\n#Indian Ocean Dipole (IOD) derived from HADISST 1.1 SSTs.\n')
    f.write('#The first column is days since January 1st, 1870 (1870,1,1)\n#The second column is the index value in degrees Celsius.\n')
    f.write('#DAYS,SST\n')
    for i in range(len(DMI_Timeseries)):
        f.write(f'{(DMI_Dates[i] - ref_date).days}, {DMI_Timeseries[i]}\n')
    f.close()

    return None

def process_and_write_data() -> None:
    '''
        Puts everything together to produce the climate change corrected DMI
        data and save it to a file in the ./DATA/ folder.

        Parameters:
            None

        Inputs:
            None
    '''

    #bring in and load the data
    h_dates,h_lat,h_lon,h_sst = load_HadISST_data()
    #correct the sst data
    corrected_sst = correction_loop(h_sst)
    #get the inds for the monthly dates
    monthly_date_inds = date_grouper(h_dates)
    #get the monthly means
    monthly_means = get_monthly_means(monthly_date_inds,corrected_sst)
    #get the anomalous sst
    anomalous_sst = get_anomalous_ssts(monthly_date_inds,monthly_means,corrected_sst)
    #put the data into IOD East and IOD West boxes
    IOD_EAST_LATS,IOD_EAST_LONS,IOD_EAST_DATA = box_data(IOD_EAST_BOX,h_lat,h_lon,anomalous_sst)
    IOD_WEST_LATS,IOD_WEST_LONS,IOD_WEST_DATA = box_data(IOD_WEST_BOX,h_lat,h_lon,anomalous_sst)
    #take the area average in each box
    AA_IOD_EAST = area_average(IOD_EAST_LATS,IOD_EAST_LONS,IOD_EAST_DATA)
    AA_IOD_WEST = area_average(IOD_WEST_LATS,IOD_WEST_LONS,IOD_WEST_DATA)
    # take the difference to produce the timeseires
    DMI_Timeseries = AA_IOD_WEST - AA_IOD_EAST
    # now write the data to the file
    write_data(DMI_Timeseries,h_dates)

    return None

###### Main Handling #######
if __name__ == '__main__':
    pass

