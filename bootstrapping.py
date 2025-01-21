'''
 This script handles the bootstrapping experiments for Hudson and Miller 2025:
 A Curated Database of Milky Seas Since 1600.
'''

# Imports go here
# Standard Libraries
import datetime as dt
from IPython.display import clear_output
import time
# Public Libraries
import numpy as np
#import matplotlib.pyplot as plt
import datetime as dt
import numpy.ma as ma
from netCDF4 import Dataset
import pandas as pd

###### PATHS GO HERE ######
db_file = './DATA/Supplemental_2.tsv'
dmi_file = './DATA/CC_CORRECTED_DMI.txt'
wy_file = './DATA/WY_ERA5_850hPA_200hPa_winds.nc'
au_file = './DATA/AUSMI_ERA5_850hPA_winds.nc'
n34_file = './DATA/N34_Index.txt'

###### OTHER GLOBAL VARIABLES GO HERE ######
JAVA_BOX = [90,-20,120,-3] #left,bottom,right,top
BANDA_BOX = [120.5,-15,140,0]
AS_BOX = [40,-4,75,30]

###### Event Class Goes Here ######
class Event:
    '''
        Defines a class that combines observations int milky sea events, whether that be from
        eyewitness accounts, satellite observations, or both.

    '''
    def __init__(self,start_date,end_date,start_lat,start_lon,start_code,start_area,max_day,max_dist) -> None:
        '''
            Initializes an Event object using a single observation from the database and the
            thresholds for adding new observations to the event.

            Parameters:
                start_date (dt.datetime): The timestamp corresponding to the start of the observation
                end_date (dt.datetime): The timestampe corresponding to the end of the observation if it exists
                start_lat (float): The latitude corresponding to the observation [-90 -> 90]
                start_lon (float): The longitude corresponding to the observation [-180 -> 180]
                start_code (int): The confidence code for the observation from the database
                start_area (float): The satellite derived area if it exists in the database
                max_day (int): The threshold for the maximum number of days for new observations to be added to the event
                max_dist (float): The threshold for the maximum distance for new observations to be added to the event
            
            Returns:
                None
        '''
        # set the start/end date, if no end date exists just use the start date
        self.start_date = start_date
        if type(end_date) == dt.datetime:
            self.end_date = end_date
        else:
            self.end_date = start_date
        # make an array containing all dates associated with the event
        self.dates = [self.start_date,self.end_date]
        # arrays for holding all locations associated with the events
        self.lons = [start_lon]
        self.lats = [start_lat]
        # array for the confidence codes assoced with all obs in an event
        self.ev_codes = [start_code]
        # track how many obs have been merged into this event
        self.num_obs = 1
        # set the thresholds for grouping new observations into the event
        self.max_day = max_day
        self.max_dist = max_dist
        # if there is a satellite derived area add it to the event, otherwise use a placeholder
        if start_area == 'N/A':
            self.area = -9999
            self.area_type = None
        else:
            self.area = float(start_area)
            self.area_type = 'Satellite'
    
    def distance(self,new_lat,new_lon) -> float:
        '''
            Calculates the distance in degrees between a new observation and the original observation
            used to create the event

            Parameters:
                new_lat (float): The latitude of the new observation
                new_lon (float): The longitude of the new observation
            
            Returns:
                dist (float): The distance between the new observation and the original observation
        '''
        return np.sqrt( (self.lons[-1] - new_lon)**2 + (self.lats[-1] - new_lat)**2)
    
    def consider_new_ob(self,new_start,new_end,new_lat,new_lon,ev_code,area) -> bool:
        '''
            Considers whether or not a new observation should be merged into this event.

            Parameters:
                new_start (dt.datetime): The timestamp corresponding to the start of the new observation
                new_end (dt.datetime): The timestamp corresponding to the end of the new observation
                new_lat (float): The latitude of the new observation
                new_lon (float): The longitude of the new observation
                ev_code (int): The confidence code of the new observation
                area (float): The satellite derived area of the new observation if it exists

            Returns:
                should_be_added (bool): Whether or not the new observation should be added to the event
        '''
        if np.abs((new_start - self.start_date).days) <= self.max_day and self.distance(new_lat,new_lon) <= self.max_dist:
            self.add_new_obs_details(new_start,new_end,new_lat,new_lon,ev_code,area)
            return True
        else:
            return False

    def add_new_obs_details(self,new_start,new_end,new_lat,new_lon,ev_code,area) -> None:
        '''
            Merges a new observation with the event.

            Parameters:
                new_start (dt.datetime): The timestamp corresponding to the start of the new observation
                new_end (dt.datetime): The timestamp corresponding to the end of the new observation
                new_lat (float): The latitude of the new observation
                new_lon (float): The longitude of the new observation
                ev_code (int): The confidence code of the new observation
                area (float): The satellite derived area of the new observation if it exists
            
            Returns:
                None
        '''
        # Add the new dates to self.dates
        self.dates.append(new_start)
        if type(new_end) == dt.datetime:
            self.dates.append(new_end)
        # Add the new lat/lons to self.lats and self.lons
        self.lons.append(new_lon)
        self.lats.append(new_lat)
        # Add the new ev_codes to self.ev_codes
        self.ev_codes.append(ev_code)
        # Increase the number of observations in the event
        self.num_obs += 1
        # Sort the dates by time and make the latest date the new self.end_date
        self.dates = sorted(self.dates)
        self.end_date = self.dates[-1]
        # If a satellite area is found and one isn't currently being used make it the new area.
        if area != 'N/A' and self.area == -9999:
            self.area = float(area)

###### Year Class Goes Here ######
class Year:
    def __init__(self,ob_months):
        self.ob_months = ob_months
        if len(ob_months) == 0:
            self.ob_months = [-9999]

###### FUNCTIONS GO HERE ######
def load_database() -> tuple[np.ndarray]:
    db = pd.read_csv(db_file, sep = '\t')

    db_start_dates = db['Observation Start Date'].values
    db_end_dates = db['Observation End Date'].values
    db_lats = db['Approximate Lat'].values
    db_lons = db['Approximate Lon'].values
    db_confidence = db['Confidence'].values
    db_area = db['Area KM2'].values

    return db_start_dates,db_end_dates,db_lats,db_lons,db_confidence,db_area

def convert_separate_monthly_data(years,jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec):
    start_year = int(years[0])
    end_year = int(years[-1] + 1)

    #make an array for the dates
    date_array =[]
    for i in range(end_year-start_year):
        for j in range(12):
            date_array.append(dt.datetime(start_year+i,j+1,1))
    
    date_array = np.array(date_array)
    #convert the monthly data arrays to one long array
    all_data_array = []
    for i in range(len(jan)):
        all_data_array.append(jan[i])
        all_data_array.append(feb[i])
        all_data_array.append(mar[i])
        all_data_array.append(apr[i])
        all_data_array.append(may[i])
        all_data_array.append(jun[i])
        all_data_array.append(jul[i])
        all_data_array.append(aug[i])
        all_data_array.append(sep[i])
        all_data_array.append(oct[i])
        all_data_array.append(nov[i])
        all_data_array.append(dec[i])
    all_data_array = np.array(all_data_array)
  
    return date_array,all_data_array

def load_n34_data() -> tuple[np.ndarray]:
    n34_year,n34_jan,n34_feb,n34_mar,n34_apr,n34_may,n34_jun,n34_jul,n34_aug,n34_sep,n34_oct,n34_nov,n34_dec = np.loadtxt(n34_file,unpack = True)
    n34_dates,n34_data = convert_separate_monthly_data(n34_year,n34_jan,n34_feb,n34_mar,n34_apr,n34_may,n34_jun,n34_jul,n34_aug,n34_sep,n34_oct,n34_nov,n34_dec)

    return n34_dates,n34_data

def load_dmi_data() -> tuple[np.ndarray]:
    dmi_days_since_jan_1_1870,dmi_data = np.loadtxt(dmi_file,delimiter=',',skiprows=5,unpack=True)
    dmi_dates = np.array([dt.datetime(1870,1,1) + dt.timedelta(days = d) for d in dmi_days_since_jan_1_1870])
    # change to being the first of month for compliance with rest of code
    dmi_dates = np.array([dt.datetime(did.year,did.month,1) for did in dmi_dates])

    return dmi_dates,dmi_data

def load_wy_data() -> tuple[np.ndarray]:

    wy_lons = ma.getdata(Dataset(wy_file).variables['longitude'][:])
    wy_lats = ma.getdata(Dataset(wy_file).variables['latitude'][:])
    wy_time = ma.getdata(Dataset(wy_file).variables['valid_time'][:]) #hours since 1900-01-01
    wy_u200 = ma.getdata(Dataset(wy_file).variables['u'][:,0,:,:])
    wy_u850 = ma.getdata(Dataset(wy_file).variables['u'][:,1,:,:])
    wy_dates = np.array([dt.datetime.fromtimestamp(time) + dt.timedelta(hours=13) for time in wy_time])
    wy_dates = np.array([dt.datetime(wyd.year,wyd.month,1) for wyd in wy_dates])
    wy_index = np.nanmean(wy_u850 - wy_u200,axis = (1,2))

    return wy_index,wy_dates,wy_lats,wy_lons

def load_ausmi_data() -> tuple[np.ndarray]:

    au_lons = ma.getdata(Dataset(au_file).variables['longitude'][:])
    au_lats = ma.getdata(Dataset(au_file).variables['latitude'][:])
    au_time = ma.getdata(Dataset(au_file).variables['valid_time'][:]) #hours since 1900-01-01
    au_u850 = ma.getdata(Dataset(au_file).variables['u'][:,:,:])
    au_dates = np.array([dt.datetime.fromtimestamp(time) + dt.timedelta(hours=13) for time in au_time])
    au_dates = np.array([dt.datetime(aud.year,aud.month,1) for aud in au_dates])
    au_index = np.nanmean(au_u850[:],axis = (1,2))

    return au_index,au_dates,au_lats,au_lons

def string_parser(loc_string:str) -> list:
    '''Parses a location string and returns the degrees, minutes, and direction that make it up.
    If the location string is S or W it will be converted to negative'''
    
    #0 is the number of degrees
    #1 is the word 'deg'
    #2 is the number of minutes followed by a '
    #3 is the N,S,E,W signifier

    split_string = loc_string.split(' ')

    num_degrees = float(split_string[0])
    num_minutes = float(split_string[2][:-1]) #this takes out the ' to signify minutes
    dir_signifier = split_string[3]

    return [num_degrees,num_minutes,dir_signifier]

def convert_minutes_to_degrees(num_degrees:float) -> float:
    '''Given a number of minutes convert that to the fraction of degrees it represents
    and then return that value as a float.'''

    if num_degrees >= 60:
        raise ValueError("Number of Minutes in a degree can't be greater than 60.")
    elif num_degrees < 0:
        raise ValueError("Can't have negative minutes in a degree.")

    return num_degrees/60.

def get_direction_sign(dir_signifier:str) -> int:
    '''Based on the direction (N,S,E,W) return +/- 1 to multiply final result by.'''

    if dir_signifier == 'W' or dir_signifier == 'S':
        return -1
    else:
        return 1

def get_loc_in_degrees(loc_string:str) -> float:
    '''Takes in the locaction string and returns the coordinate on earth as a decimal -180 -> 180.'''
    parsed_string_list = string_parser(loc_string)
    minutes_in_decimal = convert_minutes_to_degrees(parsed_string_list[1])
    direction_sign = get_direction_sign(parsed_string_list[2])

    return direction_sign * (parsed_string_list[0]+minutes_in_decimal)

def convert_dates(start_date_array:np.ndarray,end_date_array:np.ndarray) -> list:
    '''
        Converts dates from the database into a list of datetime objects
        to make them more usable.

        Parameters:
            start_date_array (np.ndarray): Numpy array containing the start dates as strings
            end_date_array (np.ndarray): Numpy array containing the end dates as strings
        
        Returns:
            start_conv_dates (list): The start dates converted to datetime objects
            end_conv_dates (list): The end dates converted to datetime objects
    '''

    # fill in missing end dates with the start date
    missing_end_dates = np.where(end_date_array == '(?)')
    end_date_array[missing_end_dates] = start_date_array[missing_end_dates]

    # convert the dates
    start_conv_dates = np.array([dt.datetime.strptime(d,"%m/%d/%Y") for d in start_date_array])
    end_conv_dates = np.array([dt.datetime.strptime(d,'%m/%d/%Y') for d in end_date_array])

    return start_conv_dates,end_conv_dates

def get_box_obs(box:list,combined_lats:np.ndarray,combined_lons:np.ndarray) -> np.ndarray:
    '''
        Returns the indices needed to divide the observations into those
        that do and do not belong within a specified box

        Parameters:
            box (list): A list defining the box (west, south, east, north)
            combined_lats (np.ndarray): The lats of all obs combined
            combined_lons (np.ndarray): The lons of all obs combined

        Returns:
            loc_match (np.ndarray): A numpy array containing the indicies which
                correspond to which obs go in the box
    '''
    #get the overlaps of the lats and lons with the box and return the matching inds
    west_match = np.where(combined_lons >= box[0])
    east_match = np.where(combined_lons <= box[2])
    south_match = np.where(combined_lats >= box[1])
    north_match = np.where(combined_lats <= box[3])

    lon_match = np.intersect1d(west_match,east_match)
    lat_match = np.intersect1d(south_match,north_match)

    loc_match = np.intersect1d(lon_match,lat_match)

    return loc_match

def limit_ROI_data_to_TOI(lats:np.ndarray,lons:np.ndarray,start_dates:np.ndarray,
                          end_dates:np.ndarray,codes:np.ndarray,areas:np.ndarray,
                          season:np.ndarray) -> tuple[np.ndarray]:
    '''
        Limits data refined to a region of interest (ROI) to our specific
        time of interest (TOI) for the bootstrapping

        Parameters:
            lats (np.ndarray): The latitudes for the observations
            lons (np.ndarray): The longitudes for the observations
    '''
    new_lats = []
    new_lons = []
    new_start_dates = []
    new_end_dates = []
    new_codes = []
    new_areas = []
    for i in range(len(start_dates)):
        if start_dates[i].year >= 1950 and start_dates[i].year < 2000:
            if season == 'summer' and start_dates[i].month >= 6 and start_dates[i].month <= 10:
                new_lats.append(lats[i])
                new_lons.append(lons[i])
                new_start_dates.append(start_dates[i])
                new_end_dates.append(end_dates[i])
                new_codes.append(codes[i])
                new_areas.append(areas[i])
            elif season == 'winter' and (start_dates[i].month == 12 or start_dates[i].month <= 3):
                new_lats.append(lats[i])
                new_lons.append(lons[i])
                new_start_dates.append(start_dates[i])
                new_end_dates.append(end_dates[i])
                new_codes.append(codes[i])
                new_areas.append(areas[i])
    
    return np.array(new_lats),np.array(new_lons),np.array(new_start_dates),np.array(new_end_dates),np.array(new_codes), np.array(new_areas)

def classify_n34_months(n34_data:np.ndarray) -> np.ndarray:

    n34_month_classes = np.zeros(len(n34_data))
    streak_counter = 0
    i_counter = 0
    while i_counter < len(n34_data):
        if np.abs(n34_data[i_counter]) >= 0.5:
            streak_counter = 1
            i_sign = np.sign(n34_data[i_counter])
            for j in range(i_counter+1,len(n34_data)):
                if np.sign(n34_data[j]) == i_sign and np.abs(n34_data[j]) >= 0.5:
                    streak_counter += 1
                elif streak_counter >= 4:
                    n34_month_classes[i_counter:i_counter+streak_counter] = i_sign
                    i_counter = j-1
                    streak_counter = 0
                    break
                else:
                    streak_counter = 0
                    i_counter = j-1
                    break
        i_counter+=1

    return n34_month_classes

def classify_dmi_months(dmi_data:np.ndarray) -> np.ndarray:

    dmi_month_classes = np.zeros(len(dmi_data))
    for i in range(len(dmi_data)):
        if dmi_data[i] >= 0.25:
            dmi_month_classes[i] = 1
        elif dmi_data[i] <= -0.25:
            dmi_month_classes[i] = -1

    return dmi_month_classes

def classify_monsoon_years(wy_index,wy_dates,au_index,au_dates) -> tuple[np.ndarray]:

    wy_index_JJAS_mean_yearly = []
    wy_index_DJFM_mean_yearly = []
    au_index_JJAS_mean_yearly = []

    wy_years = np.arange(1950,2024,1)
    au_years = np.arange(1950,2024,1)

    for i in range(len(wy_years)):
        wy_year_list_sum = []
        wy_year_list_wint = []
        au_year_list = []
        for j in range(len(wy_index)):
            if (wy_dates[j].month >= 6 and wy_dates[j].month <=9) and (wy_dates[j].year == wy_years[i]):
                wy_year_list_sum.append(wy_index[j])
            if (wy_dates[j].month == 12 and wy_dates[j].year == wy_years[i]-1) or (wy_dates[j].month <= 3 and wy_dates[j].year == wy_years[i]):
                wy_year_list_wint.append(wy_index[j])
            if (au_dates[j].month >= 6 and au_dates[j].month <=9) and (au_dates[j].year == au_years[i]):
                au_year_list.append(au_index[j])  
        wy_year_mean_sum = np.mean(wy_year_list_sum)
        wy_year_mean_wint = np.mean(wy_year_list_wint)
        au_year_mean = np.mean(au_year_list)
        wy_index_JJAS_mean_yearly.append(wy_year_mean_sum)
        wy_index_DJFM_mean_yearly.append(-1 * wy_year_mean_wint) #WY is negative in boreal winter so multiply by -1 here to make life easier
        au_index_JJAS_mean_yearly.append(-1*au_year_mean)#AUSMI is negative in boreal summer so multiply by -1 here to make life easier

    threshold = 0.2
    wy_sum_bottom_quintile = np.quantile(wy_index_JJAS_mean_yearly,threshold)
    wy_sum_top_quintile = np.quantile(wy_index_JJAS_mean_yearly,1-threshold)
    wy_wint_bottom_quintile = np.quantile(wy_index_DJFM_mean_yearly,threshold)
    wy_wint_top_quintile = np.quantile(wy_index_DJFM_mean_yearly,1-threshold)
    au_top_quintile = np.quantile(au_index_JJAS_mean_yearly,1-threshold)
    au_bottom_quintile = np.quantile(au_index_JJAS_mean_yearly,threshold)

    wy_yearly_classification_sum = np.zeros(len(wy_years))
    wy_yearly_classification_wint = np.zeros(len(wy_years))
    au_yearly_classification = np.zeros(len(au_years))
    for i in range(len(wy_index_JJAS_mean_yearly)):
        if wy_index_JJAS_mean_yearly[i] <= wy_sum_bottom_quintile:
            wy_yearly_classification_sum[i] = -1
        elif wy_index_JJAS_mean_yearly[i] >= wy_sum_top_quintile:
            wy_yearly_classification_sum[i] = 1
        
        if wy_index_DJFM_mean_yearly[i] <= wy_wint_bottom_quintile:
            wy_yearly_classification_wint[i] = -1
        elif wy_index_DJFM_mean_yearly[i] >= wy_wint_top_quintile:
            wy_yearly_classification_wint[i] = 1

        if au_index_JJAS_mean_yearly[i] <= au_bottom_quintile:
            au_yearly_classification[i] = -1
        elif au_index_JJAS_mean_yearly[i] >= au_top_quintile:
            au_yearly_classification[i] = 1

    return wy_yearly_classification_sum,wy_yearly_classification_wint,au_yearly_classification


def random_order(year_holder):
    rng = np.random.default_rng(np.random.randint(0,10000000,size = 1))
    #make the inds
    inds = np.arange(0,len(year_holder),1)
    #shuffle them using the Fisher-Yates Shuffle (Uniform)
    rng.shuffle(inds)

    return inds

def scramble_year_holder(year_holder,random_inds):
    new_year_holder = np.empty(len(year_holder),dtype=object)
    for i in range(len(new_year_holder)):
        new_year_holder[i] = year_holder[random_inds[i]]
    return new_year_holder

# The ones designed to take the state of the environmental state given an event's date
def get_year_counts_ENSO(year_order,year_holder,dates,month_classes):
    #make variables to hold the counts
    POS_Counts = 0
    NEG_Counts = 0
    NEU_Counts = 0
    #loop over the year holder and compare the dates against n34_dates
    for i in range(len(year_holder)):
        if not year_holder[i].ob_months[0] == -9999:
            ob_classes = np.empty(5)
            for k in range(5):
                if k <= 1:
                    ob_date = dt.datetime(int(year_order[i])-1,12-k,1) #look at the state at the start of the year
                    ob_classes[k] = month_classes[np.where(dates == ob_date)[0][0]]
                else:
                    ob_date = dt.datetime(int(year_order[i]),k-1,1)
                    ob_classes[k] = month_classes[np.where(dates == ob_date)[0][0]]

            for j in range(len(year_holder[i].ob_months)):
                ob_class = np.mean(ob_classes)
                if ob_class >= 0.5:
                    POS_Counts+=1
                elif ob_class <= -0.5:
                    NEG_Counts+=1
                else:
                    NEU_Counts+=1
    return POS_Counts,NEU_Counts,NEG_Counts

def get_end_year_counts_ENSO(year_order,year_holder,dates,month_classes):
    #make variables to hold the counts
    POS_Counts = 0
    NEG_Counts = 0
    NEU_Counts = 0
    #loop over the year holder and compare the dates against n34_dates
    for i in range(len(year_holder)):
        if not year_holder[i].ob_months[0] == -9999:
            ob_classes = np.empty(5)
            for k in range(5):
                if k <= 1:
                    ob_date = dt.datetime(int(year_order[i]),12-k,1) #look at the state at the end of the year
                    ob_classes[k] = month_classes[np.where(dates == ob_date)[0][0]]
                else:
                    ob_date = dt.datetime(int(year_order[i])+1,k-1,1) #the start of the next year
                    ob_classes[k] = month_classes[np.where(dates == ob_date)[0][0]]

            for j in range(len(year_holder[i].ob_months)):
                ob_class = np.mean(ob_classes)
                if ob_class >= 0.5:
                    POS_Counts+=1
                elif ob_class <= -0.5:
                    NEG_Counts+=1
                else:
                    NEU_Counts+=1
    return POS_Counts,NEU_Counts,NEG_Counts

def get_year_counts_DMI(year_order,year_holder,dates,month_classes):
    #make variables to hold the counts
    POS_Counts = 0
    NEG_Counts = 0
    NEU_Counts = 0
    #loop over the year holder and compare the dates against n34_dates
    for i in range(len(year_holder)):
        if not year_holder[i].ob_months[0] == -9999:
            ob_classes = np.empty(5)
            for k in range(5):
                ob_date = dt.datetime(int(year_order[i]),k+6,1) #look at the state at the state between June-Oct
                ob_classes[k] = month_classes[np.where(dates == ob_date)[0][0]]
            ob_class = np.mean(ob_classes)
            for j in range(len(year_holder[i].ob_months)):
                if ob_class >= 0.5:
                    POS_Counts+=1
                elif ob_class < 0.5 and ob_class > -0.5:
                    NEU_Counts+=1
                elif ob_class <= -0.5:
                    NEG_Counts+=1
    return POS_Counts,NEU_Counts,NEG_Counts

def get_year_counts_WY_and_AUSMI(year_holder,yearly_classifications):
    #make variables to hold the counts
    STR_Counts = 0
    AVE_Counts = 0
    WEAK_Counts = 0
    #loop over the year holder and get the counts
    for i in range(len(year_holder)):
        if not year_holder[i].ob_months[0] == -9999:
            for j in range(len(year_holder[i].ob_months)):
                ob_class = yearly_classifications[i]
                if ob_class == 1:
                    STR_Counts += 1
                elif ob_class == -1:
                    WEAK_Counts += 1
                else:
                    AVE_Counts += 1

    return STR_Counts,AVE_Counts,WEAK_Counts

#now to make the loops which do these repeatedly based on the N I set
def ENSO_and_DMI_loop(N,year_order,year_holder,dates,month_classes,event):
    #the counts for ENSO or IOD depending on which month_classes I slot in
    POS_Counts = np.empty(N)
    NEU_Counts = np.empty(N)
    NEG_Counts = np.empty(N)

    for i in range(N):
        random_order_inds = random_order(year_holder)
        randomized_year_holder = scramble_year_holder(year_holder,random_order_inds)
        if event == 'ENSO':
            POS,NEU,NEG = get_year_counts_ENSO(year_order,randomized_year_holder,dates,month_classes)
        elif event == 'ENSO_END':
            POS,NEU,NEG = get_end_year_counts_ENSO(year_order,randomized_year_holder,dates,month_classes)
        elif event == 'IOD':
            POS,NEU,NEG = get_year_counts_DMI(year_order,randomized_year_holder,dates,month_classes)
        POS_Counts[i] = POS
        NEU_Counts[i] = NEU
        NEG_Counts[i] = NEG
    
    return POS_Counts,NEU_Counts,NEG_Counts

def Monsoon_loop(N,year_holder,yearly_classification):
    STR_Counts = np.empty(N)
    AVE_Counts = np.empty(N)
    WEAK_Counts = np.empty(N)

    for i in range(N):
        random_order_inds = random_order(year_holder)
        randomized_year_holder = scramble_year_holder(year_holder,random_order_inds)
        STR,AVE,WEAK = get_year_counts_WY_and_AUSMI(randomized_year_holder,yearly_classification)
        STR_Counts[i] = STR
        AVE_Counts[i] = AVE
        WEAK_Counts[i] = WEAK
    
    return STR_Counts,AVE_Counts,WEAK_Counts

#now the function that takes the data and turns it into events
def Event_Grouper(ob_lats,ob_lons,ob_start_dates,ob_end_dates,ob_codes,ob_areas,day_dist,deg_dist):
    event_list = []
    event_coding = np.zeros(len(ob_start_dates))
    for i in range(len(ob_start_dates)):
        if i == 0:
            new_ev = Event(ob_start_dates[i],ob_end_dates[i],ob_lats[i],ob_lons[i],ob_codes[i],ob_areas[i],day_dist,deg_dist)
            event_coding[i] = 1
            event_list.append(new_ev)
        else:
            if event_coding[i] == 0:
                new_ev = Event(ob_start_dates[i],ob_end_dates[i],ob_lats[i],ob_lons[i],ob_codes[i],ob_areas[i],day_dist,deg_dist)
                event_coding[i] = 1
                event_list.append(new_ev)
        for j in range(i+1,len(ob_start_dates)):
            if event_coding[j] == 0:
                if event_list[-1].consider_new_ob(ob_start_dates[j],ob_end_dates[j],ob_lats[j],ob_lons[j],ob_codes[j],ob_areas[j]):
                    event_coding[j] = 1
    
    return event_list

def Events_to_Year_Objects(event_list):
    def_year_order = np.arange(1950,2000,1)
    def_year_holder = np.empty(len(def_year_order),dtype = object)
    for i in range(len(def_year_order)):
        year_list = []
        for j in range(len(event_list)):
            if event_list[j].dates[0].year == def_year_order[i] or (event_list[j].dates[0].month == 12 and event_list[j].dates[0].year == def_year_order[i]-1):
                year_list.append(event_list[j].dates[0].month)
        def_year_holder[i] = Year(year_list)
    
    return def_year_order,def_year_holder

#the function the records how significant a result is:
def Get_Truth_Sig_From_Distribution(truth_val,distribution):
    '''
        Given the bootstrap distribution determines how far away the observed
        truth value froms from the distribution and records a code of 1,2,3, or 4
        if it exceeds the 80%,85%,90%, or 95% confidence interval. Codes are negative
        if it is below the mean. Will return a 0 if the result is not significant at any
        of those levels.

        Inputs:
        truth_val(float): The value actually observed in reality
        distribution (np.ndarray): The distribution gained through the bootstrapping.

        Output:
        event_code (int): How significant the truth_val is
    '''
    dist_mean = np.nanmean(distribution)
    dist_std = np.nanstd(distribution)

    #make the thresholds
    #first the ones for significantly more
    above_80 = dist_mean + 1.28*dist_std
    above_85 = dist_mean + 1.44*dist_std
    above_90 = dist_mean + 1.68*dist_std
    above_95 = dist_mean + 1.96*dist_std
    #now for significantly below
    below_80 = dist_mean - 1.28*dist_std
    below_85 = dist_mean - 1.44*dist_std
    below_90 = dist_mean - 1.68*dist_std
    below_95 = dist_mean - 1.96*dist_std

    #now to see what the truth_val falls
    if truth_val > dist_mean:
        if truth_val > above_95:
            return 4
        elif truth_val > above_90:
            return 3
        elif truth_val > above_85:
            return 2
        elif truth_val > above_80:
            return 1
        else:
            return 0
    elif truth_val <= dist_mean:
        if truth_val < below_95:
            return -4
        elif truth_val < below_90:
            return -3
        elif truth_val < below_85:
            return -2
        elif truth_val < below_80:
            return -1
        else:
            return 0

    return -9999

def bootstrap_loop(N):
    # Load in the data for the boostrapping
    db_start_dates,db_end_dates,db_lats,db_lons,db_confidence,db_area = load_database()
    n34_dates,n34_data = load_n34_data()
    dmi_dates,dmi_data = load_dmi_data()
    wy_index,wy_dates,wy_lats,wy_lons = load_wy_data()
    au_index,au_dates,au_lats,au_lons = load_ausmi_data()
    # transform the db data into the forms I want
    db_conv_lons = np.array([get_loc_in_degrees(dbl) for dbl in db_lons])
    db_conv_lats = np.array([get_loc_in_degrees(dbl) for dbl in db_lats])
    db_conv_start_dates,db_conv_end_dates = convert_dates(db_start_dates,db_end_dates)
    #split the db_data into the boxes I defined at the start
    java_ins = get_box_obs(JAVA_BOX,db_conv_lats,db_conv_lons)
    as_all_ins = get_box_obs(AS_BOX,db_conv_lats,db_conv_lons)
    banda_ins = get_box_obs(BANDA_BOX,db_conv_lats,db_conv_lons)
    #split the as_all_ins into as_sum_ins and as_wint_ins
    #boreal summer is June-Oct, boreal winter is Dec-Mar
    as_sum_ins = []
    as_wint_ins = []
    for i in range(len(as_all_ins)):
        ind_date = db_conv_start_dates[as_all_ins[i]]
        if ind_date.month >= 6 and ind_date.month <= 10:
            as_sum_ins.append(as_all_ins[i])
        elif ind_date.month == 12 or ind_date.month <= 3:
            as_wint_ins.append(as_all_ins[i])
    as_sum_ins = np.array(as_sum_ins)
    as_wint_ins = np.array(as_wint_ins)
    #now that I have these indices make arrays of lats/lons/dates/ev_codes for each ROI
    #first the AS Summer
    as_sum_lats = db_conv_lats[as_sum_ins]
    as_sum_lons = db_conv_lons[as_sum_ins]
    as_sum_start_dates = db_conv_start_dates[as_sum_ins]
    as_sum_end_dates = db_conv_end_dates[as_sum_ins]
    as_sum_codes = db_confidence[as_sum_ins]
    as_sum_areas = db_area[as_sum_ins]
    #next AS Winter
    as_wint_lats = db_conv_lats[as_wint_ins]
    as_wint_lons = db_conv_lons[as_wint_ins]
    as_wint_start_dates = db_conv_start_dates[as_wint_ins]
    as_wint_end_dates = db_conv_end_dates[as_wint_ins]
    as_wint_codes = db_confidence[as_wint_ins]
    as_wint_areas = db_area[as_wint_ins]
    #now Banda
    banda_lats = db_conv_lats[banda_ins]
    banda_lons = db_conv_lons[banda_ins]
    banda_start_dates = db_conv_start_dates[banda_ins]
    banda_end_dates = db_conv_end_dates[banda_ins]
    banda_codes = db_confidence[banda_ins]
    banda_areas = db_area[banda_ins]
    #finally Java
    java_lats = db_conv_lats[java_ins]
    java_lons = db_conv_lons[java_ins]
    java_start_dates = db_conv_start_dates[java_ins]
    java_end_dates = db_conv_end_dates[java_ins]
    java_codes = db_confidence[java_ins]
    java_areas = db_area[java_ins]

    #now limit to my TOI
    bs_as_sum_lats,bs_as_sum_lons,bs_as_sum_start_dates,bs_as_sum_end_dates,bs_as_sum_codes,bs_as_sum_areas = limit_ROI_data_to_TOI(as_sum_lats,as_sum_lons,as_sum_start_dates,as_sum_end_dates,as_sum_codes,as_sum_areas,'summer')
    bs_as_wint_lats,bs_as_wint_lons,bs_as_wint_start_dates,bs_as_wint_end_dates,bs_as_wint_codes,bs_as_wint_areas = limit_ROI_data_to_TOI(as_wint_lats,as_wint_lons,as_wint_start_dates,as_wint_end_dates,as_wint_codes,as_wint_areas,'winter')
    bs_banda_lats,bs_banda_lons,bs_banda_start_dates,bs_banda_end_dates,bs_banda_codes,bs_banda_areas = limit_ROI_data_to_TOI(banda_lats,banda_lons,banda_start_dates,banda_end_dates,banda_codes,banda_areas,'summer')
    bs_java_lats,bs_java_lons,bs_java_start_dates,bs_java_end_dates,bs_java_codes,bs_java_areas = limit_ROI_data_to_TOI(java_lats,java_lons,java_start_dates,java_end_dates,java_codes,java_areas,'summer')

    #now classify months and years
    n34_month_classes = classify_n34_months(n34_data)
    dmi_month_classes = classify_dmi_months(dmi_data)
    wy_yearly_classification_sum,wy_yearly_classification_wint,au_yearly_classification = classify_monsoon_years(wy_index,wy_dates,au_index,au_dates)

    #Define how fine a setting grid I use
    day_distances = np.array([2,3,4,5,7,10,15,20,30,40,50,60,70,80,90])
    deg_distances = np.array([0.25,0.5,0.75,1,1.5,2,3,4,5,6,7,8])
    #Define the arrays that hold the codes
    # 0 is for Arabian Sea Summer, 1 is for Arabian Sea Winter
    # 2 is for Banda Sea, 3 is for Java
    #ENSO CODES
    el_nino_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    neu_enso_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    la_nina_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    #ENSO end of year codes
    en_end_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    ne_end_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    ln_end_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    #IOD CODES
    pos_iod_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    neu_iod_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    neg_iod_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    #SW MONSOON CODES
    str_sw_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    ave_sw_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    weak_sw_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    #NE MONSOON CODES
    str_ne_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    ave_ne_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    weak_ne_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    #SE MONSOON CODES
    str_se_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    ave_se_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    weak_se_monsoon_codes = np.zeros((4,len(day_distances),len(deg_distances)))
    #Define arrays that have the lats,lons,dates, and event_codes for all 4 regions
    #so that data can be accessed with similar indexing as the significance code arrays above
    regional_lats = np.empty(4,dtype = np.ndarray)
    regional_lons = np.empty(4,dtype = np.ndarray)
    regional_start_dates = np.empty(4,dtype = np.ndarray)
    regional_end_dates = np.empty(4,dtype = np.ndarray)
    regional_codes = np.empty(4,dtype = np.ndarray)
    regional_areas = np.empty(4,dtype = np.ndarray)
    #fill these arrays
    regional_lats[0] = bs_as_sum_lats[:]
    regional_lats[1] = bs_as_wint_lats[:]
    regional_lats[2] = bs_banda_lats[:]
    regional_lats[3] = bs_java_lats[:]

    regional_lons[0] = bs_as_sum_lons[:]
    regional_lons[1] = bs_as_wint_lons[:]
    regional_lons[2] = bs_banda_lons[:]
    regional_lons[3] = bs_java_lons[:]

    regional_start_dates[0] = bs_as_sum_start_dates[:]
    regional_start_dates[1] = bs_as_wint_start_dates[:]
    regional_start_dates[2] = bs_banda_start_dates[:]
    regional_start_dates[3] = bs_java_start_dates[:]

    regional_end_dates[0] = bs_as_sum_end_dates[:]
    regional_end_dates[1] = bs_as_wint_end_dates[:]
    regional_end_dates[2] = bs_banda_end_dates[:]
    regional_end_dates[3] = bs_java_end_dates[:]

    regional_codes[0] = bs_as_sum_codes[:]
    regional_codes[1] = bs_as_wint_codes[:]
    regional_codes[2] = bs_banda_codes[:]
    regional_codes[3] = bs_java_codes[:]

    regional_areas[0] = bs_as_sum_areas[:]
    regional_areas[1] = bs_as_wint_areas[:]
    regional_areas[2] = bs_banda_areas[:]
    regional_areas[3] = bs_java_areas[:]

    timer = time.time()
    for i in range(len(day_distances)): #my day settings for the event definition
        for j in range(len(deg_distances)): #my deg settings for the event definition
            for k in range(4): #my regions of interest
                #define the events for the region
                event_list = Event_Grouper(regional_lats[k],regional_lons[k],regional_start_dates[k],regional_end_dates[k],
                                            regional_codes[k],regional_areas[k],day_distances[i],deg_distances[j])
                #group them into year objects
                year_order,year_holder = Events_to_Year_Objects(event_list)
                #get the truth values for this event definition and region
                bs_elnino_truth,bs_neutral_enso_truth,bs_lanina_truth = get_year_counts_ENSO(year_order,year_holder,n34_dates,n34_month_classes)
                bs_el_end_truth,bs_neu_end_truth,bs_la_end_truth = get_end_year_counts_ENSO(year_order,year_holder,n34_dates,n34_month_classes)
                bs_posdmi_truth,bs_neutral_dmi_truth,bs_negdmi_truth = get_year_counts_DMI(year_order,year_holder,dmi_dates,dmi_month_classes)
                bs_str_sw_truth,bs_ave_sw_truth,bs_weak_sw_truth = get_year_counts_WY_and_AUSMI(year_holder[4:],wy_yearly_classification_sum)
                bs_str_ne_truth,bs_ave_ne_truth,bs_weak_ne_truth = get_year_counts_WY_and_AUSMI(year_holder[4:],wy_yearly_classification_wint)
                bs_str_se_truth,bs_ave_se_truth,bs_weak_se_truth = get_year_counts_WY_and_AUSMI(year_holder[4:],au_yearly_classification)
                #get the bootstrap definitions for this event definition and region
                elnino_dist,neutral_enso_dist,lanina_dist = ENSO_and_DMI_loop(N,year_order,year_holder,
                                                                            n34_dates,n34_month_classes,event = 'ENSO')
                elend_dist,neutral_enso_end_dist,laend_dist = ENSO_and_DMI_loop(N,year_order,year_holder,
                                                                            n34_dates,n34_month_classes,event = 'ENSO_END')
                posdmi_dist,neutral_dmi_dist,negdmi_dist = ENSO_and_DMI_loop(N,year_order,year_holder,
                                                                            dmi_dates,dmi_month_classes,event = 'IOD')
                str_sw_dist,ave_sw_dist,weak_sw_dist = Monsoon_loop(N,year_holder[4:],wy_yearly_classification_sum)
                str_ne_dist,ave_ne_dist,weak_ne_dist = Monsoon_loop(N,year_holder[4:],wy_yearly_classification_wint)
                str_se_dist,ave_se_dist,weak_se_dist = Monsoon_loop(N,year_holder[4:],au_yearly_classification)
                #slot the truth values into the code arrays
                #ENSO
                el_nino_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_elnino_truth,elnino_dist)
                la_nina_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_lanina_truth,lanina_dist)
                neu_enso_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_neutral_enso_truth,neutral_enso_dist)
                #ENSO END OF YEAR
                en_end_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_el_end_truth,elend_dist)
                ne_end_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_neu_end_truth,neutral_enso_end_dist)
                ln_end_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_la_end_truth,laend_dist)
                #IOD
                pos_iod_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_posdmi_truth,posdmi_dist)
                neu_iod_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_neutral_dmi_truth,neutral_dmi_dist)
                neg_iod_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_negdmi_truth,negdmi_dist)
                #SW Monsoon
                str_sw_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_str_sw_truth,str_sw_dist)
                ave_sw_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_ave_sw_truth,ave_sw_dist)
                weak_sw_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_weak_sw_truth,weak_sw_dist)
                #NE Monsoon
                str_ne_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_str_ne_truth,str_ne_dist)
                ave_ne_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_ave_ne_truth,ave_ne_dist)
                weak_ne_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_weak_ne_truth,weak_ne_dist)
                #SE Monsoon
                str_se_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_str_se_truth,str_se_dist)
                ave_se_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_ave_se_truth,ave_se_dist)
                weak_se_monsoon_codes[k,i,j] = Get_Truth_Sig_From_Distribution(bs_weak_se_truth,weak_se_dist)
            clear_output(wait=True)
            print(f'{(((i)/len(day_distances))+(.0666666*((j+1)/len(deg_distances))))*100:.2f}' +' percent done')
            print(f'Elapsed Time: {time.time() - timer:.0f} Seconds')

    clear_output(wait=True)
    print('100 percent done')
    print(f'Elapsed Time: {time.time() - timer:.0f} Seconds')

    np.save("./DATA/el_nino_codes.npy",el_nino_codes)
    np.save("./DATA/la_nina_codes.npy",la_nina_codes)
    np.save("./DATA/neutral_enso_codes.npy",neu_enso_codes)

    np.save("./DATA/eoy_el_nino_codes.npy",en_end_codes)
    np.save("./DATA/eoy_la_nina_codes.npy",ln_end_codes)
    np.save("./DATA/eoy_neutral_enso_codes.npy",ne_end_codes)

    np.save("./DATA/pos_iod_codes.npy",pos_iod_codes)
    np.save("./DATA/neu_iod_codes.npy",neu_iod_codes)
    np.save("./DATA/neg_iod_codes.npy",neg_iod_codes)

    np.save("./DATA/str_sw_monsoon_codes.npy",str_sw_monsoon_codes)
    np.save("./DATA/ave_sw_monsoon_codes.npy",ave_sw_monsoon_codes)
    np.save("./DATA/weak_sw_monsoon_codes.npy",weak_sw_monsoon_codes)

    np.save("./DATA/str_ne_monsoon_codes.npy",str_ne_monsoon_codes)
    np.save("./DATA/ave_ne_monsoon_codes.npy",ave_ne_monsoon_codes)
    np.save("./DATA/weak_ne_monsoon_codes.npy",weak_ne_monsoon_codes)

    np.save("./DATA/str_se_monsoon_codes.npy",str_se_monsoon_codes)
    np.save("./DATA/ave_se_monsoon_codes.npy",ave_se_monsoon_codes)
    np.save("./DATA/weak_se_monsoon_codes.npy",weak_se_monsoon_codes)

    return None

###### MAIN GOES HERE ######
if __name__ == '__main__':
    pass