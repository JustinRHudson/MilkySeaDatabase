'''
    This script downloads the ERA5 data used in the paper Hudson 
    and Miller 2025: A Curated Database of Milky Seas Since 1600.

    This script assumes you have set up the Copernicus Data Store API on your
    machine. If you haven't please follow the instructions at the link below:
    https://cds.climate.copernicus.eu/how-to-api

    This script also needs to be run in an environment with the cdsapi package
    available. The cdsapi_env.yml file can be used to create such an
    environment following the directions at if one is needed:
    https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file
'''

###### IMPORTS ######
# Standard Libraries #
import os
# Public Libraries #
import cdsapi

###### FUNCTIONS GO HERE ######
def retrieve_WY_winds() -> None:
    '''
        Uses a request to the ERA5 Copernicus Climate Data Store API
        to download the winds needed to create the Webster Yang Monsoon Index.

        The data is saved to ./DATA/WY_ERA5_850hPa_200hPa_winds.nc
    '''

    dataset = "reanalysis-era5-pressure-levels-monthly-means"
    request = {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": ["u_component_of_wind"],
        "pressure_level": ["200", "850"],
        "year": [
            "1950", "1951", "1952",
            "1953", "1954", "1955",
            "1956", "1957", "1958",
            "1959", "1960", "1961",
            "1962", "1963", "1964",
            "1965", "1966", "1967",
            "1968", "1969", "1970",
            "1971", "1972", "1973",
            "1974", "1975", "1976",
            "1977", "1978", "1979",
            "1980", "1981", "1982",
            "1983", "1984", "1985",
            "1986", "1987", "1988",
            "1989", "1990", "1991",
            "1992", "1993", "1994",
            "1995", "1996", "1997",
            "1998", "1999", "2000",
            "2001", "2002", "2003",
            "2004", "2005", "2006",
            "2007", "2008", "2009",
            "2010", "2011", "2012",
            "2013", "2014", "2015",
            "2016", "2017", "2018",
            "2019", "2020", "2021",
            "2022", "2023"
        ],
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "time": ["00:00"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [20, 40, 0, 110]
    }
    target = './DATA/WY_ERA5_850hPA_200hPa_winds.nc'
    client = cdsapi.Client()
    client.retrieve(dataset, request, target)

    return None

def retrieve_AUSMI_winds() -> None:
    '''
        Uses a request to the ERA5 Copernicus Climate Data Store API
        to download the winds needed to create the AUSMI Monsoon Index.

        The data is saved to ./DATA/AUSMI_ERA5_850hPa_winds.nc
    '''

    dataset = "reanalysis-era5-pressure-levels-monthly-means"
    request = {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": ["u_component_of_wind"],
        "pressure_level": ["850"],
        "year": [
            "1950", "1951", "1952",
            "1953", "1954", "1955",
            "1956", "1957", "1958",
            "1959", "1960", "1961",
            "1962", "1963", "1964",
            "1965", "1966", "1967",
            "1968", "1969", "1970",
            "1971", "1972", "1973",
            "1974", "1975", "1976",
            "1977", "1978", "1979",
            "1980", "1981", "1982",
            "1983", "1984", "1985",
            "1986", "1987", "1988",
            "1989", "1990", "1991",
            "1992", "1993", "1994",
            "1995", "1996", "1997",
            "1998", "1999", "2000",
            "2001", "2002", "2003",
            "2004", "2005", "2006",
            "2007", "2008", "2009",
            "2010", "2011", "2012",
            "2013", "2014", "2015",
            "2016", "2017", "2018",
            "2019", "2020", "2021",
            "2022", "2023"
        ],
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "time": ["00:00"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [-5, 110, -15, 130]
    }
    target = './DATA/AUSMI_ERA5_850hPA_winds.nc'
    client = cdsapi.Client()
    client.retrieve(dataset, request, target)

    return None

def retrieve_NWIO_Environmental_Winds() -> None:
    '''
        Uses a request to the ERA5 Copernicus Climate Data Store API
        to download the winds needed to create the environmental winds used in
        Figure 8.

        The data is saved to ./DATA/NWIO_Environmental_ERA5_1000hPa_winds.nc
    '''

    dataset = "reanalysis-era5-pressure-levels-monthly-means"
    request = {
        "product_type": ["monthly_averaged_reanalysis"],
        "variable": ["u_component_of_wind","v_component_of_wind"],
        "pressure_level": ["1000"],
        "year": [
            "1941", "1942", "1943",
            "1944", "1945", "1948",
            "1947", "1948", "1949",
            "1950", "1951", "1952",
            "1953", "1954", "1955",
            "1956", "1957", "1958",
            "1959", "1960", "1961",
            "1962", "1963", "1964",
            "1965", "1966", "1967",
            "1968", "1969", "1970",
            "1971", "1972", "1973",
            "1974", "1975", "1976",
            "1977", "1978", "1979",
            "1980", "1981", "1982",
            "1983", "1984", "1985",
            "1986", "1987", "1988",
            "1989", "1990", "1991",
            "1992", "1993", "1994",
            "1995", "1996", "1997",
            "1998", "1999", "2000",
            "2001", "2002", "2003",
            "2004", "2005", "2006",
            "2007", "2008", "2009",
            "2010", "2011", "2012",
            "2013", "2014", "2015",
            "2016", "2017", "2018",
            "2019", "2020", "2021",
            "2022", "2023"
        ],
        "month": [
            "01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12"
        ],
        "time": ["00:00"],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [20, 40, 0, 60]
    }
    target = './DATA/NWIO_Environmental_ERA5_1000hPa_winds.nc'
    client = cdsapi.Client()
    client.retrieve(dataset, request, target)

    return None

def main() -> None:
    if not os.path.exists('./DATA/WY_ERA5_850hPA_200hPa_winds.nc'):
        retrieve_WY_winds()
    else:
        print('Webster-Yang Winds Found. Moving to next file...')
    if not os.path.exists('./DATA/AUSMI_ERA5_850hPA_winds.nc'):
        retrieve_AUSMI_winds()
    else:
        print('AUSMI Winds Found. Moving to next file...')
    if not os.path.exists('./DATA/NWIO_Environmental_ERA5_1000hPa_winds.nc'):
        retrieve_NWIO_Environmental_Winds()
    else:
        print('Environmental Winds Found...')
    return None


###### MAIN GOES HERE ######
if __name__ == '__main__':
    main()
    print('\n****************************\n* ERA5 DOWNLOADS COMPLETE! *\n****************************\n')