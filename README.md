# From Sailors to Satellites: A Curated Database of Bioluminescent Milky Seas Spanning 1600 - Present
#### By Justin Hudson, and Steven D. Miller
------------

This repository contains the code, data, and environment files used to create the figures seen in our manuscript. This repository also contains the database both as a human readable PDF and a machine readable tab separated values file (.tsv).

The code is split into several different scripts and a Jupyter notebook. To replicate the figures only one script needs to be run (Download_ERA5_Data.py), the rest is handled by the jupyter notebook Figure_Gen.ipynb.

The python script Download_ERA5_Data.py uses the Copernicus Climae Date Store API to download the ERA5 data used in this manuscript and places the folders in the DATA folder. The size of the ERA5 files prevents them from being included in this GitHub repository so a script to download them had been provided instead.

The jupyter notebook Figure_Gen.ipynb processes the data and generates the figures for the paper. To create Figures 9-12, the code in bootstrapping.py must be called, this however can take several hours to run for 5,000 iterations, it is currently set to only run 100 times. Comment out or change that line in Figure_Gen.ipynb in order to run it for all 5,000 iterations.

## Steps to run the code

1. Set up the necessary python environments. There are two .yml files within this repository: cdsapi_env.yml and figgen_env.yml.
    - cdsapi_env.yml sets up the python environment needed to download the ERA5 data using the Copernicus Climate Date Store API
    - figgen_env.yml sets up the python environment needed to run the notebook Figure_Gen.ipynb
    - [Here is a link on how to set up a python environment from a yml file using conda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)
2. Download the necessary data
    - Download the HadISST SST data [here](https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html), download the file named HadISST_sst.nc.gz near the bottom of the page. Unzip the file and place the HadISST_sst.nc file in the DATA folder
    - Download the the WSTMS data [here](https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density), download the Global Ship Density folder, unzip the folder and place the entire folder within the DATA folder
3. Run the Download_ERA5_Data.py file
    - This assumes you have set up the correct API and have an account with Copernicus Climate Data Store
    - [Here are the set-up instructions for the API](https://cds.climate.copernicus.eu/how-to-api)
    - Once you have the API set up run the script inside the cdsapi_env environment, this should download the ERA5 data needed and place it into the DATA folder
4. Run the Figure_Gen.ipynb notebook
    - The kernel/environment used to run this notebook is the figgen_env environment. The notebook should run from start to finish producing the figures as it goes.
    - NOTE: The bootstrapping code takes a long time to run, for demonstration purposes it is set to only use 100 iterations by default, please change that cell to 5,000 iterations which is what was used in the paper if you want to truly recreate the figures.
    - The notebook should display the figures but all figures will be placed into the FIGURES folder

## Corresponding Author
Justin Hudson: justin.hudson@colostate.edu