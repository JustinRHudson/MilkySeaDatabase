'''
    This script coarsens the World Seaborne Trade Monitoring System (WSTMS)
    data from 500-meter resolution to ~.1 degree resolution, coarsening by a
    factor of 20.
'''

###### IMPORTS GO HERE ######
## Standard Libraries ##
from PIL import Image # type: ignore
import gc
import os
## Public Libraries ##
import numpy as np

###### FUNCTIONS GO HERE ######
def load_tif(tif_path:str) -> Image.Image:
    '''
        Loads in the WSTMS tif file and returns it as an PIL Image

        Parameters:
            tif_path (str): The local filepath to the WSTMS tif file
        
        Returns:
            gsd_image (Image.Image): A PIL Image of global shipping density
    '''
    if os.path.exists(tif_path):
        Image.MAX_IMAGE_PIXELS = None
        gsd_image = Image.open(tif_path)
    else:
        raise ValueError('Provided path not found')
    
    return gsd_image

def return_pixel_coordinates() -> tuple[np.ndarray]:
    '''
        Returns two numpy arrays which contain the coordinates for the pixels
        within the image we're rescaling

        Parameters:
            None
        
        Returns:
            pixel_lats (np.ndarray): The latitudes for each pixel
            pixel_lons (np.ndarray): The longitudes for each pixel
    '''

    pixel_lats = np.arange(-84.49,85.0,0.005)[::-1]
    pixel_lons = np.arange(-180.0128112750,180.0128112750,0.005)

    return pixel_lats,pixel_lons

def split_image(image:Image) -> tuple[Image,Image,float,float]:
    '''
        Divides the image into two halves to reduce the memory requirements for
        processing.

        Parameters:
            image (Image): The image to be split into two halves
        
        Returns:
            left_half (Image): The left half of the image
            right_half (Image): The right half of the image
            im_width (float): The width (in pixels) of the original image
            im_height (float): The height (in pixels) of the original image
    '''

    #get the height and width of the image
    im_width = image.size[0]
    im_height = image.size[1]
    #divide into two halves
    left_half = image.crop((0,0, int(im_width/2), int(im_height)))
    right_half = image.crop((int(im_width/2),0,int(im_width),int(im_height)))
    
    return left_half,right_half,im_width,im_height

def rescale_image_half(image_half:Image,im_width:int,im_height:int) -> np.ndarray:
    '''
        Rescales halves of the image to our desired resolution of ~0.1 degrees.

        Parameters:
            image_half (Image): The half of the image we're rescaling
            im_width (int): The original width of the image
            im_height (int): The original height of the image
        
        Returns:
            rescaled_image (np.ndarray): The rescaled half of the image
    '''

    #convert from Image to numpy array
    image_half_array = np.reshape(np.asarray(image_half),(im_height,int(im_width/2)))
    #make a new array to hold the rescaled array
    rescaled_image = np.empty((int(im_height/20),int(im_width/40)))
    #fill the new array
    for i in range(rescaled_image.shape[0]-1):
        for j in range(rescaled_image.shape[1]-1):
            rescaled_image[i,j] = np.mean(image_half_array[i*20:(i+1)*20,j*20:(j+1)*20])

    return rescaled_image

def rescale_lat_lons(lats:np.ndarray,lons:np.ndarray,im_width:int,im_height:int) -> tuple[np.ndarray]:
    '''
        Rescales the pixel lats/lons to match the resolution of the final image

        Parameters:
            lats (np.ndarray): The latitudes to be rescaled
            lons (np.ndarray): The longitudes to be rescaled
            im_width (int): The original width of the image
            im_height (int): The original height of the image

        Returns:
            rescaled_lats (np.ndarray): The rescaled latitudes
            rescaled_lons (np.ndarray): The rescaled longitudes
    '''

    #make arrays to hold the scaled lats/lons
    rescaled_lats = np.empty((int(im_height/20)))
    rescaled_lons = np.empty((int(im_width/20)))
    #fill the arrays
    #start with lats
    for i in range(len(rescaled_lats) - 1):
        rescaled_lats[i] = np.mean(lats[i*20:(i+1)*20])
    #now lons
    for i in range(len(rescaled_lons) - 1):
        rescaled_lons[i] = np.mean(lons[i*20:(i+1)*20])
    #fill in the final values
    rescaled_lats[-1] = rescaled_lats[-2] - 0.1
    rescaled_lons[-1] = rescaled_lons[-2] + 0.1

    return rescaled_lats,rescaled_lons

def log_scale_image(unlogged_image:np.ndarray) -> np.ndarray:
    '''
        Scales the image using the base 10 log

        Parameters:
            unlogged_image (np.ndarray):  The image to be log scaled
        
        Returns:
            log_scaled_image (np.ndarray): The log-scaled image
    '''

    log_scaled_image = np.log10(unlogged_image)

    return log_scaled_image

def recombine_image_halves(left_half:np.ndarray,right_half:np.ndarray,im_width:int,im_height:int) -> np.ndarray:
    '''
        recombines the halves of the image into a single image.

        Parameters:
            left_half (np.ndarray): The left half of the image
            right_half (np.ndarray): The right half of the image
            im_width (int): The original width of the image
            im_height (int): The original height of the image
        
        Returns:
            recombined_image (np.ndarray): The recombined halves of the image
    '''

    #make an empty array for the image
    recombined_image = np.empty((int(im_height/20),int(im_width/20)))
    #fill it
    recombined_image[:,0:int(im_width/40)] = left_half
    recombined_image[:,int(im_width/40):] = right_half

    return recombined_image

def save_results(save_location:str,log_scaled_image:np.ndarray,image_lats:np.ndarray,image_lons:np.ndarray) -> None:
    '''
        Saves the image and the lat/lons as .npy objects so they can be used to
        make plots with.

        Parameters:
            save_location (str): Where to save the reprocessed image and lat/lons
            log_scaled_image (np.ndarray): The reprocessed imge
            image_lats (np.ndarray): the reprocessed lats
            image_lons (np.ndarray): the reprocessed lons

        Returns:
            None
    '''

    #check if the save location is a valid location
    if not os.path.exists(save_location):
        raise ValueError('save_location does not exist.')
    else:
        if save_location[-1] == '/':
            np.save(save_location + 'gsd_lats',image_lats)
            np.save(save_location + 'gsd_lons',image_lons)
            np.save(save_location + 'gsd_coarse',log_scaled_image)
        else:
            np.save(save_location + '/gsd_lats',image_lats)
            np.save(save_location + '/gsd_lons',image_lons)
            np.save(save_location + '/gsd_coarse',log_scaled_image)

    return None

def reprocess_image_and_save(image_location:str,save_location:str) -> None:
    '''
        Reprocesses the WSTMS image and returns it coarsened by a factor of
        twenty so that it's more usable for global applications.

        Parameters:
            image_location (str): Where the WSTMS image is stored locally
            save_location (str): Where to save the coarsened image/data
        
        Returns:
            None
    '''

    #check if the locations provided are valid paths
    if not os.path.exists(image_location):
        raise ValueError('image_location path does not exist')
    elif not os.path.exists(save_location):
        raise ValueError('save_location path does not exist')

    #load the image
    gsd_image = load_tif(image_location)
    #split the image in half so its easier to work with
    left_half,right_half,im_width,im_height = split_image(gsd_image)
    # remove the image from memory
    del gsd_image
    gc.collect()
    #rescale each half of the image
    left_rescaled = rescale_image_half(left_half,im_width,im_height)
    right_rescaled = rescale_image_half(right_half,im_width,im_height)
    # remove the original left and right halves from memory
    del left_half
    del right_half
    gc.collect()
    #recombine the halves
    recombined_image = recombine_image_halves(left_rescaled,right_rescaled,im_width,im_height)
    #now get appropriate lat/lons for the rescaled image
    pixel_lats,pixel_lons = return_pixel_coordinates()
    rescaled_lats,rescaled_lons = rescale_lat_lons(pixel_lats,pixel_lons,im_width,im_height)
    #log scale the image
    rescaled_image = log_scale_image(recombined_image)
    #save the results of the reprocessing
    save_results(save_location,rescaled_image,rescaled_lats,rescaled_lons)

    return None

###### MAIN GOES HERE ######
if __name__ == '__main__':
    pass