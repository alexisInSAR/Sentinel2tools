#! /usr/bin/env python3
# -*- coding: iso-8859-1 -*-

##########################################################################################
# Header information 

"""S2_to_geotif.py: Scripts to convert the Sentinel-2 images into .tif images."""

__author__ = "Alexis Hrysiewicz"
__copyright__ = "Copyright 2022"
__credits__ = ["Alexis Hrysiewicz"]
__license__ = "GPLV3"
__version__ = "1.0.0"
__maintainer__ = "Alexis Hrysiewicz"
__status__ = "Release 1.0 Beta"
__date__ = "Jun. 2022"

# Update
# 03/09/2021: - major error with the RGB images
#             - removal of tmp.tif
#             - error during the cloud mask due to a false option.

##########################################################################################
# Python packages
import os
import os.path
import optparse
import sys
import shutil
import rasterio
import fiona
import rasterio.mask
import numpy as np
from skimage import exposure
from rasterio.warp import reproject, Resampling
from zipfile import ZipFile
from rasterio.enums import ColorInterp

###########################################################################
# Class definition for the user options 

class OptionParser (optparse.OptionParser):

    def check_required(self, opt):
        option = self.get_option(opt)
        if getattr(self.values, option.dest) is None:
            self.error("%s option not supplied" % option)

###########################################################################
# Definition of options 

if len(sys.argv) < 5:
    prog = os.path.basename(sys.argv[0])
    print("example: python3 %s -f S2.list -o ./Save -m crop.shp -c y -b 'B01,B02,B03' -r 'RGB,IR' -i 'NDVI,NDWI'" %
          sys.argv[0])
    sys.exit(-1)
else:
    usage = "usage: %prog [options] "
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--file", dest="file", action="store", type="string", default='S2.list')                    #The list of Sentinel-2 images
    parser.add_option("-o", "--output", dest="output", action="store", type="string", default='./Save')                 #The output directory
    parser.add_option("-m", "--crop", dest="crop", action="store", type="string", default='')                           #The cropping option: give a .shp file
    parser.add_option("-c", "--cloud", dest="cloud", action="store", type="string", default='n')                        #The mask of clouds from ESA processing: y or n
    parser.add_option("-b", "--bands", dest="bands", action="store", type="string", default='all')                      #The saved band: all or 'B01,B02,B04' (a list of the selected bands) or ''
    parser.add_option("-r", "--rgb", dest="rgb", action="store", type="string", default='')                             #The saved false-color images: RGB, IR or 'RGB,IR' or ''
    parser.add_option("-i", "--index", dest="index", action="store", type="string", default='')                         #The saved index images: NDVI, NDWI or 'NDVI,NDWI' or ''
    (options, args) = parser.parse_args()

###########################################################################
# Definition of functions

#Function to compute the index image
def computation_index(path_image,index):
    if index == 'NDVI':
        #Selection of correct bands
        name_bi1 = path_image+'_B04.tif'
        name_bi2 = path_image+'_B08.tif'

        #Read the bands
        bi1 = rasterio.open(name_bi1).read(1).astype(np.float32) 
        bi1_meta = rasterio.open(name_bi1).meta
        bi2 = rasterio.open(name_bi2).read(1).astype(np.float32) 

        #Computation 
        rastermatrix = (bi2 - bi1) / (bi2 + bi1)

        #Update the metadata
        out_meta = bi1_meta
        out_meta.update({"driver": "GTiff","dtype":'float32','nodata': np.nan})

    elif index == 'NDWI':
        #Selection of correct bands
        name_bi1 = path_image+'_B08.tif'
        name_bi2 = path_image+'_B12.tif'
        #Read the bands
        bi1 = rasterio.open(name_bi1)
        bi2 = rasterio.open(name_bi2)
        bi2_meta = rasterio.open(name_bi1).meta

        #Interpolation
        bi1res = np.zeros(bi2.read(1).astype(np.float32).shape)
        reproject(bi1.read(1).astype(np.float32), bi1res, src_transform=bi1.transform, src_crs=bi1.crs, dst_transform=bi2.transform, dst_crs=bi2.crs, resampling=Resampling.nearest)

        #Computation 
        rastermatrix = (1*bi1res.astype(np.float32) - 1*bi2.read(1).astype(np.float32)) / (1*bi1res.astype(np.float32) + 1*bi2.read(1).astype(np.float32))
        
        #Update the metadata
        out_meta = bi2_meta
        out_meta.update({"driver": "GTiff","dtype":'float32','nodata': np.nan})

    return rastermatrix, out_meta

###########################################################################
###########################################################################
# MAIN
###########################################################################
###########################################################################

###########################################################################
#Check the user parameters
print('The file of S2 is %s' % options.file)
if not os.path.isfile(options.file):
    print("ERROR: Please enter a valid name")
    sys.exit()

print('The output directy is %s' % options.output)
if not os.path.isdir(options.output):
    print("The output directory does not exist, it will be created.")

if options.crop == '':
    print('No cropping')
else:
    if not os.path.isfile(options.crop):
        print("ERROR: Please enter a valid name")
        sys.exit()
    else:
        print('The .shp used for cropping is %s' % options.crop)

if options.cloud == 'y':
    print('Masking of clouds: YES')
elif options.cloud == 'n':
    print('Masking of clouds: NO')
else:
    print('ERROR: Bad option for cloud masking')
    sys.exit()

if not options.bands == '':
    if options.bands == 'all':
        print('\tThe all bands will be processed')
        listbandproc = ['B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','TCI']
    else:
        listbandproc = options.bands.split(',')
        for bi in listbandproc:
            if bi == 'B01' or bi == 'B02' or bi == 'B03' or bi == 'B04' or bi == 'B05' or bi == 'B06' or bi == 'B07' or bi == 'B08' or bi == 'B09' or bi == 'B10' or bi == 'B11' or bi == 'B12' or bi == 'TCI':
                print('\tThe band %s will be processed' % bi)
            else:
                print("ERROR: Please enter a valid band")
                sys.exit()
else:
    listbandproc = ''

if not options.rgb == '':
    listrgb = options.rgb.split(',')
    for rgbi in listrgb:
        if rgbi == 'RGB' or rgbi == 'IR':
            print('\tThe %s image will be saved' % rgbi)
        else:
            print("ERROR: Please enter a valid RGB displaying")
            sys.exit()
else:
    listrgb = ''

if not options.index == '':
    listindex = options.index.split(',')
    for index in listindex:
        if index == 'NDVI' or index == 'NDWI':
            print('\tThe %s image will be saved' % index)
        else:
            print("ERROR: Please enter a valid index")
            sys.exit()
else:
    listindex = ''

###########################################################################
# Main processing
nbi = 1      
nb_images = 1          
fid = open(options.file,'r')
for file in fid:
    nb_images = nb_images + 1
fid.close() 
nb_images = nb_images -1                                                                                                    #To have a progress bar...
                                                                                    
fid = open(options.file,'r')
for file in fid:
    file = file.split('\n')[0]                                                                                      #We need to remove the carriage return
    ###########################################################################
    print('Processing for the %s file...' % (file))

    ###########################################################################
    #Unzip the images
    unzip_step = False                                                                                              #Here, we save if we need to unzip the file
    if '.zip' in file:
        print('\tWe have a .zip file: we need to unzip it. But it is possible that the file already exists')
        #Check the unzipped file
        name=file.split('.')[0]                                                                    
        if os.path.exists(name+'.SAFE'):
            print('\tThe file is already unzipped.')
        else:
            print('\tThe file is not unzipped.')
            unzip_step = True 
            #Unzipping using a Python package
            path=file.split('/')[0] 
            if path == file:
                path = '.'
            with ZipFile(name+'.zip', 'r') as zipObj:
                zipObj.extractall(path)   

    ###########################################################################
    #Path of band
    pathim = file.split('.')[0]+'.SAFE'
    listdir0 = [f for f in os.listdir(pathim) if os.path.dirname(os.path.join(pathim, f))]                                  #We explore the directory of the image to find the bands
    for f in listdir0:
        if 'GRANULE' in f:
            pathim = pathim+'/GRANULE'
            listdir1 = [f for f in os.listdir(pathim) if os.path.dirname(os.path.join(pathim, f))]
            for fbis in listdir1:
                if 'L1C' in fbis:
                    pathim = pathim+'/'+fbis
                    listdir2 = [f for f in os.listdir(pathim) if os.path.dirname(os.path.join(pathim, f))]
                    for fter in listdir2:
                        if 'IMG_DATA' in fter:
                            pathim = pathim+'/'+fter
    listband= [f for f in os.listdir(pathim) if os.path.dirname(os.path.join(pathim, f))]  

    ###########################################################################
    #Path of cloud mask
    pathmask = file.split('.')[0]+'.SAFE'
    listdir0 = [f for f in os.listdir(pathmask) if os.path.dirname(os.path.join(pathmask, f))]                              #We explore the directory of the image to find the mask of cloud 
    for f in listdir0:
        if 'GRANULE' in f:
            pathmask = pathmask+'/GRANULE'
            listdir1 = [f for f in os.listdir(pathmask) if os.path.dirname(os.path.join(pathmask, f))]
            for fbis in listdir1:
                if 'L1C' in fbis:
                    pathmask = pathmask+'/'+fbis
                    listdir2 = [f for f in os.listdir(pathmask) if os.path.dirname(os.path.join(pathmask, f))]
                    for fter in listdir2:
                        if 'QI_DATA' in fter:
                            pathmask = pathmask+'/'+fter
    pathmask = pathmask+'/MSK_CLOUDS_B00.gml'

    ###########################################################################
    #Creation of directory
    if not os.path.isdir(options.output):
        os.makedirs(options.output)                                                                                         #Here, we create the output directory

    name_dir = listband[0].split('.')[0].split('_')[1]
    if not os.path.isdir(options.output+'/'+name_dir):                                  
        os.makedirs(options.output+'/'+name_dir)                                                                            #Here, we create the image directory

    ###########################################################################
    #Save the images following the conditions from the user
    if not options.crop == '':
        with fiona.open(options.crop, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]

        for band in listband:
            print('\tProcessing for the %s band...' % (band))
            with rasterio.open(pathim+'/'+band) as src:
                out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
                out_meta = src.meta
                out_meta.update({"driver": "GTiff","height": out_image.shape[1],"width": out_image.shape[2],"transform": out_transform})
                
            if options.cloud == 'y':
                name_save = options.output+'/'+name_dir+'/'+'tmp.tif'
                name_save_bis = options.output+'/'+name_dir+'/'+band.split('.')[0]+'.tif'
                with rasterio.open(name_save, "w", **out_meta) as dest:     
                    dest.write(out_image)

                cmd = "ogr2ogr -f 'GeoJSON' %s %s -overwrite" %(options.output+'/'+name_dir+'/mask_cloud.input',pathmask)
                os.system(cmd)
                cmd = 'rio mask %s %s --invert --geojson-mask %s --overwrite' %(name_save,name_save_bis,options.output+'/'+name_dir+'/mask_cloud.input')
                os.system(cmd)

                os.remove(options.output+'/'+name_dir+'/mask_cloud.input')
                os.remove(name_save)

            else:
                name_save = options.output+'/'+name_dir+'/'+band.split('.')[0]+'.tif'
                with rasterio.open(name_save, "w", **out_meta) as dest:     
                    dest.write(out_image)
    else:
        for band in listband:     
            print('\tProcessing for the %s band...' % (band))           
            if options.cloud == 'y':
                name_save = options.output+'/'+name_dir+'/'+'tmp.tif'
                name_save_bis = options.output+'/'+name_dir+'/'+band.split('.')[0]+'.tif'

                cmd = 'rio convert %s %s -f GTiff --overwrite' % (pathim+'/'+band,name_save)
                os.system(cmd)
                cmd = "ogr2ogr -f 'GeoJSON' %s %s -overwrite" %(options.output+'/'+name_dir+'/mask_cloud.input',pathmask)
                os.system(cmd)
                cmd = 'rio mask %s %s --invert --geojson-mask %s --overwrite' %(name_save,name_save_bis,options.output+'/'+name_dir+'/mask_cloud.input')
                os.system(cmd)

                os.remove(options.output+'/'+name_dir+'/mask_cloud.input')
                os.remove(name_save)
            else:
                name_save = options.output+'/'+name_dir+'/'+band.split('.')[0]+'.tif'
                cmd = 'rio convert %s %s -f GTiff --overwrite' % (pathim+'/'+band,name_save)
                os.system(cmd)

    ###########################################################################
    # Save the full path of band without number
    name_input = listband[0].split('.')[0].split('_')[0]+'_'+listband[0].split('.')[0].split('_')[1]
    
    ###########################################################################
    #Save the false color images 
    if not listrgb == '':
        for rgbi in listrgb:
            print('\tProcessing for the %s false-color image...' % (rgbi))

            #Reading of the bands according to the false-colour image
            if rgbi == 'RGB':
                name_save = options.output+'/'+name_dir+'/'+listband[0].split('.')[0].split('_')[0]+'_'+listband[0].split('.')[0].split('_')[1]+'_RGB.tif'
                bandred = rasterio.open(options.output+'/'+name_dir+'/'+name_input+'_B04.tif')
                bandgreen = rasterio.open(options.output+'/'+name_dir+'/'+name_input+'_B03.tif')
                bandblue = rasterio.open(options.output+'/'+name_dir+'/'+name_input+'_B02.tif')
            elif rgbi == 'IR':
                name_save = options.output+'/'+name_dir+'/'+listband[0].split('.')[0].split('_')[0]+'_'+listband[0].split('.')[0].split('_')[1]+'_IR.tif'
                bandred = rasterio.open(options.output+'/'+name_dir+'/'+name_input+'_B08.tif')
                bandgreen = rasterio.open(options.output+'/'+name_dir+'/'+name_input+'_B04.tif')
                bandblue = rasterio.open(options.output+'/'+name_dir+'/'+name_input+'_B03.tif')

            #Merging of bands
            RGBimage = np.array([bandred.read(1), bandgreen.read(1), bandblue.read(1)]).transpose(1,2,0)
            p2, p98 = np.percentile(RGBimage, (2,98))
            RGBimage = exposure.rescale_intensity(RGBimage, in_range=(p2, p98)) /10000
            RGBimage = RGBimage.transpose(2,0,1)

            #Saving of images
            out_meta = bandred.meta
            out_meta.update({"driver": "GTiff","count": 3, "dtype": rasterio.uint8})
            nameRGB = band.split('.')[0].split('_')[1]
            with rasterio.open(name_save, "w", **out_meta) as dest:     
                dest.colorinterp = [ColorInterp.red, ColorInterp.green, ColorInterp.blue]
                dest.write((RGBimage[0,:,:]/np.max(RGBimage[0,:,:]))*255,1)
                dest.write((RGBimage[1,:,:]/np.max(RGBimage[1,:,:]))*255,2)
                dest.write((RGBimage[2,:,:]/np.max(RGBimage[2,:,:]))*255,3)

    ###########################################################################
    #Save the index images
    if not listindex == '': 
        for index in listindex:
            print('\tProcessing for the %s index image...' % (index))

            rastermatrix, out_meta = computation_index(options.output+'/'+name_dir+'/'+name_input,index)                                            #Computation of indexes using the defined function
            with rasterio.open(options.output+'/'+name_dir+'/'+name_input+'_'+index+'.tif', "w", **out_meta) as dest:                               #Save the index image
                dest.write(rastermatrix, indexes=1)

    ###########################################################################
    #Finalisation
    print('\tRemove the tmp files')
    if os.path.isfile(options.output+'/'+name_dir+'/tmp.tif'):
        os.remove(options.output+'/'+name_dir+'/tmp.tif') 
    if not options.bands == 'all':
        for bi in listband:
            biname = bi.split('.')[0].split('_')[-1]
            if not biname in listbandproc:
                os.remove(options.output+'/'+name_dir+'/'+name_input+'_'+biname+'.tif') 

    if unzip_step:
        shutil.rmtree(file.split('.')[0]+'.SAFE')          
    
    print('\tProcessing: %d on %d image(s).' % (nbi,nb_images))
    nbi = nbi + 1

fid.close()