# Set of Python script to process and manipulate Sentinel-2 images

**Requirements:**
- Python3 
- Linux/Mac (not tested on Windows) 
- GDAL

**Installation of required Python packages**

```bash
pip3 install os optparse sys shutil rasterio fiona numpy scikit-image zipfile 
```
or
```bash
conda install os optparse sys shutil rasterio fiona numpy scikit-image zipfile 
```

***

## *S2_to_geotif.py*

Script to extract Sentinel-2 images and save images in geotiff. Several options for the bands, index selections are provided, with the otpions for masking the clouds. 

Firstly, a list of Sentinel-2 images must be created. For example: 

```bash
ls ../Data_Sentinel2/*.zip > list_S2.txt
```

i.e., 

```txt
../Data_Sentinel2/S2B_MSIL1C_20210317T063509_N0209_R134_T40KCB_20210318T154922.zip
../Data_Sentinel2/S2B_MSIL1C_20210406T063509_N0300_R134_T40KCB_20210406T075858.zip
```

then: 

```bash
S2_to_geotiff.py -f list_S2.txt -o ./Save -m crop.shp -c y -b 'B01,B02,B03' -r 'RGB,IR' -i 'NDVI,NDWI'
```

With: 
- -f: the list of Sentinel-2 images [mandatory]; 
- -o: the directory of outputs [mandatory];
- -m: the .shp file to crop the images [optional]; 
- -c: cloud masking (y or n) [optional]; 
- -b: the selected stored bands [optional]; 
- -r: the selected composite images [optional]; 
- -u: the selected index images [optional]. 

***

## Author

*Alexis Hrysiewicz,* Postdoctoral Researcher, School of Earth Sciences, University College Dublin, iCRAG
