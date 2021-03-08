#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Title: Intial Spatial Analysis
#Coder: C. Nathan Jones (cnjones7@ua.edu)
#Date: 3/7/2021
#Purpose: Examine hydrogeomorphic features across Chickasawhatchee Swamp
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#1.0 Setup workspace -----------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear memory
rm(list=ls(all=TRUE))

#Defin relevant working directories
data_dir<-"data/"
workspace_dir<-"data/workspace/"

#Download packages 
library(stars)
library(whitebox)
library(igraph)
library(sf)
library(raster)
library(tidyverse)
library(parallel)

#define master proj
p<-"+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#Download relevant data 
dem<-raster(paste0(data_dir,"dem2.tif"))
streams<-st_read(paste0(data_dir, "flowlinest.shp"))

#Define watershed area threshold
threshold=107639  #1 ha

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2.0 Create stream layer -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Export DEM to workspace
writeRaster(dem, paste0(workspace_dir,"dem.tif"), overwrite=T)

#fill single cell depressions
wbt_fill_single_cell_pits(
  dem = "dem.tif",
  output  = "dem_pits.tif",
  wd = workspace_dir
)

#Smooth DEM
wbt_gaussian_filter(
  input = "dem_pits.tif", 
  output = "dem_smoothed.tif",
  wd = workspace_dir)

#breach depressions
wbt_breach_depressions(
  dem =    "dem_smoothed.tif",
  output = "dem_breached.tif",
  fill_pits = F,
  wd = workspace_dir)

#Flow direction raster
wbt_d8_pointer(
  dem= "dem_breached.tif",
  output ="fdr.tif",
  wd = workspace_dir
)

#Flow accumulation raster
wbt_d8_flow_accumulation(
  input = "dem_breached.tif",
  output = "fac.tif",
  wd = workspace_dir
)

#Create Stream Raster
wbt_extract_streams(
  flow_accum = "fac.tif",
  output = "stream.tif",
  threshold = threshold,
  wd = workspace_dir
)

#Convert stream to vector
wbt_raster_streams_to_vector(
  streams = "stream.tif",
  d8_pntr = "fdr.tif",
  output = "streams.shp",
  wd = workspace_dir)

#Read streams layer in 
streams<-st_read(paste0(workspace_dir,"streams.shp"), crs=st_crs(dem@crs))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.0 Define Depressions -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#3.1 Break raster up into meaninful subbasins for wetland delineation-----------
#Define watershed area threshold
threshold_lr=10763900  #1000 ha

#Create Stream Raster
wbt_extract_streams(
  flow_accum = "fac.tif",
  output = "stream_lr.tif",
  threshold = threshold_lr,
  wd = workspace_dir
)

#Break dem into subshed
wbt_subbasins(
  d8_pntr = "fdr.tif",
  streams = 'stream_lr.tif',
  output = "subbasin.tif",
  wd = workspace_dir,
)

#Read subsheds into R and develop unique id
sheds<-raster(paste0(workspace_dir,"subbasin.tif"))
sheds_id<-raster::unique(sheds)

#3.2 Create function to estimate depressions for each subshed ------------------
#Create function
fun<-function(n){

  #Define subbasin
  wbt_equal_to(
    input1 = "subbasin.tif",
    input2 = paste0(sheds_id[n]),
    output = "temp_shed.tif", 
    wd= workspace_dir
  )
  
  #Set 0 to NA
  wbt_set_nodata_value(
    input = "temp_shed.tif", 
    output = "crop.tif",
    back_value = "0", 
    wd= workspace_dir
  )
  
  #Create vector
  wbt_raster_to_vector_polygons(
    input  = "crop.tif",
    output = "crop.shp",
    wd     = workspace_dir
  )
  
  #Crop dem
  dem_temp<-raster(paste0(workspace_dir,"dem_smoothed.tif"))
  crop<-st_read(paste0(workspace_dir,"crop.shp"))
  dem_temp<-crop(dem_temp, crop)
  dem_temp<-mask(dem_temp, crop)
  writeRaster(dem_temp, paste0(workspace_dir,"dem_cropped.tif"), overwrite=T)
  
  #Define depressions
  set.seed(100)
  wbt_stochastic_depression_analysis(dem = "dem_cropped.tif", 
                                     output = "giws.tif", 
                                     rmse = 0.18, 
                                     range = 6, 
                                     iterations = 100, 
                                     wd = workspace_dir)
  
  #3.2 Define depression based on threshold of occurence in stochastic procedure----
  #Reclass raster (any depression that was delineated less than 80% of the time is out!)
  wbt_reclass(
    input = paste0(workspace_dir,"giws.tif"), 
    output = paste0(workspace_dir,"reclass.tif"), 
    reclass_vals = "'0;0;0.8'")
  wbt_reclass(
    input = paste0(workspace_dir,"reclass.tif"), 
    output = paste0(workspace_dir,"reclass.tif"), 
    reclass_vals = "'1;0.8;1'")
  wbt_majority_filter(
    input = "reclass.tif",
    output = "reclass_filter.tif",
    filterx = 10,
    filtery = 10, 
    wd = workspace_dir
  )
  
  #Convert 0 to NA
  giws<-raster(paste0(workspace_dir,"reclass_filter.tif")) 
  giws<-raster::clump(giws)
  
  #Export as polygon
  giws[giws==1]<-NA
  giws<- giws %>% st_as_stars(.) %>% st_as_sf(., merge = TRUE)
  
  #Write polygon shapes to workspace
  st_write(giws, 
           paste0(workspace_dir, "giws.shp"), 
           delete_layer=TRUE)
  
  #Estimate perimeter and area
  wbt_polygon_perimeter(
    input="giws.shp", 
    wd = workspace_dir)
  wbt_polygon_area(
    input="giws.shp", 
    wd = workspace_dir)
  
  #Filter by area and P:A Ratio
  giws<-st_read(paste0(workspace_dir, "giws.shp"))
  giws<-giws %>%
    #Remove small depressions
    filter(AREA>250) %>%
    #Remove oddly large depressions
    filter(AREA<1e5) %>%
    #Remove ditched depressions
    mutate(p_a_ratio = AREA/PERIMETER) %>%
    filter(p_a_ratio>2)
  
  #export giws shape
  giws
}

#3.3 Apply fun -----------------------------------------------------------------
#Create error catching function
execute<-function(m){
  tryCatch(
    fun(m), 
    error = function(e) NA)
}

#Apply function
t0<-Sys.time()
x<-lapply(
  X = seq(1, length(sheds_id)), 
  FUN = execute)
tf<-Sys.time()
tf-t0

#pull results
df<-bind_rows(x)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4.0 Export Map! -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
