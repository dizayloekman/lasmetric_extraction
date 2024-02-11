############
## UAV-LiDAR Wall-to-wall Metrics Extraction
############

##Input data
## 1. Normalized lasfile data (.laz) (file naming of lasfile is following the petak name, e.g. TPGE010701 (petak keyid in Riau))
## 2. AOI/petak shapefile

##Output data
## 1. folder to store the output by petak (e.g. TPGE010701_grid_metric)
## 2. .csv file of lasmetrics by grid (e.g. TPGE010701_metric.csv)
## 3. .tif file of lasmetrics by grid (e.g. TPGE010701_metric.tif)

##Packages required
##!!install all required packages if not available yet!!
##!!lidRmetrics .tar.gz file in the directory folder (or any folder) (lidRmetrics is not available in R package repository)

##install lidRmetrics from .tar.gz file source
lidRmetrics <- "C:/Users/00026152/OneDrive - PT. Global Komunikasi Mandiri/DATA/04_Etc/R/ptompalski-lidRmetrics-6ba3b13.tar.gz" #lidRmetrics .tar.gz file path
install.packages(lidRmetrics, repos=NULL, type="source")

devtools::install_github("ptompalski/lidRmetrics")
library(lidRmetrics)

##install required packages from CRAN
install.packages("raster","rgdal","rgeos","reshape2","tidyverse","readr","ggstatsplot","moments","foreign","dplyr","plyr","ggplot2","readxl","tidyr","xlsx","glmnet", "mlbench", "kernlab", "mlbench","RWeka", "party", "gbm", "Cubist", "earth", "caret","readxl","writexl","randomForest", "gstat", "mlbench","cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "rnaturalearth", "rnaturalearthdata","GGally","rpart","ipred","nnet","rlang","varImp", "caret","randomForest","geometry","quadmesh","Rvcg","sp", "lidR", "sf", "Lmoments", "lidRmetrics", "terra")

##Call all required packages
Packages <- c("raster","rgdal","rgeos","reshape2","tidyverse","readr","ggstatsplot","moments","foreign","dplyr","plyr","ggplot2","readxl","tidyr","xlsx","glmnet", "mlbench", "kernlab", "mlbench","RWeka", "party", "gbm", "Cubist", "earth", "caret","readxl","writexl","randomForest", "gstat", "mlbench","cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "rnaturalearth", "rnaturalearthdata","GGally","rpart","ipred","nnet","rlang","varImp", "caret","randomForest","geometry","quadmesh","Rvcg","sp", "lidR", "sf", "Lmoments", "lidRmetrics", "terra")
lapply(Packages, require, character.only = TRUE)

##Setting up the working directory
setwd("C:/Users/00026152/OneDrive - PT. Global Komunikasi Mandiri/Documents/temp/20230615/New folder/") #adjust the working directory according to where the file is stored

##Get all file names with .las or .laz extension from wd 
data_frame_names <- list.files(pattern = "*.laz") #adjust file extention to .laz or .las

#Directory to save data
#dir<-"//172.30.178.92/argis/prod/UAV_LiDAR/Input/20230612/EPEL_36_RI/" 
dir<- "C:/Users/00026152/OneDrive - PT. Global Komunikasi Mandiri/Documents/temp/lasmetric/20230615/EPEL_36_RI"

##Landuse or AOI shapefile
#landuse <- st_read("//172.30.178.92/argis/prod/UAV_LiDAR/Pilot/05 Shapefile/AoI/AOI_LiDAR_2023_47N.shp") #adjust AOI shapefile location
landuse <- st_read("C:/Users/00026152/OneDrive - PT. Global Komunikasi Mandiri/DATA/09_Pilot_Digital_Forest_Inventory/Data/AOI_LiDAR_2023_Revised/AOI_LiDAR_2023_47N.shp")

##set of metric fuction
##need to run metrics_echo due to some issue in lidRmetric package
metrics_echo <- function(ReturnNumber, NumberOfReturns, z=NULL, zmin=NA) {
  
  if (!is.na(zmin) & is.null(z)) {warning("Both z and zmin parameters are required to apply zmin filter. zmin threshold not applied.")}
  
  if (!is.na(zmin) & !is.null(z))  {
    zfilter <- z>zmin
    ReturnNumber <- ReturnNumber[zfilter]
    NumberOfReturns <- NumberOfReturns[zfilter]
  }
  
  
  filter <- ReturnNumber > 0 & NumberOfReturns > 0 & ReturnNumber <= NumberOfReturns
  ReturnNumber <- ReturnNumber[filter]
  NumberOfReturns <- NumberOfReturns[filter]
  
  n = length(ReturnNumber)
  
  first = sum(ReturnNumber == 1)/n * 100
  intermediate = sum(ReturnNumber > 1 & ReturnNumber < NumberOfReturns)/n * 100
  last = sum(ReturnNumber == NumberOfReturns & ReturnNumber > 1)/n * 100
  
  
  #single/multiple
  single <- sum(NumberOfReturns==1) / n * 100
  multiple <- sum(NumberOfReturns > 1) / n * 100
  
  
  out <- list(pfirst = first, pintermidiate = intermediate, plast = last, psingle=single, pmultiple=multiple)
  
  return(out)
}


metrics_set4  <- function(x, y, z, i,
                          ReturnNumber,NumberOfReturns,
                          zmin=NA, 
                          threshold = c(2,5), 
                          dz=1, 
                          interval_count=10, 
                          zintervals=c(0, 0.15, 2, 5, 10, 20, 30),
                          pixel_size=1,
                          vox_size=1) {
  
  m_set1    <- metrics_set1(z = z, zmin = zmin, threshold = threshold, dz = dz, interval_count = interval_count, zintervals = zintervals)
  m_rumple  <- metrics_rumple(x = x, y = y, z = z, pixel_size = pixel_size)
  m_vox     <- metrics_voxels(x = x, y = y, z = z, vox_size = vox_size, zmin = zmin)
  m_kde     <- metrics_kde(z = z, zmin = zmin)
  m_echo    <- metrics_echo(z = z, ReturnNumber = ReturnNumber, NumberOfReturns=NumberOfReturns)
  m_HOME    <- metrics_HOME(z = z, i = i, zmin = zmin)
  m_shape <- lidR::stdshapemetrics(x, y, z)
  m_tree <- lidR::stdtreemetrics(x, y, z)
  
  m <- c(m_set1, m_rumple, m_vox, m_kde, m_echo, m_HOME, m_shape, m_tree)
  
  return(m)
  
}

##Loop for deriving lasmetrics
##call each data to process
system.time(for (i in data_frame_names) {
  #try to finish all loop ignoring the error
  tryCatch({
    #Read data without scan angle
    print(paste("Read Lasfile Started", sep=" "))
    print(Sys.time())
    dat <- readALSLAS(i, select = "* -a", filter = "-drop_z_below 0 -drop_z_above 40") #read lasfile
    print(paste("Lasfile Loaded", sep=" "))
    print(Sys.time())
    
    #clip lasfile by petak boundary
    ##call lasfile by petak name
    fld<- basename(data_frame_names)
    name <- tools::file_path_sans_ext(fld) #call each lasfile (by petak)
    
    ##Get list of petak within flight mission
    len<-lengths(strsplit(name, '_')) #to identify/count petak name from lazfile name
    ptk_list<-list()
    j=1
    while (j<=len){
      ptk_list[j]<-scan(text=name, sep="_", what="", quiet=TRUE)[j]
      j <- j + 1
    }
    
    ##Clipping Process
    for(k in ptk_list){
      newdir <- paste0(dir, k) #set name of folder
      path <- dir.create(newdir) #create new folder to save lasmetric output each petak
      petak<-subset(landuse,landuse$Keyid==k) #select each petak #replace id with petak name
      dat_clip <- clip_roi(dat, petak) #clip las by petak
      print(paste(k, "Clipped", sep=" "))
      lasname<-paste0(newdir, "/", k,".laz")
      writeLAS(dat_clip,lasname) #export las by petak
      
      #Read clipped lazfile
      clippedlas <- readALSLAS(lasname)
      
      #Extract lasmetrics 
      print(paste("Extract Metric Started", sep=" "))
      print(Sys.time())
      m1 <- pixel_metrics(clippedlas, ~metrics_set4(X, Y, Z, Intensity, ReturnNumber, NumberOfReturns, zmin=NA, threshold = c(2,5), dz=1, interval_count=10, zintervals=c(0, 0.15, 2, 5, 10, 20, 30), pixel_size=1, vox_size=1), 15, na.rm=TRUE) #extract user-defined lasmetric parameter by pixel_metrics function #output=raster
      metric <- as.data.frame(m1) #get raster value as dataframe
      #metric <- as.tibble(m1) #use this insted if error occured
      #metric <- as_tibble(m1) #use this insted if error occured
      xy <- coordinates(raster(m1)) #generate centroid xy coordinate each grid/cell
      #xy <- xyFromCell(m1, 1:ncell(m1)) #use this insted if error occured
      metric$GridId <- seq.int(nrow(metric)) #set grid id
      metric$KeyId <- k #set keyid
      metric <- merge(metric, xy, by = "row.names") #merge lasmetric & xy coordinate
      metric <- metric[c(109, 110, 2:108, 111:112)] #arrange the column
      file.name = paste0(k,"_metric",".csv") #set the name of .csv file
      #write.csv(metric,paste(path, file.name, sep = ""), row.names = FALSE, quote = FALSE)
      m2 <- write.table(metric, paste0(newdir, "/", file.name, sep = ""), row.names = FALSE, quote = FALSE, sep="|") #write .csv file of lasmetric
      #file.name2 = paste0(k,"_metric",".tif") #set name of .tif file
      #m3 <- writeRaster(m1, paste0(path, "/", file.name2, sep = ""), overwrite=TRUE) #write .tif file of lasmetric
      print(paste(name, "Lasmetric Extracted", sep=" "))
      print(Sys.time())
      #rm(dat_clip) #remove clipped lasfile from memory
      #gc() #clear memory
      print(paste(name, "Clipped Lasfile Removed from Memory", sep=" "))
      print(paste(name, "Process Finished", sep=" "))
      print(Sys.time())
    }
    }, 
  #Print error if any 
  error=function(e){cat(i,conditionMessage(e), "/n")})
})
