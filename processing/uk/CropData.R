library(rgdal)
setwd("~/Downloads/")
fgdb <- "~/Downloads/Download_Crops2017_1399393/lccm-2017_3308728/lccm-2017_3308728.gdb"

# List all feature classes in a file geodatabase
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)
fc <- readOGR(dsn=fgdb,layer="crop_map_2017")
fc@data
mode_func <- function(x, na.rm){
  as.numeric(names(which.max(table(x))))
}

fc@data$cc <- as.numeric(fc@data$crop_code)

library(raster)
r <- raster(extent(c(0, 7e5, 0, 1.25e6)), res= c(1000,1000))
rr <- rasterize(fc, r, field = "cc", fun = mode_func)
writeRaster(rr, "~/Documents/PhD/GIT/UKplantSim/data/Crop2017", format = "GTiff", overwrite = T)
