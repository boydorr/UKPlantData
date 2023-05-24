library(rgdal)
library(raster)

mode_func <- function(x, na.rm){
  as.numeric(names(which.max(table(x))))
}

soils <- readOGR("../data/Hutton_Soils250K_OpenData/qmsoils_UCSS_v1.2.shp")
r <- raster(extent(c(0, 7e5, 0, 1.25e6)), res= c(1000,1000))
rr <- rasterize(soils, r, field = "SERCDE1", fun = mode_func)
writeRaster(rr, "~/Documents/PhD/GIT/UKplantSim/data/HuttonSoils", format = "GTiff", overwrite = T)
