using RCall
R"
library(rgdal)
library(raster)

ukpop = raster('data/raw/UK_density_2011.asc')
uk = raster('data/raw/UK.tif')
ukpop_new = extend(ukpop, uk)
writeRaster(ukpop_new, 'data/final/UK_density.tif', format = 'GTiff')
"