library(raster)
library(rasterVis)
library(rgdal)
library(gridExtra)
library(sp)
library(latticeExtra)
library(ggplot2)
tpi = raster("CF_TPI_full.tif")
ele = raster("CF_elevation_full.tif")

tpi_norm = raster(crs = crs(tpi), ext = extent(ele), 
                  res = c(10, 10))
tpi_corrected = resample(tpi, tpi_norm)
plot(tpi_corrected)
writeRaster(tpi_corrected, "CF_TPI_corrected.tif", overwrite=T)
ele = resample(ele, tpi_corrected)
writeRaster(ele, "CF_elevation_full.tif", overwrite = T)
ele_crop <- raster("CF_elevation.tif")
ele_crop = resample(ele_crop, tpi_corrected)
writeRaster(ele_crop, "CF_elevation.tif", overwrite = T)

lcm <- raster("data/LCM.tif")
colours = c("red", "darkgreen", "brown", "green", "lightgreen", "green3", "yellow4", "yellow",
            "purple", "pink", "turquoise3", "lavender", "navy", "blue", "gold", "gold", "lightyellow","lightyellow", "lightblue", "black", "grey")
labels = c("Broadleaved woodland", "Coniferous woodland", "Arable",
"Improved grassland", "Neutral grassland", "Calcareous grassland", "Acid grassland",
"Fen, Marsh, Swamp", "Heather",
"Heather grassland", "Bog",
"Inland Rock", "Saltwater",
"Freshwater", "Supra-littoral rock",
"Supra-littoral sediment", "Littoral rock",
"Littoral sediment", "Saltmarsh", "Urban",
"Suburban")
lcm_crop <- crop(lcm, tpi_norm)
lcm_crop2 <- mask(lcm_crop, cf)
# get rid of inconsistencies inside CF boundary
weird_cells <- which(lcm_crop2@data@values > 11)
lcm_crop[weird_cells] = 8
plot(lcm_crop, col = colours, legend = T, legend.lab = labels)
writeRaster(lcm_crop, "LCM.tif", overwrite = T)

cf <- readOGR("CorsFochno.shp")
lcm_crop2[weird_cells] = 8
plot(lcm_crop2, col = colours[1:11], legend = F)
plot(cf, add = T)
cf <- rasterize(cf, tpi_norm)
plot(cf)
writeRaster(cf, "CF_outline.tif", overwrite = T)

lcm <- readAll(raster("LCM.tif"))
lcm <- crop(lcm, cf)
lcm <- ratify(lcm)
rat <- levels(lcm)[[1]]
rat$landcover <- labels[rat[[1]]]
levels(lcm) <- rat
ditches <- readAll(raster("Ditches_full.tif"))
ditches <- ratify(ditches)
rat = levels(ditches)[[1]]
rat$ditch = c("Primary Ditch", "Secondary Ditch")
levels(ditches) <- rat
pdf("../plots/LCM.pdf", paper = "a4", height = 11, width = 8)
l1 <- levelplot(lcm, col.regions = colours[rat$ID], main = list("A", x = 0.1))
l2 <- levelplot(ditches, col.regions = c("blue", "purple"), main = list("B", x = 0.1), maxpixels = 1e6, par.settings=list(panel.background=list(col="lightgreen")))
grid.arrange(l1, l2, nrow=2)
dev.off()


library(rgrass)
