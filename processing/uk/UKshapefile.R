r = raster("~/Documents/PhD/Data/CEH/CEH_landcover_2015.tif")
r2 = raster("Documents/PhD/Data/CEH/CEH_landcover_2015_ireland.tif")

r = projectRaster(r, crs = crs(r2), res = c(1000, 1000))
r = extend(r, extent)
r2 = extend(r2, extent)
plot(r2)
plot(r, add = T)
merge(r, r2)

r <- crop(r, c(0, 7e+05, 0, 1250000))
r@data@values[r@data@values > 0] = 250
writeRaster(r, "~/Documents/COVID/UK.tif", format = 'GTiff', overwrite = T)

r = raster("~/Documents/COVID/scotland_run/test_scotland_Hosp30.tif")
plot(r)
mask = raster("~/Documents/PhD/GIT/EcoSISTEM/test/examples/ScotlandDensity2011.tif")
mask = crop(mask, extent(r))
rnew = mask(r, mask)
cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))
plot(rnew, col = c("lightgrey", cols(20)))

library(rasterVis)
library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))
levelplot(rnew,
          layout=c(1, 1),
          zlim = c(0, 200),
          col.regions=c("lightgrey", cols(20)),
          scales=list(draw=FALSE ))
gplot(rnew)+ geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradient(low = 'white', high = 'blue') +
  coord_equal()


all_Hosp <- list.files("~/Documents/COVID/scotland_run/", full.names = TRUE, pattern = "Hosp")
replace = sprintf("%02d", as.numeric(str_extract(all_Hosp, "[0-9].")))
new_all_Hosp <- str_replace(all_Hosp, "[0-9].", replace)
order = order(new_all_Hosp)
# Create a time series raster stack
Hosp_stack <- stack(all_Hosp[order][seq(0, 30, 5)])
Hosp_stack = mask(Hosp_stack, mask)
rasterNames  <- paste("Day", str_extract(names(Hosp_stack), "[0-9].*"))
library(rasterVis)
library(RColorBrewer)
cols <- colorRampPalette(brewer.pal(9,"YlOrRd"))
levelplot(Hosp_stack,
          layout=c(3, 2),
          zlim = c(0, 200),
          names.attr=rasterNames,
          col.regions=c("lightgrey", cols(20)),
          scales=list(draw=FALSE ))

library(raster)
library(rgdal)
library(ncdf4)

ncfname = "~/Documents/PhD/GIT/EcoSISTEM/test/examples/temp/ukv_daily_t1o5m_mean_20200101.nc"
temp = nc_open(ncfname)
tmpin <- raster(ncfname)
crs(tmpin) = '+proj=ob_tran +o_proj=longlat +lon_0=357.5 +o_lon_p=0 +o_lat_p=37.5 +a=6371229 +b=6371229 +to_meter=0.0174532925199 +wktext'
datazones = readOGR("data-raw/datazone_shapefile/SG_DataZone_Bdry_2011.shp")
dz_r = raster(extent(datazones), res = c(1000, 1000))
dz_r = rasterize(datazones, dz_r)
plot(datazones)
plot(dz_r)
reproj1 = projectRaster(from = dz_r, crs = crs(tmpin))
plot(reproj1)
extent(tmpin)[1] = extent(tmpin)[1] - 360
extent(tmpin)[2] = extent(tmpin)[2] - 360
new_tmpin = crop(tmpin, reproj1)
reproj2 = resample(reproj1, new_tmpin)
plot(mask(new_tmpin, reproj2))



reproj = projectRaster(from = tmpin, crs = crs(dz_r))
resampled = resample(reproj, dz_r)
plot(resample(reproj, dz_r))


ext = extent(5513, 470513, 530301.5, 1220302)
proj = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +towgs84=446.448,-125.157,542.06,0.15,0.247,0.842,-20.489 +units=m +no_defs'
res = c(1000, 1000)
r = raster(ext = ext, resolution = res, crs = proj)
reproj = projectRaster(from = tmpin, crs = proj)
reproj = resample(reproj, r)
vals = reproj@data@vals

