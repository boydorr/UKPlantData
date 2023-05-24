using NetCDF
using EcoSISTEM.ClimatePref
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Peatland

times = uconvert.(months, (unique(ncread("data/CHESSLandMonHydEn2015.nc", "time")) .* s) .- 54years)
eastings = unique(ncread("data/CHESSLandMonHydEn2015.nc", "eastings")) .* m
northings = unique(ncread("data/CHESSLandMonHydEn2015.nc", "northings")) .* m
xs = 261000.0m:1000m:266000.0m; ys = 289000.0m:1000m:293000.0m
readCHESS("data/CHESSLandMonHydEn2015.nc", "smcl", times, eastings, northings)