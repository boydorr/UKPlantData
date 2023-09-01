using UKPlantData
using DataFrames
using DataFrameMacros
using Unitful
using Unitful.DefaultSymbols
using UKPlantData.Units
using AxisArrays
using StatsBase
using Statistics
using CSV
using JLD2
using Distributions
using BritishNationalGrid
using Plots
using DataPipeline

handle = DataPipeline.initialise()

gbif = CSV.read(link_read!(handle, "UKModel/GBIF-records"), DataFrame)

pa = CSV.read(link_read!(handle, "UKModel/PLANTATT-data"), DataFrame, normalizenames=true)
spp_names = pa.Taxon_name
gbif = filter(g -> !ismissing(g.species), gbif)
ukspecies = filter(g-> g.species ∈ spp_names, gbif)
uk = filter(g -> (g.decimalLatitude > 50.0) & (g.decimalLatitude < 58.7) & (g.decimalLongitude > -7.6) & (g.decimalLongitude < 1.7), ukspecies)


# Load UK gbif data and transform to national grid
uk = @transform(uk, :east = BNGPoint(lon = :decimalLongitude, lat = :decimalLatitude).e, :north = BNGPoint(lon = :decimalLongitude, lat = :decimalLatitude).n)

# Create reference for UK grid
ref = createRef(1000.0m, 1000.0m, 7e5m, 1000.0m, 1.3e6m)
uk = @transform(uk, :refval = extractvalues(:east * m, :north * m, ref))

# Read in landcover 2015 data
file = link_read!(handle, "UKModel/LCM")
zippath = joinpath(splitpath(file)[1:end-1])
run(`tar xvf $path -C $zippath`)
newfile = zippath
lc = readLC(newfile)

uk = @transform(uk, :lc = lc.array[:refval])
LC_counts = @groupby(uk, :species)
LC_counts = @combine(LC_counts, :lc = countmap(:lc))
file = link_write!(handle, "LCM-prefs")
JLD2.@save file LC_counts

namean(x) = mean(x[.!isnan.(x)])
nastd(x) = std(x[.!isnan.(x)])

dir = link_read!(handle, "UKModel/HadUK-tas")
times = collect(2021year:1month:2021year+11month)
tas = readHadUK(dir, "tas", times)
dir = link_read!(handle, "UKModel/HadUK-rain")
rainfall = readHadUK(dir, "rainfall", times)
dir = link_read!(handle, "UKModel/HadUK-sun")
sun = readHadUK(dir, "sun", times)

# Take means of 2015 (same as for LC)
meantas2015 = mapslices(mean, tas.array./K, dims = 3)[:, :, 1]
meanrainfall2015 = mapslices(mean, rainfall.array./mm, dims = 3)[:, :, 1]
meansun2015 = mapslices(mean, (uconvert.(kJ, 1km^2 .* sun.array .* 1000*(W/m^2)))./kJ, dims = 3)[:, :, 1]

# Extract reference values
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
uk = @transform(uk, :refval = extractvalues(:east * m, :north * m, ref))
uk = @transform(uk, :tas = meantas2015[:refval], :rainfall = meanrainfall2015[:refval], :sun =  meansun2015[:refval])

# Calculate averages per species and plot as histogram
had_counts = @groupby(uk, :species)
had_counts = @combine(had_counts, :tas = namean(:tas), :rainfall = namean(:rainfall), :sun = namean(:sun),
:tas_st = nastd(:tas), :rain_st = nastd(:rainfall))
file = link_write!(handle, "Had-prefs")
JLD2.@save file had_counts

DataPipeline.finalise(handle)