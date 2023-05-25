using UKPlantData
using DataFrames
using DataFrameMacros
using Unitful
using Unitful.DefaultSymbols
using UKPlantData.Units
using AxisArrays
using StatsBase
using CSV
using JLD2
using Distributions
using Plots

# Read in landcover 2015 data
file = "data/raw/CEH_landcover_2015.tif"
lc = readLC(file)

# Read in national plant monitoring scheme and attributes
pr = CSV.read("data/raw/records-2019-10-21/records-2019-10-21.csv", DataFrame, normalizenames=true)
pa = CSV.read("data/raw/PLANTATT_19_Nov_08.csv", DataFrame, normalizenames=true)

# Check for overlap between species
spp_pr = unique(pr.Scientific_name)
spp_pa = pa.Taxon_name
cross_species = spp_pr ∩ spp_pa

# Filter for overlap and convert OSGR into northing/easting
pr = filter(p -> p.Scientific_name ∈ cross_species, pr)
pr = filter(p -> p.OSGR_1km != "", pr)
@transform!(pr, :OSGR_1km = String.(:OSGR_1km))
pa = filter(p -> p.Taxon_name ∈ cross_species, pa)
pr = @transform(pr, :east = OSGR_eastnorth(:OSGR_1km)[1], :north = OSGR_eastnorth(:OSGR_1km)[2])

# Extract grid values as reference
ref = createRef(1000.0m, 1000.0m, 7e5m, 1000.0m, 1.3e6m)
pr = @transform(pr, :refval = extractvalues(:east * m, :north * m, ref))

# Extract lc values
pr = @transform(pr, :lc = lc.array[:refval])
LC_counts = @groupby(pr, :Scientific_name)
LC_counts = @combine(LC_counts, :lc = countmap(:lc))

namean(x) = mean(x[.!isnan.(x)])
nastd(x) = std(x[.!isnan.(x)])
# Load HadUK temperature and rainfall
dir = "data/raw/HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "data/raw/HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "data/raw/HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Take means of 2015 (same as for LC)
meantas2015 = mapslices(namean, tas.array[:, :, 2015year..2015year+11months]./K, dims = 3)[:, :, 1]
meanrainfall2015 = mapslices(namean, rainfall.array[:, :, 2015year..2015year+11months]./mm, dims = 3)[:, :, 1]
meansun2015 = mapslices(mean, (uconvert.(kJ, 1km^2 .* sun.array[:, :, 2015year..2015year+11months] .* 1000*(W/m^2)))./kJ, dims = 3)[:, :, 1]

# Extract reference values
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
pr = @transform(pr, :refval = extractvalues(:east * m, :north * m, ref))
pr = @transform(pr, :tas = meantas2015[:refval], :rainfall = meanrainfall2015[:refval], :sun =  meansun2015[:refval])

# Calculate averages per species and plot as histogram
had_counts = @groupby(pr, :Scientific_name)
had_counts = @combine(had_counts, :tas = namean(:tas), :rainfall = namean(:rainfall), :sun = namean(:sun),
:tas_st = nastd(:tas), :rain_st = nastd(:rainfall))
JLD2.@save "data/final/NBN_atlas_prefs_UK.jld2" had_counts

sppDict = Dict(zip(cross_species, 1:length(cross_species)))
@transform!(pr, :SppID = sppDict[:Scientific_name])
start = startingArray(pr, length(cross_species), 10)

# Correct for collection effort - population density
start2 = Float64.(start)
file = "data/final/UK_density.tif"
density = Array{Float64, 2}(readfile(file, 0.0m, 7e5m, 5e5m, 1.25e6m))
density[isnan.(density)] .= 0
density[density .== 0] .= 1
file = "data/raw/UK.tif"
outline = Array{Float64, 2}(readfile(file, 0.0m, 7e5m, 5e5m, 1.25e6m))
density[outline .== 0] .= NaN
for i in 1:size(start, 1)
    start2[i, :] ./= density[1:end]
end
start2[isnan.(start2)] .= 0
start2 = rand.(Poisson.(start2))

abun = sum(start2, dims = 1)[1, :]
abun = reshape(abun, 700, 1250)
heatmap(abun', background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)

JLD2.@save "data/final/StartArrayNBN.jld2" start2

## Alternatively extract from GBIF on workstation

# For workstation - filter GBIF records to UK
gbif = CSV.read("data/raw/0025546-190918142434337.csv", DataFrame)
JLD2.@save "data/final/UKgbif.jld2"

pa = CSV.read("data/raw/PLANTATT_19_Nov_08.csv", DataFrame, normalizenames=true)
spp_names = pa.Taxon_name
ukspecies = filter(g-> g.species in spp_names, gbif)
JLD2.@save "UKspecies.jld2" ukspecies

# Load UK gbif data and transform to national grid
uk = JLD2.load("UKspecies.jld2", "ukspecies")
uk = @transform(uk, (:east = BNGPoint(lon = :decimallongitude, lat = :decimallatitude).e, :north = BNGPoint(lon = :decimallongitude, lat = :decimallatitude).n))

# Create reference for UK grid
ref = createRef(1000.0m, 1000.0m, 7e5m, 1000.0m, 1.3e6m)
uk = @transform(uk, :refval = extractvalues(:east * m, :north * m, ref))

# Read in landcover 2015 data
file = "data/raw/CEH_landcover_2015.tif"
lc = readLC(file)

uk = @transform(uk, :lc = lc.array[:refval])
LC_counts = @groupby(pr, :Scientific_name)
LC_counts = @combine(LC_counts, :lc = countmap(:lc))
@save "data/final/GBIF_LC_prefs_UK.jld2" LC_counts

@everywhere namean(x) = mean(x[.!isnan.(x)])
@everywhere nastd(x) = std(x[.!isnan.(x)])
@everywhere using Statistics
dir = "data/raw/HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "data/raw/HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "data/raw/HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Take means of 2015 (same as for LC)
meantas2015 = mapslices(mean, tas.array[:, :, 2015year..2015year+11months]./K, dims = 3)[:, :, 1]
meanrainfall2015 = mapslices(mean, rainfall.array[:, :, 2015year..2015year+11months]./mm, dims = 3)[:, :, 1]
meansun2015 = mapslices(mean, (uconvert.(kJ, 1km^2 .* sun.array[:, :, 2015year..2015year+11months] .* 1000*(W/m^2)))./kJ, dims = 3)[:, :, 1]

# Extract reference values
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
uk = @transform(uk, :refval = extractvalues(:east * m, :north * m, ref))
uk = @transform(uk, :tas = meantas2015[:refval], :rainfall = meanrainfall2015[:refval], :sun =  meansun2015[:refval])

# Calculate averages per species and plot as histogram
had_counts = @groupby(UK, :Scientific_name)
had_counts = @combine(had_counts, :tas = namean(:tas), :rainfall = namean(:rainfall), :sun = namean(:sun),
:tas_st = nastd(:tas), :rain_st = nastd(:rainfall))
JLD2.@save "data/final/GBIF_had_prefs_UK.jld2" had_counts

sppDict = Dict(zip(cross_species, 1:length(cross_species)))
@transform!(uk, :SppID = sppDict[:Scientific_name])
start = startingArray(uk, length(cross_species), 10)
file = "data/final/UK_density.tif"
start2 = Float64.(start)
density = Array{Float64, 2}(readfile(file, 0.0m, 7e5m, 5e5m, 1.25e6m))
density[density .== 0] .= 1
for i in 1:size(start, 1)
    start2[i, :] ./= density[1:end]
end
start2[isnan.(start2)] .= 0
start2 = rand.(Poisson.(start2))

abun = sum(start2, dims = 1)[1, :]
abun = reshape(abun, 700, 1250)
heatmap(abun', background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)

JLD2.@save "data/final/StartArrayNBN.jld2" start