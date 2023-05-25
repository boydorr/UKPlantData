using Revise
using UKPlantData
using UKPlantData.Units
using CSV
using DataFrames
using DataFrameMacros
using BritishNationalGrid
using Unitful
using Unitful.DefaultSymbols
using AxisArrays
using Statistics
using Distributions
using Plots
using JLD2

## CALCULATE SIMULATION START POPULATIONS ##
bsbi = CSV.read("data/raw/plantdata/PlantData_87to99.txt", DataFrame)
species = CSV.read("data/raw/plantdata/PlantData_Species.txt", DataFrame)
squares = CSV.read("data/raw/plantdata/PlantData_Squares.txt", DataFrame)
bsbi = innerjoin(bsbi, squares, on = :OS_SQUARE)
bsbi = innerjoin(bsbi, species, on = :TAXONNO)
bsbi = rename(bsbi, :TAXONNO => :SppID)
pa = CSV.read("data/raw/PLANTATT_19_Nov_08.csv", DataFrame)

cross_species = bsbi.NAME ∩ pa."Taxon name"
bsbi = filter(b -> b.NAME ∈ cross_species, bsbi)

# Create reference for UK grid
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.3e6m)
bsbi = @transform(bsbi, :refval = extractvalues(:EAST* m, :NORTH * m, ref))

# Smooth out plant records across whole UK
start = startingArray(bsbi, nrow(species), 10)
JLD2.@save "data/final/StartArray.jld2" start

# Plot to check
abun = sum(start .> 0, dims = 1)[1, :]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
heatmap(abun', background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)


# Repeat for Scotland
start = startingArray(bsbi, length(species), 10, 1000.0m, 500.0m, 7e5m, 5e5m + 500m, 1.25e6m)
JLD2.@save "Soil_StartArray.jld2" start

## CALCULATE PLANT CLIMATE PREFS ##
namean(x) = mean(x[.!isnan.(x)])
nastd(x) = std(x[.!isnan.(x)])
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
bsbi = @transform(bsbi, :tas = meantas2015[:refval], :rainfall = meanrainfall2015[:refval], :sun =  meansun2015[:refval])

# Calculate averages per species and plot as histogram
bsbi_counts = @groupby(bsbi, :NAME)
bsbi_counts = @combine(bsbi_counts, :tas = namean(:tas), :rainfall = namean(:rainfall), :sun = namean(:sun), 
:tas_st = nastd(:tas), :rain_st = nastd(:rainfall))
JLD2.@save "data/final/BSBI_had_prefs_UK.jld2" bsbi_counts

lc = readLC("data/raw/CEH_landcover_2015.tif")
bsbi = @transform(bsbi, :lc = lc.array[:refval])
lc_counts = @groupby(bsbi, :NAME)
lc_counts = @combine(lc_counts, :lc = countmap(:lc))
JLD2.@save "data/final/BSBI_lc_prefs_uk.jld2" lc_counts

JLD2.@save "data/final/BSBI_prefs_UK.jld2" bsbi

soils = readSoils("data/raw/HuttonSoils.tif")
bsbi = @transform(bsbi, :soil = soils.array[:refval])
soil_counts = @groupby(bsbi, :NAME)
soil_counts = @combine(soil_counts, :soil = [unique(:soil)])
soil_counts = filter(s -> !all(isnan.(s.soil)), soil_counts)
soil_counts = @transform(soil_counts, :soil = Int.(:soil[.!isnan.(:soil)]))
JLD2.@save "data/final/BSBI_soil_prefs_uk.jld2" soil_counts
