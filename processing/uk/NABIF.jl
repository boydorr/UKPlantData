using UKplantSim
using JuliaDB
using JuliaDBMeta
using BritishNationalGrid
using Unitful
using EcoSISTEM.Units
using Unitful.DefaultSymbols
using AxisArrays
using Statistics
using EcoSISTEM
using Distributions
using Diversity
using Plots

bsbi = loadtable("data/plantdata/PlantData_87to99.txt")
species = loadtable("data/plantdata/PlantData_Species.txt")
squares = loadtable("data/plantdata/PlantData_Squares.txt")
bsbi = join(bsbi, squares, lkey = :OS_SQUARE, rkey = :OS_SQUARE)
bsbi = join(bsbi, species, lkey = :TAXONNO, rkey = :TAXONNO)
bsbi = rename(bsbi, :TAXONNO => :SppID)

bsbi = distribute(bsbi, 12)

pa = readPlantATT("data/PLANTATT_19_Nov_08.csv")

spp_bsbi = unique(collect(select(bsbi, :NAME)))
spp_pa = collect(select(pa, :Taxon_name))
cross_species = spp_bsbi ∩ spp_pa
bsbi = filter(b -> b.NAME ∈ cross_species, bsbi)

# Create reference for UK grid
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.3e6m)
bsbi = transform(bsbi, (:refval => (:EAST, :NORTH) => x -> UKplantSim.extractvalues(x[1] * m, x[2] * m, ref)))
bsbi = insertcols(bsbi, 2, :refid => fill(1, length(bsbi)))

start = startingArray(bsbi, length(species), 10)

abun = norm_sub_alpha(Metacommunity(start), 0.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white,
grid = false, color = :algae, aspect_ratio = 1)


bsbi = transform(bsbi, :cr => :refval => x -> coarsenRef(x, 700, 10))

@everywhere namean(x) = mean(x[.!isnan.(x)])
@everywhere nastd(x) = std(x[.!isnan.(x)])
@everywhere using Statistics
dir = "data/HadUK/tas/"
times = collect(2008year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "data/HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "data/HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Take means of 2015 (same as for LC)
meantas2015 = mapslices(mean, tas.array[:, :, 2015year..2015year+11months]./K, dims = 3)[:, :, 1]
meanrainfall2015 = mapslices(mean, rainfall.array[:, :, 2015year..2015year+11months]./mm, dims = 3)[:, :, 1]
meansun2015 = mapslices(mean, (uconvert.(kJ, 1km^2 .* sun.array[:, :, 2015year..2015year+11months] .* 1000*(W/m^2)))./kJ, dims = 3)[:, :, 1]

# Extract reference values
bsbi = transform(bsbi, (:tas => :refval => x -> meantas2015[x], :rainfall => :refval => x -> meanrainfall2015[x], :sun => :refval => x -> meansun2015[x]))

# Calculate averages per species and plot as histogram
bsbi_counts = collect(groupby((tas = :tas => namean, rainfall = :rainfall => namean, sun = :sun => namean, tas_st = :tas_st => nastd, rain_st = :rain_st => nastd), uk, :species))
save(bsbi_counts, "BSBI_had_prefs_UK")

lc = readLC("data/CEH_landcover_2015.tif")
bsbi = transform(bsbi, :lc => :refval => x -> lc.array[x])
LC_counts = collect(groupby(countmap, pr, :NAME, select = :lc))
save(LC_counts, "BSBI_lc_prefs_uk")

save(bsbi, "BSBI_prefs_UK")

soils = readSoils("data/HuttonSoils.tif")
bsbi = transform(bsbi, :soil => :refval => x -> soils.array[x])
soil_counts = collect(groupby(unique, bsbi, :NAME, select =:soil))
soil_counts = filter(s -> !all(isnan.(s.unique)), soil_counts)
save(soil_counts, "BSBI_soil_prefs_uk")
