using UKplantSim
using JuliaDB
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using EcoSISTEM.ClimatePref
using StatsBase
using Plots

# Read in landcover 2015 data
file = "data/CEH_landcover_2015.tif"
lc = readLC(file)

# Read in national plant monitoring scheme and attributes
pr = readNPMS("data/records-2019-10-21/records-2019-10-21.csv")
pa = readPlantATT("data/PLANTATT_19_Nov_08.csv")

# Check for overlap between species
spp_pr = unique(select(pr, :Scientific_name))
spp_pa = select(pa, :Taxon_name)
cross_species = spp_pr ∩ spp_pa

# Filter for overlap and convert OSGR into northing/easting
pr = filter(p -> p.Scientific_name ∈ cross_species, pr)
pr = filter(p -> p.OSGR_1km != "", pr)
pa = filter(p -> p.Taxon_name ∈ cross_species, pa)
pr = transform(pr, (:east => :OSGR_1km => x -> OSGR_eastnorth(x)[1], :north => :OSGR_1km => x -> OSGR_eastnorth(x)[2]))

# Extract grid values as reference
ref = createRef(1000.0m, 1000.0m, 7e5m, 1000.0m, 1.3e6m)
n = select(pr, :north) .* m; e = select(pr, :east).* m
pr = transform(pr, :refval => (:east, :north) => x -> extractvalues(x[1] * m, x[2] * m, ref))

# Extract lc values
pr = transform(pr, :lc => :refval => x -> lc.array[x])
LC_counts = groupby(countmap, pr, :Scientific_name, select = :lc)

# Plot as barplot
counts = select(LC_counts, :countmap)
plot(counts)
Plots.pdf("plots/LC_types_UK.pdf")

namean(x) = mean(x[.!isnan.(x)])
# Load HadUK temperature and rainfall
dir = "data/HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "data/HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)

# Take means of 2015 (same as for LC)
meantas2015 = mapslices(namean, tas.array[:, :, 2015year..2015year+11months]./K, dims = 3)[:, :, 1]
meanrainfall2015 = mapslices(namean, rainfall.array[:, :, 2015year..2015year+11months]./mm, dims = 3)[:, :, 1]

# Extract reference values
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
n = select(pr, :north) .* m; e = select(pr, :east).* m
pr = transform(pr, :refval => (:east, :north) => x -> extractvalues(x[1] * m, x[2] * m, ref))
pr = transform(pr, (:tas => :refval => x -> meantas2015[x], :rainfall => :refval => x -> meanrainfall2015[x]))

# Calculate averages per species and plot as histogram
had_counts = groupby((tas = :tas => namean, rainfall = :rainfall => namean), pr, :Scientific_name)
histogram(ustrip.(uconvert.(°C, select(had_counts, :tas).*K)), grid = false, layout = (@layout [a; b]), xlabel = "Mean temperature (°C)", legend = false, guidefontsize = 12, tickfontsize = 12, size = (1000, 1200))
histogram!(select(had_counts, :rainfall), grid = false, subplot = 2, xlabel = "Rainfall (mm)", legend = false, guidefontsize = 12, tickfontsize = 12, color = 2)
Plots.pdf("plots/Climate_prefs_UK.pdf")


## Alternatively extract from GBIF
using UKplantSim
using JuliaDB
using BritishNationalGrid
using Unitful
using EcoSISTEM.Units
using Unitful.DefaultSymbols
using AxisArrays
using Statistics

@everywhere using BritishNationalGrid
@everywhere using Unitful.DefaultSymbols

retrieve_HadUK("tas", 2010, 2017, "HadUK/tas")
retrieve_HadUK("rainfall", 2010, 2017, "../rainfall")
retrieve_HadUK("sun", 2010, 2017, "../sun")
retrieve_HadUK("groundfrost", 2010, 2017, "../groundfrost")
retrieve_HadUK("snowLying", 2010, 2017, "../snow")

# For workstation - filter GBIF records to UK
gbif = load("GBIF_TPL")
uk = filter(g -> (g.decimallatitude > 50.0) & (g.decimallatitude < 58.7) & (g.decimallongitude > -7.6) & (g.decimallongitude < 1.7), gbif)
save(ukspecies, "UKgbif")
uk = load("UKgbif")
pa = readPlantATT("PLANTATT_19_Nov_08.csv")
spp_names = collect(select(pa, :Taxon_name))
ukspecies = filter(g-> g.species in spp_names, uk)
save(ukspecies, "UKspecies")

# Load UK gbif data and transform to national grid
uk = load("UKspecies")
uk = transform(uk, (:east => (:decimallatitude, :decimallongitude) => x -> BNGPoint(lon = x[2], lat = x[1]).e, :north => (:decimallatitude, :decimallongitude) => x -> BNGPoint(lon = x[2], lat = x[1]).n))

# Create reference for UK grid
ref = createRef(1000.0m, 1000.0m, 7e5m, 1000.0m, 1.3e6m)
uk = transform(uk, :refval => (:east, :north) => x -> extractvalues(x[1] * m, x[2] * m, ref))

# Read in landcover 2015 data
file = "CEH_landcover_2015.tif"
lc = readLC(file)

uk = transform(uk, :lc => :refval => x-> lc.array[x])
LC_counts = groupby(countmap, pr, :Scientific_name, select = :lc)
save(LC_counts, "GBIF_LC_prefs_UK")

@everywhere namean(x) = mean(x[.!isnan.(x)])
@everywhere nastd(x) = std(x[.!isnan.(x)])
@everywhere using Statistics
dir = "HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Take means of 2015 (same as for LC)
meantas2015 = mapslices(mean, tas.array[:, :, 2015year..2015year+11months]./K, dims = 3)[:, :, 1]
meanrainfall2015 = mapslices(mean, rainfall.array[:, :, 2015year..2015year+11months]./mm, dims = 3)[:, :, 1]
meansun2015 = mapslices(mean, (uconvert.(kJ, 1km^2 .* sun.array[:, :, 2015year..2015year+11months] .* 1000*(W/m^2)))./kJ, dims = 3)[:, :, 1]

# Extract reference values
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
uk = transform(uk, :refval => (:east, :north) => x -> extractvalues(x[1] * m, x[2] * m, ref))
uk = transform(uk, (:tas => :refval => x -> meantas2015[x], :rainfall => :refval => x -> meanrainfall2015[x], :sun => :refval => x -> meansun2015[x]))

# Calculate averages per species and plot as histogram
had_counts = collect(groupby((tas = :tas => namean, rainfall = :rainfall => namean, sun = :sun => namean, tas_st = :tas_st => nastd, rain_st = :rain_st => nastd), uk, :species))
save(had_counts, "GBIF_had_prefs_UK")



LC_counts = load("data/GBIF_LC_prefs_UK")
had_counts = load("data/GBIF_had_prefs_UK")

histogram(ustrip.(uconvert.(°C, select(had_counts, :tas).*K)), grid = false, layout = (@layout [a; b]), xlabel = "Mean temperature (°C)", legend = false, guidefontsize = 12, tickfontsize = 12, size = (1000, 1200))
histogram!(select(had_counts, :rainfall), grid = false, subplot = 2, xlabel = "Rainfall (mm)", legend = false, guidefontsize = 12, tickfontsize = 12, color = 2)
Plots.pdf("plots/GBIF_climate_prefs_UK.pdf")

counts = countmap(vcat(select(LC_counts, :lc)...))
plot(counts)
Plots.pdf("plots/GBIF_LC_types_UK.pdf")
