using UKplantSim
using JuliaDB
using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using EcoSISTEM.ClimatePref
using Plots
plotlyjs()

retrieve_HadUK("tas", 2008, 2017, "/Users/claireh/Documents/PhD/GIT/UKplantSim/data/HadUK/tas")
retrieve_HadUK("rainfall", 2008, 2017, "/Users/claireh/Documents/PhD/GIT/UKplantSim/data/HadUK/rainfall")
retrieve_HadUK("sun", 2008, 2017, "/Users/claireh/Documents/PhD/GIT/UKplantSim/data/HadUK/sun")
retrieve_HadUK("groundfrost", 2008, 2017, "/Users/claireh/Documents/PhD/GIT/UKplantSim/data/HadUK/groundfrost")

dir = "data/HadUK/tas/"
times = collect(1980year:1month:1990year+11month)
tas1980 = readHadUK(dir, "tas", times)

plot(tas1980, 1980year + March, clim = (270, 290))

file = "data/CEH_landcover_2015.tif"
lc = readLC(file)

file = "data/CEH_landcover_2015_ireland.tif"
lci = readLC(file, false)

plot(lc)

pr = readNPMS("data/records-2019-10-21/records-2019-10-21.csv")
pa = readPlantATT("data/PLANTATT_19_Nov_08.csv")

spp_pr = unique(select(pr, :Scientific_name))
spp_pa = select(pa, :Taxon_name)
cross_species = spp_pr ∩ spp_pa
pr = filter(p -> p.Scientific_name ∈ cross_species, pr)
pr = filter(p -> p.OSGR_1km != "", pr)
pa = filter(p -> p.Taxon_name ∈ cross_species, pa)

using PhyloNetworks
tree = readTopology("data/Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
spp_pa ∩ tip_names
