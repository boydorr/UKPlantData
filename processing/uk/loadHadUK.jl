using Revise
using UKPlantData
using UKPlantData.Units
using DataFrames
using DataFrameMacros
using Unitful
using Unitful.DefaultSymbols
using CSV
using AxisArrays
using Plots

retrieve_HadUK("tas", 2008, 2017, "/data/raw/HadUK/tas")
retrieve_HadUK("rainfall", 2008, 2017, "/data/raw/HadUK/rainfall")
retrieve_HadUK("sun", 2008, 2017, "/data/raw/HadUK/sun")
retrieve_HadUK("groundfrost", 2008, 2017, "/data/raw/HadUK/groundfrost")

dir = "data/raw/HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)

plot(tas, 2010year + March, clim = (270, 290))

file = "data/raw/CEH_landcover_2015.tif"
lc = readLC(file)

file = "data/raw/CEH_landcover_2015_ireland.tif"
lci = readLC(file, false)

plot(lc)

pr = CSV.read("data/raw/records-2019-10-21/records-2019-10-21.csv", DataFrame, normalizenames=true)
pa = CSV.read("data/raw/PLANTATT_19_Nov_08.csv", DataFrame, normalizenames=true)

spp_pr = unique(pr.Scientific_name)
spp_pa = pa.Taxon_name
cross_species = spp_pr ∩ spp_pa
pr = filter(p -> p.Scientific_name ∈ cross_species, pr)
pr = filter(p -> p.OSGR_1km != "", pr)
@transform!(pr, :OSGR_1km = String.(:OSGR_1km))
pa = filter(p -> p.Taxon_name ∈ cross_species, pa)
pr = @transform(pr, :east = OSGR_eastnorth(:OSGR_1km)[1], :north = OSGR_eastnorth(:OSGR_1km)[2])

using PhyloNetworks
tree = readTopology("data/raw/Qian2016.tree")
tipnames = tipLabels(tree)
tip_names = join.(split.(tipnames, "_"), " ")
spp_pa ∩ tip_names
