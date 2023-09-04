using Pkg; Pkg.instantiate()
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

numspecies = length(unique(uk.species))
start = startingArray(uk, nrow(numspecies), 10)
path = link_write!(handle, "StartArray")
JLD2.@save path start

DataPipeline.finalise(handle)