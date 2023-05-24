using UKplantSim
using JuliaDB
using BritishNationalGrid
using Unitful
using EcoSISTEM.Units
using Unitful.DefaultSymbols
using AxisArrays
using Statistics
using StatsBase
using EcoSISTEM
using Distributions
using Diversity
using JLD
using Plots

bsbi = loadtable("data/plantdata/PlantData_87to99.txt")
species = loadtable("data/plantdata/PlantData_Species.txt")
squares = loadtable("data/plantdata/PlantData_Squares.txt")
bsbi = join(bsbi, squares, lkey = :OS_SQUARE, rkey = :OS_SQUARE)
bsbi = join(bsbi, species, lkey = :TAXONNO, rkey = :TAXONNO)
bsbi = rename(bsbi, :TAXONNO => :SppID)

# Create reference for UK grid
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.3e6m)
bsbi = transform(bsbi, (:refval => (:EAST, :NORTH) => x -> UKplantSim.extractvalues(x[1] * m, x[1] * m, ref)))
bsbi = insertcols(bsbi, 2, :refid => fill(1, length(bsbi)))

file = "CEH_landcover_2015.tif"
lc = readLC(file)
bsbi = transform(bsbi, :lc => lc.array[:refval])

function top_ten(x)
  count = collect(values(countmap(x)))
  species = collect(keys(countmap(x)))
  thresh = 0.9 * maximum(count)
  top_species = species[sortperm(count, rev=true)]
  return top_species[count[sortperm(count, rev=true)] .> thresh]
end
bsbi_lc = groupby((species = :NAME => top_ten), bsbi, :lc)

broadleavedwoods = filter(b -> b.lc == 1, bsbi)
filter(b -> b.NAME == "Quercus robur", broadleavedwoods)
spp = select(broadleavedwoods, :species)[1]
count = select(broadleavedwoods, :counts)[1]
spp[sortperm(count, rev=true)]
