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
using JLD
using Plots
gr()

bsbi1 = loadtable("data/plantdata/PlantData_87to99.txt")
bsbi1 = @transform bsbi1 {Recording_period = "87to99"}
bsbi2 = loadtable("data/plantdata/PlantData_70to87.txt")
bsbi2 = @transform bsbi2 {Recording_period = "70to87"}
bsbi3 = loadtable("data/plantdata/PlantData_30to70.txt")
bsbi3 = @transform bsbi3 {Recording_period = "30to70"}
bsbi = merge(bsbi1, bsbi2)
bsbi = merge(bsbi, bsbi3)
species = loadtable("data/plantdata/PlantData_Species.txt")
squares = loadtable("data/plantdata/PlantData_Squares.txt")
bsbi = join(bsbi, squares, lkey = :OS_SQUARE, rkey = :OS_SQUARE)
bsbi = join(bsbi, species, lkey = :TAXONNO, rkey = :TAXONNO)
bsbi = @transform bsbi {SppID = :TAXONNO}

# Create reference for UK grid
ref = createRef(1000.0m, 500.0m, 7e5m, 500.0m, 1.25e6m)
bsbi = @transform bsbi {refval = UKplantSim.extractvalues(:EAST * m, :NORTH * m, ref), refid = 1}

neilspecies = loadtable("data/KF_Plant traits_AC.csv")
joinspecies = join(neilspecies, bsbi, lkey = :Species, rkey = :NAME, how = :inner)
joinspecies = @select joinspecies (:Species, :Family, :TAXONNO, :OS_SQUARE, :BRC_SQUARE, :EAST, :NORTH, :TAXON, :NATIVE, :HYBRID, :Recording_period)
joinspecies = @transform joinspecies {lon = BritishNationalGrid.bng2lonlat(:EAST, :NORTH)[1], lat = BritishNationalGrid.bng2lonlat(:EAST, :NORTH)[2]}
FileIO.save("Neil_species.csv", joinspecies)
