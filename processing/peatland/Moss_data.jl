using CSV
using DataFrames
using JLD2
using DataPipeline

handle = DataPipeline.initialise()

path = link_read!(handle, "UKModel/BRYOATT-data")
moss_spp = CSV.File(path, normalizenames=true)
moss_spp = DataFrame(moss_spp)
moss_spp = filter(x -> !ismissing(x.Taxon_name), moss_spp)

path = link_read!(handle, "PeatModel/PrioritySpecies")
peat_spp = CSV.File(path, normalizenames=true)
peat_spp = DataFrame(peat_spp)

# Filter for species in the compiled peatland list
peat_spp_names = rstrip.(peat_spp[!, :Species_name_or_special_variable])
moss_spp = filter(row -> row.Taxon_name ∈ peat_spp_names, moss_spp)
# Filter for species found in Wales
moss_spp = filter(row -> !ismissing(row.W), moss_spp)
# Filter for moss
moss_spp = filter(row -> row.ML == "M", moss_spp)
# Filter for species that grow on peat
moss_spp = filter(row -> !ismissing(row.PT), moss_spp)

# Check if mosses are found in area according to NBN
path = link_read!(handle, "PeatModel/NBNspecies")
cf_bryo = CSV.File(path, normalizenames=true)
cf_bryo = DataFrame(cf_bryo)
cf_spp = cf_bryo[!, :Species_Name]
moss_spp = filter(row -> row.Taxon_name ∈ cf_spp, moss_spp)
#JLD2.save("data/Peat_30_moss.jld2", "moss_spp", moss_spp)

# Add in some extra species from Cors Fochno reports and monitoring
path = link_read!(handle, "UKModel/BRYOATT-data")
og_spp = CSV.File(path, normalizenames=true)
og_spp = DataFrame(og_spp)
og_spp = filter(x -> !ismissing(x.Taxon_name), og_spp)
path = link_read!(handle, "PeatModel/CFspecies")
cf_spp = DataFrame(CSV.File(path, normalizenames=true))
extra_spp = filter(x -> x.Taxon_name ∈ cf_spp[!, :Name], og_spp)
moss_spp = vcat(moss_spp, extra_spp)
unique!(moss_spp)

path = link_write!(handle, "Moss30spp")
@save path moss_spp

DataPipeline.finalise(handle)

# peat_spp = CSV.File("data/Peatland30spp.csv", normalizenames=true)
# peat_spp = DataFrame(peat_spp)
# plant_spp = CSV.File("data/PLANTATT_19_Nov_08.csv", normalizenames=true)
# plant_spp = DataFrame(plant_spp)
# plant_spp = filter(x -> !ismissing(x.Taxon_name), plant_spp)
# plant_spp = filter(x -> x.Taxon_name ∈ peat_spp.Species, plant_spp)

# plant_spp.Br_Habitats
