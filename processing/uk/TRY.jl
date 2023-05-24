using JuliaDB

species = loadtable("plantdata/PlantData_Species.txt")
tryspecies = loadtable("plantdata/TryAccSpecies.csv")

species_joined = join(species, tryspecies, lkey = :NAME, rkey = :AccSpeciesName)
print(select(species_joined, :AccSpeciesID)[1:200])
