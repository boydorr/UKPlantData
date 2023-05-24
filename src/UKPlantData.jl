module UKPlantData

using Unitful
using Unitful.DefaultSymbols

include("units.jl")

include("ClimateTypes.jl")

include("Read.jl")
export readLC, readCrop, readHadUK, readCHESS, readUKCP, readSoils

include("Process.jl")
export OSGR_eastnorth, get_neighbours, extractvalues, createRef, startingArray, upres

end