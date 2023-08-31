using CSV
using DataPipeline
using DataFrames
using ExcelFiles

# PLANTATT
handle = DataPipeline.initialise()
path = link_read!(handle, "UKModel/PLANTATT")
xlspath = joinpath(splitpath(path)[1:end-1])
if Sys.iswindows()
    run(`tar xvf $path -C $xlspath`)
else
    run(`unzip $path`)
end
xlspath = joinpath(xlspath, "PLANTATT_19_Nov_08.xls")
data = DataFrame(ExcelFiles.load(xlspath, "Data"))
newpath = link_write!(handle, "PLANTATT-data")
CSV.write(newpath, data)

# BRYOATT
path = link_read!(handle, "UKModel/BRYOATT")
xlspath = joinpath(splitpath(path)[1:end-1])
if Sys.iswindows()
    run(`tar xvf $path -C $xlspath`)
else
    run(`unzip $path`)
end
xlspath = joinpath(xlspath, "Bryoatt_updated_2017/Bryoatt_updated_2017.xls")
data = DataFrame(ExcelFiles.load(xlspath, "BRYOATT"))
newpath = link_write!(handle, "BRYOATT-data")
CSV.write(newpath, data)

# GBIF
path = link_read!(handle, "UKModel/GBIF-raw")
zippath = joinpath(splitpath(path)[1:end-1])
if Sys.iswindows()
    run(`tar xvf $path -C $zippath`)
else
    run(`unzip $path`)
end
data = CSV.read(joinpath(zippath, "0025546-190918142434337.csv"), DataFrame, normalizenames=true)
newpath = link_write!(handle, "GBIF-records")
CSV.write(newpath, data)

DataPipeline.finalise(handle)


