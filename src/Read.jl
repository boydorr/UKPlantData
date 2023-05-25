using Unitful
using Unitful.DefaultSymbols
using UKPlantData.Units
using AxisArrays
using NetCDF
using JLD2

import Unitful.°, Unitful.°C, Unitful.mm
import ArchGDAL
const AG = ArchGDAL

unitdict = Dict("K" => K, "m" => m, "J m**-2" => J/m^2, "m**3 m**-3" => m^3, "degC" => °C, "mm" => mm, "hour" => u"hr", "kg m-2 s-1" => kg/(m^2*s), "mm/day" => mm/day, "hours since 1970-01-01T00:00:00Z" => Unitful.hr)
"""
    read(f, filename)

Function to read raster file into julia.
"""
function read(f, filename)
    return AG.environment() do
        AG.read(filename) do dataset
            f(dataset)
        end
    end
end

"""
    readfile(file::String)

Function to import a selected file from a path string.
"""
function readfile(file::String, xmin::Unitful.Quantity{Float64} = -180.0°, 
                  xmax::Unitful.Quantity{Float64} = 180.0°,
                  ymin::Unitful.Quantity{Float64} = -90.0°, 
                  ymax::Unitful.Quantity{Float64} = 90.0°)
    txy = [Float64, Int64(1), Int64(1), Float64(1)]
    #
    read(file) do dataset
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    a = Array{txy[1], 2}(undef, txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    lat, long = size(a, 1), size(a, 2);
    step_lat = (xmax - xmin) / lat;
    step_long = (ymax - ymin) / long;

    world = AxisArray(a[:, long:-1:1],
                           Axis{:latitude}(xmin:step_lat:(xmax-step_lat/2.0)),
                           Axis{:longitude}(ymin:step_long:(ymax-step_long/2.0)));

    if txy[1] <: AbstractFloat
        world[isapprox.(world, txy[4])] *= NaN;
    end;
    world
end

"""
    searchdir(path,key)

Function to search a directory `path` using a given `key` string.
"""
searchdir(path,key) = filter(x->occursin(key, x), readdir(path))

"""
    readlc(file::String)

Function to import a selected CEH land cover file from a path string. Optional arguments for extent if not Great Britain.
"""
function readLC(file::String, GB::Bool=true)
    if GB
        xmin = 0; xmax = 7e5; ymin = 0; ymax = 1.3e6
    else
        xmin = 1.8e5; xmax = 3.7e5; ymin = 3e5; ymax = 4.6e5
    end
    txy = [Float64, Int64(1), Int64(1), Float64(1)]
    #
    read(file) do dataset
        #txy[1] = AG.getdatatype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    a = Array{txy[1], 2}(undef, txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    lat, long = size(a, 1), size(a, 2);
    step = abs(xmin - xmax) / lat;
    latitude = (xmin+ step):step:xmax
    longitude = (ymin+ step):step:ymax
    size(longitude)
    lc = AxisArray(a[:, end:-1:1], Axis{:easting}(latitude * m), Axis{:northing}(longitude * m));

    if txy[1] <: AbstractFloat
        lc[lc .== lc[1]] *= NaN;
    end;
    return LandCover(lc)
end

function readCrop(file::String)
    xmin = 0; xmax = 7e5; ymin = 0; ymax = 1.25e6
    txy = [Float64, Int64(1), Int64(1), Float64(1)]
    #
    read(file) do dataset
        #txy[1] = AG.getdatatype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    a = Array{txy[1], 2}(undef, txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    lat, long = size(a, 1), size(a, 2);
    step = abs(xmin - xmax) / lat;
    latitude = (xmin+ step):step:xmax
    longitude = (ymin+ step):step:ymax
    size(longitude)
    lc = AxisArray(a[:, end:-1:1], Axis{:easting}(latitude * m), Axis{:northing}(longitude * m));

    if txy[1] <: AbstractFloat
        lc[lc .== lc[1]] *= NaN;
    end;
    return CropCover(lc)
end

"""
    readHadUK(file::String)

Function to import HadUK data into Julia from particular parameter.
"""
function readHadUK(dir::String, param::String, times::Vector{T}) where T<: Unitful.Time
    files = searchdir(dir, ".nc")
    lat = ncread(joinpath(dir, files[1]), "projection_y_coordinate")
    lon = ncread(joinpath(dir, files[1]), "projection_x_coordinate")
    units = ncgetatt(joinpath(dir, files[1]), param, "units")
    units = unitdict[units]
    array = map(files) do f
        ncread(joinpath(dir, f), param)
    end
    array = cat(dims = 3, array...)
    array[array .≈ ncgetatt(joinpath(dir, files[1]), param, "_FillValue")] .= NaN
    array = array * 1.0 * units

    # If temperature param, need to convert from Kelvin
    if typeof(units) <: Unitful.TemperatureUnits
        array = uconvert.(K, array)
    end
    uk = AxisArray(array, Axis{:easting}(lon * m), Axis{:northing}(lat * m), Axis{:month}(times))
    return HadUK(uk[0.0m..1e6m, 0.0m..1.25e6m, :])
end

"""
    readCHESS(file::String)

Function to import CHESS data into Julia from particular parameter.
"""
function readCHESS(dir::String, param::String, times::Vector{T}) where T<: Unitful.Time
    files = searchdir(dir, ".nc")
    lat = ncread(joinpath(dir, files[1]), "y")
    lon = ncread(joinpath(dir, files[1]), "x")
    units = ncgetatt(joinpath(dir, files[1]), param, "units")
    units = unitdict[units]
    array = map(files) do f
        mapslices(mean, ncread(joinpath(dir, f), param), dims = 3)[:, :, 1]
    end
    array = cat(dims = 3, array...)
    array[array .≈ ncgetatt(joinpath(dir, files[1]), param, "_FillValue")] .= NaN
    array = array * 1.0 * units

    # If temperature param, need to convert from Kelvin
    if typeof(units) <: Unitful.TemperatureUnits
        array = uconvert.(K, array)
    end
    uk = AxisArray(array, Axis{:easting}(lon * m), Axis{:northing}(lat * m), Axis{:month}(times))
    return CHESS(uk[0.0m..1e6m, 0.0m..1.25e6m, :])
end
function readCHESS(file::String, param::String)
    x = reverse(ncread(file, "x"))
    y = ncread(file, "y")
    time = ncread(file, "time") * 1.0
    timeunits = ncgetatt(file, "time", "units")
    time = uconvert.(years, time .* NEWUNITDICT[timeunits] .+ 1970years)
    units = ncgetatt(file, param, "units")
    units = NEWUNITDICT[units]
    array = ncread(file, param) * 1.0
    array[array .≈ ncgetatt(file, param, "_FillValue")] .= NaN
    
    return AxisArray(array * units, Axis{:x}(x * m), Axis{:y}(y * m),
                        Axis{:time}(time))
end
function readCHESS(file::String, param::String, times::Vector{T}, xs::Vector{L}, ys::Vector{L}) where {T <: Unitful.Time, L <: Unitful.Length}
    x = reverse(ncread(file, "x")) * m
    y = ncread(file, "y") * m
    select_xs = findall((x .>= xs[1]) .& (x .<= xs[end]))
    select_ys = findall((y .>= ys[1]) .& (y .<= ys[end]))
    time = ncread(file, "time") * 1.0
    timeunits = ncgetatt(file, "time", "units")
    time = uconvert.(years, time .* NEWUNITDICT[timeunits] .+ 1970years)
    select_times = findall((time .>= times[1]) .& (time .<= times[end]))
    units = ncgetatt(file, param, "units")
    units = NEWUNITDICT[units]
   
    array = NetCDF.readvar(NetCDF.open(file), param, start=[minimum(select_xs), minimum(select_ys), minimum(select_times)],count = [length(select_xs),length(select_ys), length(select_times)])
    array[array .≈ ncgetatt(file, param, "_FillValue")] .= NaN
    array = array * 1.0
    
    return AxisArray(array * units, Axis{:x}(x[select_xs]), Axis{:y}(y[select_ys]),
                        Axis{:time}(time[select_times]))
end

"""
    readUKCP(file::String)

Function to import UKCP18 data into Julia from particular parameter.
"""
function readUKCP(file::String, param::String, times::Vector{T}, active::String) where T <: Unitful.Time
    lat = ncread(file, "projection_y_coordinate")
    lon = ncread(file, "projection_x_coordinate")
    units = ncgetatt(file, param, "units")
    units = unitdict[units]
    array = ncread(file, param)
    array[array .≈ ncgetatt(file, param, "_FillValue")] .= NaN
    addarray = zeros(Float64, 82, 2, 1200, 1)
    array = cat(array, addarray, dims = 2)
    array = array * 1.0 * units
    append!(lat, [lat[end] + 12000.0, lat[end] + 24000.0])

    # If temperature param, need to convert from Kelvin
    if typeof(units) <: Unitful.TemperatureUnits
        array = uconvert.(K, array)
    end
    uk = AxisArray(array[:, :, :, 1], Axis{:easting}(lon * m), Axis{:northing}(lat * m), Axis{:month}(times))
    uk = upres(uk, 12)
    uk = uk[1000.0m..7e5m, 1000.0m..1.25e6m, :]
    active_grid = JLD.load(active, "active") .== 0
    uk[active_grid, :] *= NaN
    return UKCP(uk)
end

"""
    readSoils(file::String)

Function to read in Hutton soil data from file.
"""
function readSoils(file::String)
    xmin = 0; xmax = 7e5; ymin = 0; ymax = 1.25e6
    txy = [Float64, Int64(1), Int64(1), Float64(1)]
    #
    read(file) do dataset
        #txy[1] = AG.getdatatype(AG.getband(dataset, 1))
        txy[2] = AG.width(AG.getband(dataset, 1))
        txy[3] = AG.height(AG.getband(dataset, 1))
        txy[4] = AG.getnodatavalue(AG.getband(dataset, 1))
        print(dataset)
    end

    a = Array{txy[1], 2}(undef, txy[2], txy[3])
    read(file) do dataset
        bd = AG.getband(dataset, 1);
        AG.read!(bd, a);
    end;
    lat, long = size(a, 1), size(a, 2);
    step = abs(xmin - xmax) / lat;
    latitude = (xmin+ step):step:xmax
    longitude = (ymin+ step):step:ymax
    size(longitude)
    soils = AxisArray(a[:, end:-1:1], Axis{:easting}(latitude * m), Axis{:northing}(longitude * m));

    if txy[1] <: AbstractFloat
        soils[soils .== soils[1]] *= NaN;
    end;
    return Soils(soils)
end
