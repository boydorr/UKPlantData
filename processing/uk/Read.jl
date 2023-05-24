using Unitful
using Unitful.DefaultSymbols
using EcoSISTEM.Units
using AxisArrays
using NetCDF
using JuliaDB
using JLD

import Unitful.°, Unitful.°C, Unitful.mm
import ArchGDAL
import EcoSISTEM.ClimatePref.read
import EcoSISTEM.ClimatePref.searchdir
const AG = ArchGDAL

unitdict = Dict("K" => K, "m" => m, "J m**-2" => J/m^2, "m**3 m**-3" => m^3, "degC" => °C, "mm" => mm, "hour" => u"hr", "kg m-2 s-1" => kg/(m^2*s), "mm/day" => mm/day)

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

function upres(aa::AxisArray{T, 3} where T, rescale::Int64)
    grid = size(aa)
    grid = (grid[1] .* rescale, grid[2] .* rescale, grid[3])
    array = Array{typeof(aa[1]), 3}(undef, grid)
    map(1:grid[3]) do time
        for x in 1:size(aa, 1)
            for y in 1:size(aa, 2)
        array[(rescale*x-(rescale-1)):(rescale*x),
            (rescale*y-(rescale - 1)):(rescale*y), time] .= aa[x, y, time]
            end
        end
    end
    lon = aa.axes[1].val
    step = lon[2] - lon[1]
    smallstep = (lon[2] - lon[1]) / rescale
    newlon = collect(lon[1]+smallstep:smallstep:(lon[end]+step))
    lat = aa.axes[2].val
    smallstep = (lat[2] - lat[1]) / rescale
    newlat = collect(lat[1]+smallstep:smallstep:(lat[end]+step))
    return AxisArray(array,
        Axis{:longitude}(newlon),
        Axis{:latitude}(newlat),
        Axis{:time}(aa.axes[3].val))
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
    readPlantATT(file::String)

Function to import PlantATT data as a JuliaDB table.
"""
function readPlantATT(file::String)
    coldict = Dict(:BRC_code => String, :Taxon_name => String, :Fam => String, :FamA => String, :OrdA => String, :NS => String, :CS => String, :RS => String, :Chg => Float64, :Hght => Float64, :Len => Float64, :P1 => String, :P2 => String, :LF1 => String, :LF2 => String, :W => String, :Clone1 => String, :Clone2 => String, :E1 => Int64, :E2 => Int64, :C => String, :NBI => Int64, :NEur => String, :SBI => Int64, :SEur => String, :Origin => String, :GB => Int64, :IR => Int64, :CI => Int64, :Tjan => Float64, :Tjul => Float64, :Prec => Int64, :Co => String, :Br_Habitats => Vector{Int64}, :L => Int64, :F => Int64, :R => Int64, :N => Int64, :S => Int64)
    return loadtable(file, colparsers = coldict)
end

"""
    readNPMS(file::String)

Function to import National Plant Monitoring Scheme data as a JuliaDB table.
"""
function readNPMS(file::String)
    return loadtable(file)
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
