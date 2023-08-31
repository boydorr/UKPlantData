# using Revise
using CSV
using DataFrames
# using Dates
using Plots
using JLD2
using Unitful
using Unitful.DefaultSymbols
using UKPlantData
using UKPlantData.Units
using Statistics
using DataPipeline
using NetCDF

handle = DataPipeline.initialise()

# rainfall = CSV.read("../../data/raw/RainfallCF.csv", DataFrame)
# missing_rows = ismissing.(rainfall[!, "AWS CF"])
# rainfall[missing_rows, "AWS CF"] .= rainfall[missing_rows, "CF correction (88.3%)"]

# rainfall[!, :DATE] .= Date.(rainfall[!, :DATE], "dd/mm/yyyy")
# rainfall[!, :month] = Dates.month.(rainfall[!, :DATE])
# rainfall[!, :year] = Dates.year.(rainfall[!, :DATE])
# gdf = groupby(rainfall, [:year, :month])
# monthly_rain = combine(gdf, "AWS CF" => (x -> sum(skipmissing(x))) => :sum)
# plot(monthly_rain[!, :sum])

file = "data/raw/CF_outline.tif"
cf = readfile(file, 261000.0m, 266000.0m, 289000.0m, 293000.0m)
grid = size(cf)
# prec_array = fill(0.0mm, grid[1], grid[2], nrow(monthly_rain))
# for i in 1:nrow(monthly_rain)
#     prec_array[:, :, i] .= fill(monthly_rain[i, :sum], grid) * mm
# end

# bud = WaterTimeBudget(prec_array, 1)
# @save "data/final/RainfallBudget.jld2" bud
path = link_read!(handle, "PeatModel/rainfall")
param = "pr"
times = 2010years+1month:1Units.month:2080year
xs = 261000.0m:1m:266000.0m; ys = 289000.0m:1m:293000.0m
pr = readCHESS(path, param, collect(times), collect(xs), collect(ys))

# 1kg water per 1m^2 is 1mm thick
pr = uconvert.(mm, (ustrip.(pr) .* mm/s) .* 1month)
pr = mean(pr, dims = (1, 2))[1, 1, :]
# plot!(1:843, ustrip.(pr))

prec_array2 = fill(0.0mm, grid[1], grid[2], length(pr))
for i in eachindex(pr)
    prec_array2[:, :, i] .= pr[i]
end

path = link_write!(handle, "RainfallFuture")
bud = prec_array2
@save path bud

path = link_write!(handle, "RainfallBurnin")
bud = prec_array2[:, :, 1:120]
@save path bud

DataPipeline.finalise(handle)