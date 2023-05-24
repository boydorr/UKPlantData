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
using Plots
gr()

import UKplantSim.LandCover

crop = readCrop("data/Crop2017.tif")

times = collect(2008year:1month:2017year+11month)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)
active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))

cc = Array(crop.array[:, 1000m .. 1.25e6m])
cc[isnan.(cc)] .= 0
cc[.!active] .= NaN
crop = CropCover(AxisArray(cc, Axis{:northing}(sun.array.axes[1].val), Axis{:easting}(sun.array.axes[2].val)))
plot(crop)
Plots.pdf("CropCover2017.pdf")


lc = readLC("CEH_landcover_2015.tif")
plot(lc)
Plots.pdf("LandCover2017.pdf")

function combineLC(lc::LandCover, cc::CropCover)
    lc = lc.array[:, 0m .. 1.25e6m]
    cropland = findall(cc.array .> 0)
    agland = findall(lc .== 3)
    inter = agland ∩ cropland
    lc[inter] += (cc.array[inter] .+ 18)
    return LandCover(lc)
end

newlc = combineLC(lc, crop)
plot(newlc)
Plots.pdf("LandCropCombo.pdf")

start = JLD.load("StartArray.jld", "start")
start = reshape(start, size(start, 1), 700, 1250)
ag = findall(lc.array .== 3)
start[:, ag] .= 0
crops = zeros(11, 700, 1250)
for i in 1:11
    locs = findall((newlc.array .== 3) .| (newlc.array .== (i+21)))
    crops[i, locs] .= 1e3
end
start = cat(start, crops, dims = 1)
abun = start[end, :, :]
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)
Plots.pdf("StartAbun.pdf")

abun = reshape(start, size(start, 1), 700 * 1250)
abun = norm_sub_alpha(Metacommunity(abun), 1.0)[:diversity]
abun = reshape(abun, 700, 1250)
abun[isnan.(abun)] .= 0
abun[.!active] .= NaN
heatmap(transpose(abun), background_color = :lightblue, background_color_outside = :white, grid = false, color = :algae, aspect_ratio = 1)
Plots.pdf("StartCropDiv.pdf")

using JLD
dir = "HadUK/tas/"
times = collect(2008year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "HadUK/sun/"
sun = readHadUK(dir, "sun", times)

# Take means of 2015 (same as for LC)
meantas2015 = mapslices(mean, tas.array[:, :, 2015year..2015year+11months]./K, dims = 3)[:, :, 1]
meanrainfall2015 = mapslices(mean, rainfall.array[:, :, 2015year..2015year+11months]./mm, dims = 3)[:, :, 1]
meansun2015 = mapslices(mean, (uconvert.(kJ, 1km^2 .* sun.array[:, :, 2015year..2015year+11months] .* 1000*(W/m^2)))./kJ, dims = 3)[:, :, 1]
namean(x) = mean(x[.!isnan.(x)])
nastd(x) = std(x[.!isnan.(x)])

trt_means = map(1:11) do i
    locs = findall(newlc.array .== (i+21))
    return [namean(meantas2015[locs]), namean(meanrainfall2015[locs]), namean(meansun2015[locs])]
end
trt_means = hcat(trt_means...)
JLD.save("Crop_trait_means.jld", "trt_means", trt_means)

trt_stds = map(1:11) do i
    locs = findall(newlc.array .== (i+21))
    return [nastd(meantas2015[locs]), nastd(meanrainfall2015[locs]), nastd(meansun2015[locs])]
end
trt_stds = hcat(trt_stds...)
JLD.save("Crop_trait_stds.jld", "trt_stds", trt_stds)


locs = findall(newlc.array .== (22))
beet_means = [namean(meantas2015[locs]), namean(meanrainfall2015[locs]), namean(meansun2015[locs])]
beet_stds = [nastd(meantas2015[locs]), nastd(meanrainfall2015[locs]), nastd(meansun2015[locs])]

locs = findall((newlc.array .== 29) .| (newlc.array .== 31))
barley_means = [namean(meantas2015[locs]), namean(meanrainfall2015[locs]), namean(meansun2015[locs])]
barley_stds = [nastd(meantas2015[locs]), nastd(meanrainfall2015[locs]), nastd(meansun2015[locs])]


traits = JuliaDB.load("BSBI_had_prefs_UK")
traits = filter(t-> !isnan(t.sun) & !isnan(t.rainfall) & !isnan(t.tas_st) & !isnan(t.rain_st), traits)
traits = filter(t -> (t.rain_st > 0) & (t.tas_st > 0), traits)

push!(rows(traits), (NAME = "Beta vulgaris", tas = beet_means[1], rainfall = beet_means[2], sun = beet_means[3], tas_st = beet_stds[1], rain_st = beet_stds[2]))
push!(rows(traits), (NAME = "Hordeum vulgare", tas = barley_means[1], rainfall = barley_means[2], sun = barley_means[3], tas_st = barley_stds[1], rain_st = barley_stds[2]))
table(traits, pkey = :NAME, copy = false)
JuliaDB.save(traits, "Crop_had_prefs_UK")

crop_names = ["Beta vulgaris", "Vicia faba", "Grass", "Zea mays", "Brassica napus", "Other", "Solanum tuberosum", "Hordeum vulgare", "Triticum aestivum", "Hordeum vulgare", "Triticum aestivum"]
cropids = [findall(x .== select(traits, :NAME)) for x in crop_names]
