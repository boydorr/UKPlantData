using Revise
using UKPlantData
using DataFrames
using DataFrameMacros
using BritishNationalGrid
using Unitful
using Unitful.DefaultSymbols
using UKPlantData.Units
using AxisArrays
using Statistics
using Distributions
using Plots
using JLD2

crop = readCrop("data/raw/Crop2017.tif")
lc = readLC("data/raw/CEH_landcover_2015.tif")
newlc = combineLC(lc, crop)
plot(newlc)

JLD2.@load "data/final/StartArray.jld2"
start = reshape(start, size(start, 1), 700, 1250)
ag = findall(lc.array .== 3)
start[:, ag] .= 0
crops = zeros(11, 700, 1250)
for i in 1:11
    locs = findall((newlc.array .== 3) .| (newlc.array .== (i+21)))
    crops[i, locs] .= 1e3
end
start = cat(start, crops, dims = 1)
@JLD2.save "data/final/StartArrayCrops.jld2" start

dir = "data/raw/HadUK/tas/"
times = collect(2010year:1month:2017year+11month)
tas = readHadUK(dir, "tas", times)
dir = "data/raw/HadUK/rainfall/"
rainfall = readHadUK(dir, "rainfall", times)
dir = "data/raw/HadUK/sun/"
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
JLD2.@save "data/final/Crop_trait_means.jld2" trt_means

trt_stds = map(1:11) do i
    locs = findall(newlc.array .== (i+21))
    return [nastd(meantas2015[locs]), nastd(meanrainfall2015[locs]), nastd(meansun2015[locs])]
end
trt_stds = hcat(trt_stds...)
JLD2.@save "data/final/Crop_trait_stds.jld2" trt_stds


locs = findall(newlc.array .== (22))
beet_means = [namean(meantas2015[locs]), namean(meanrainfall2015[locs]), namean(meansun2015[locs])]
beet_stds = [nastd(meantas2015[locs]), nastd(meanrainfall2015[locs]), nastd(meansun2015[locs])]

locs = findall((newlc.array .== 29) .| (newlc.array .== 31))
barley_means = [namean(meantas2015[locs]), namean(meanrainfall2015[locs]), namean(meansun2015[locs])]
barley_stds = [nastd(meantas2015[locs]), nastd(meanrainfall2015[locs]), nastd(meansun2015[locs])]


traits = JLD2.load("data/final/BSBI_had_prefs_UK.jld2", "bsbi_counts")
traits = filter(t-> !isnan(t.sun) & !isnan(t.rainfall) & !isnan(t.tas_st) & !isnan(t.rain_st), traits)
traits = filter(t -> (t.rain_st > 0) & (t.tas_st > 0), traits)

append!(traits, DataFrame(NAME = "Beta vulgaris", tas = beet_means[1], rainfall = beet_means[2], sun = beet_means[3], tas_st = beet_stds[1], rain_st = beet_stds[2]))
append!(traits, DataFrame(NAME = "Hordeum vulgare", tas = barley_means[1], rainfall = barley_means[2], sun = barley_means[3], tas_st = barley_stds[1], rain_st = barley_stds[2]))
JLD2.@save "data/final/Crop_had_prefs_UK.jld2" traits
