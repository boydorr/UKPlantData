using UKplantSim
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Plots
gr()
cd("/home/claireh/Documents/UK/")

# dir = "/home/claireh/Documents/UK/HadUK/sun/"
# times = collect(2008year:1month:2017year+11month)
# sun = readHadUK(dir, "sun", times)
# active = Array{Bool, 2}(.!isnan.(sun.array[:, :, 1]))
# JLD.save("/home/claireh/Documents/UK/UK_active.jld", "active", Int.(active))

file = "UKCP18/pr_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc"

times = collect(1980year+1month:1month:2080year)
prec = readUKCP(file, "pr", times, "UK_active.jld")
rainfall = uconvert.(mm/month, prec.array) .* month

anim = @animate for i=12:12:1200
    decrain = ustrip.(rainfall[:, :, (i-11):i])
    decrain = transpose(mapslices(mean, decrain, dims =3)[:, :, 1])
    heatmap(decrain, background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, clim = (0, 500))
end
gif(anim, "rainfall.gif", fps = 5)

file = "UKCP18/tas_rcp85_land-rcm_uk_12km_01_mon_198012-208011.nc"
times = collect(1980year+1month:1month:2080year)
tas = readUKCP(file, "tas", times, "UK_active.jld")

anim = @animate for i=12:12:1200
    years = 1981:2080
    yr = years[Int(i/12)]
    dectas = ustrip.(tas.array[:, :, (i-11):i])
    dectas = transpose(mapslices(mean, dectas, dims =3)[:, :, 1])
    heatmap(dectas, background_color = :lightblue, background_color_outside = :white, grid = false, aspect_ratio = 1, clim = (275, 290), title = "Year $yr")
end
gif(anim, "temp.gif", fps = 5)
