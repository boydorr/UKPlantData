using NetCDF
using UKplantSim
using EcoSISTEM.Units
using Unitful
using Unitful.DefaultSymbols
using Plots
gr()
cd("/home/claireh/Documents/UK")
ncinfo("CHESS/prec/chess_precip_201001.nc")


times = collect(2010year:1month:2011year+9month)
dir = "CHESS/prec1/"
prec = readCHESS(dir, "precip", times)

prec = uconvert.(kg/km^2, prec.array .* 1month)

heatmap(transpose(ustrip.(prec[:,:,1])), aspect_ratio =1)
