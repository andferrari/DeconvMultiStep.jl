# reconstruction using only short bases

using DeconvMultiStep
using FITSIO

root = "/root"

psf = read(FITS(joinpath(root, ARGS[2]))[1])
dirty = read(FITS(joinpath(root, ARGS[3]))[1])

i_lowres = fista(psf, dirty, parse(Float64, ARGS[1]), 100)

f = FITS(joinpath(root, ARGS[4]), "w")
write(f, i_lowres)
close(f)
