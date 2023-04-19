# reconstruction using only short bases

using DeconvMultiStep
using FITSIO

root = "/root"

psf = read(FITS(joinpath(root, "tmp_psf_low.fits"))[1])
dirty = read(FITS(joinpath(root, "tmp_dirty_low.fits"))[1])

i_lowres = fista(psf, dirty, parse(Float64, ARGS[1]), 100)

f = FITS(joinpath(root, "ferrari_output_low.fits"), "w")
write(f, i_lowres)
close(f)
