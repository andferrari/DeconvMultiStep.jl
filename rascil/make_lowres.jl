using DeconvMultiStep
using FITSIO

root = "/root"

psf = read(FITS(joinpath(root, "tmp_psf_low.fits"))[1])
dirty = read(FITS(joinpath(root, "tmp_dirty_low.fits"))[1])

i_lowres = fista(psf, dirty, 0.5, 50)

f = FITS(joinpath(root, "ferrari_output_low.fits"), "w")
write(f, i_lowres)
close(f)