# reconstruction using all bases, 
# i.e. full telescope resolution

using DeconvMultiStep
using FITSIO

root = "/root"

psf = read(FITS(joinpath(root, "tmp_psf_full.fits"))[1])
dirty = read(FITS(joinpath(root, "tmp_dirty_full.fits"))[1])

i_fullres = fista(psf, dirty, parse(Float64, ARGS[1]), 100)

f = FITS(joinpath(root, "ferrari_output_fullres.fits"), "w")
write(f, i_fullres)
close(f)
