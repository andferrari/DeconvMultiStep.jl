# reconstruction using only short bases

using DeconvMultiStep
using FITSIO

root = "/root"

psf = read(FITS(joinpath(root, "tmp_psf_low.fits"))[1])
dirty = read(FITS(joinpath(root, "tmp_dirty_low.fits"))[1])

# sky is needed to compute the regularisation parameter
sky = read(FITS(joinpath(root, "sky.fits"))[1])
λ = comp_λ(psf, dirty, sky)

i_lowres = fista(psf, dirty, λ, 200)

f = FITS(joinpath(root, "ferrari_output_low.fits"), "w")
write(f, i_lowres)
close(f)