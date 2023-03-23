# reconstruction using all bases, 
# i.e. full telescope resolution

using DeconvMultiStep
using FITSIO

root = "/root"

psf = read(FITS(joinpath(root, "tmp_psf_full.fits"))[1])
dirty = read(FITS(joinpath(root, "tmp_dirty_full.fits"))[1])

# sky is needed to compute the regularisation parameter
sky = read(FITS(joinpath(root, "sky.fits"))[1])
λ = comp_λ(psf, dirty, sky)

i_fullres = fista(psf, dirty, λ, 200)

f = FITS(joinpath(root, "ferrari_output_fullres.fits"), "w")
write(f, i_fullres)
close(f)
