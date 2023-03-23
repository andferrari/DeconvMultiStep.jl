using DeconvMultiStep
using FITSIO

root = "/root"

psf = read(FITS(joinpath(root, "tmp_psf_high.fits"))[1])
dirty = read(FITS(joinpath(root, "tmp_dirty_high.fits"))[1])
i_lowres = read(FITS(joinpath(root, "ferrari_output_low.fits"))[1])

# make LP and HP filters
ℓ = 61
δ = 10
n_pix, _ = size(psf)

G = make_filters(ℓ, δ, n_pix)

# sky is needed to compute the regularisation parameter
sky = read(FITS(joinpath(root, "sky.fits"))[1])
λ = comp_λ(psf, dirty, sky, G, i_lowres)

# reconstruction
i_multistep = fista(psf, dirty, λ, 200; G=G, ip = i_lowres)

f = FITS(joinpath(root, "ferrari_output_multistep.fits"), "w")
write(f, i_multistep)
close(f)