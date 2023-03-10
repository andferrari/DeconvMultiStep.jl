using DeconvMultiStep
using FITSIO

psf = read(FITS("/root/tmp_psf_high.fits")[1])
dirty = read(FITS("/root/tmp_dirty_high.fits")[1])
i_low = read(FITS("/root/ferrari_output_low.fits")[1])
n_pix, _ = size(psf)

G = make_filters(ℓ, δ, n_pix)
i_multistep, mse = fista(psf, dirty, 1e-4, 200; G=G, ip = i_low)

f = FITS("/root/ferrari_output_full.fits", "w")
write(f, i_e)
close(f)
