using DeconvMultiStep
using FITSIO

root = "/root"

psf = read(FITS(joinpath(root, "tmp_psf_high.fits"))[1])
dirty = read(FITS(joinpath(root, "tmp_dirty_high.fits"))[1])
i_lowres = read(FITS(joinpath(root, "ferrari_output_low.fits"))[1])
n_pix, _ = size(psf)

G = make_filters(5, 50, n_pix)
i_fullres = fista(psf, dirty, 1e-4, 200; G=G, ip = i_lowres)

f = FITS(joinpath(root, "ferrari_output_full.fits"), "w")
write(f, i_fullres)
close(f)
