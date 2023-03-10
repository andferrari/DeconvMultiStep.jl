using DeconvMultiStep
using FITSIO

psf = read(FITS("/root/tmp_psf.fits")[1])
dirty = read(FITS("/root/tmp_dirty.fits")[1])

i_e = fista(psf, dirty, 0.5, 50)

f = FITS("/root/ferrari_output.fits", "w")
write(f, i_e)
close(f)