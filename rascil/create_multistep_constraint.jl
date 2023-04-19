using DeconvMultiStep
using FITSIO
using ImageFiltering

root = "/root"

i_prev_lowres = read(FITS(joinpath(root, ARGS[1]))[1])
i_prev_deconv = read(FITS(joinpath(root, ARGS[2]))[1])

# make LP and HP filters
ℓ = 61
δ = 10
n_pix, _ = size(i_prev_lowres)

G = make_filters(ℓ, δ, n_pix)

# make lowres constraint
i_constraint = i_prev_lowres - imfilter(i_prev_deconv, G.LowPass)

f = FITS(joinpath(root, ARGS[3]), "w")
write(f, i_constraint)
close(f)