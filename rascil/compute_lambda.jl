using DeconvMultiStep
using FITSIO

sky = read(FITS(joinpath(root, ARGS[1]))[1])
λ = comp_λ(psf, dirty, sky)

print(λ)