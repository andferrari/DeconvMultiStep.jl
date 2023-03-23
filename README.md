# DeconvMultiStep.jl

Split the deconvolution in 2 step:
1. Reconstruct a low frequency image using short bases visibilities.
2. Add high frenquencies to the low frequency image using large bases visibilities.

An example is available in `notebooks/main_meerkat.ipynb`

`DeconvMultiStep.jl/rascil` folder contains the files for Rascil minor loop:
- `make_fullres.jl`: minor loop when using all the baselines
- `make_lowres.jl`: minor loop when using only *short* baselines
- `make_multistep.jl`: minor loop when using only *long* baselines and output of `make_lowres.jl`
