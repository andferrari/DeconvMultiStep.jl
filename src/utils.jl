# plot functions

function imshow_uv(x; title="") 
    nx, ny  =size(x)
    heatmap(LinRange(-nx/2, nx/2, nx), LinRange(-ny/2, ny/2, ny), 
    x, c=:grays, xlims = (-nx/2, nx/2), ylims=(-ny/2, ny/2), title = title, cbar = false, yflip=true, ratio=1)
end

function imshow_psf(x; zoom=1, title="") 
    nx, ny  =size(x)
    heatmap(LinRange(-nx/2, nx/2, nx), LinRange(-ny/2, ny/2, ny), 
    x, xlims=(-zoom, zoom), ylims=(-zoom, zoom), title = title, cbar = true, yflip=true, ratio=1)
end

imshow(x; title="") = heatmap(x, c=:grays, border = :none, title = title, cbar = false, yflip=true, ratio=1)
snr(x, xₑ) = round(20*log10(norm(x)/norm(x-xₑ)), digits=3)

# structures

struct PSF{T <: Matrix{<:Real}}
    full::T 
    low::T
    high::T
    n_pix::Int
    
    function PSF(full::T, low::T, high::T, n_pix::Int) where {T <: Matrix{<:Real}}
        
        @assert size(full) == size(low) === size(high) "PSF's must have the same size"
        (m, n) = size(full)
        @assert m == n == n_pix "PSF's must be square with size n_pix"

        new{T}(full, low, high, n)
    end

end


struct Dirty{T <: Matrix{<:Real}}
    full::T
    low::T
    high::T 
    n_pix::Int
    
    function Dirty(full::T, low::T, high::T, n_pix::Int) where {T <: Matrix{<:Real}}
        
        @assert size(full) == size(low) === size(high) "Dirty images must have the same size"
        (m, n) = size(full)
        @assert m == n == n_pix "PSF's must be square with size n_pix"
        
        new{T}(full, low, high, n)
    end
end

struct UV
    full::Matrix{Bool} 
    low::Matrix{Bool} 
    high::Matrix{Bool} 
    n_pix::Int
    ℓ::Float64
    
    function UV(full::Matrix{Bool}, low::Matrix{Bool}, high::Matrix{Bool}, n_pix::Int, ℓ::Float64)
        
        @assert size(full) == size(low) === size(high) "UV coverages must have the same size"
        (m, n) = size(full)
        @assert m == n == n_pix "PSF's must be square with size n_pix"
        
        new(full, low, high, n, ℓ)
    end
end

# make data

"""
    make_bases(n_ants::Int, n_pix::Int; compress::Float64=1.0)

makes the synthetic bases. The bases are grided inside the square 
    (-n\\_pix/2 ... n\\_pix/2) × (-n\\_pix/2 ... n\\_pix/2) 

- n\\_ants: number of antennas
- n\\_pix: number of pixels 
- compress coefficient shrinks all the bases by a factor compress
"""
function make_bases(n_ants::Int, n_pix::Int; compress::Float64=1.0) 

    @assert compress ≤ 1.0 "Compress coefficient must be ≤ 1"
    # make bases from random gridded antennas

    pos_ants = randn(n_ants, 2)

    # |antennas position| < compress*n_pix/4

    pos_ants_grid = round.(Int, compress*n_pix*pos_ants/maximum(abs.(pos_ants))/4)

    bases = zeros(Int, Int((n_ants^2 - n_ants)/2), 2)
    q = 0
    for k in 1:n_ants, n in k+1:n_ants
        q += 1
        bases[q, :] = pos_ants_grid[n, :] - pos_ants_grid[k, :]
    end

    return bases
end

function make_bases(filename::String, n_pix::Int; compress::Float64=1.0) 


    @assert compress ≤ 1.0 "Compress coefficient must be ≤ 1"
    # make bases from random gridded antennas

    bases = readdlm(filename, skipstart=6)
    round.(Int, compress*n_pix*bases/maximum(abs.(bases))/2)

end

"""
    make_psf(bases::Matrix{Int}, n_pix::Int,  ℓ::Float64)

makes three PSFs and the related uv coverage:
- full: using all baselines
- low: using baselines shorter that ℓ
- high: using baselines larger that ℓ

- n_ants: number of antennas
- n_pix: number of pixels (no gridding!)
- compress coefficient shrinks all the bases by a factor compress
"""
function make_psf(bases::Matrix{Int}, n_pix::Int,  ℓ::Float64) 

    # split bases

    ind_high = findall(x -> (norm(x) > ℓ), eachrow(bases))
    ind_low = findall(x -> (norm(x) ≤ ℓ), eachrow(bases))
    bases_high = bases[ind_high, :]
    bases_low = bases[ind_low, :]

    # make uv plans

    uv = UV(make_uv(bases, n_pix), make_uv(bases_low, n_pix), make_uv(bases_high, n_pix), n_pix, ℓ)

    psf_full = uv2psf(uv.full)
    psf_low = uv2psf(uv.low)
    psf_high = uv2psf(uv.high)

    psf = PSF(psf_full, psf_low, psf_high, n_pix)

    return psf, uv
end

"""
    make_uv(bases::Matrix{Int}, n_pix::Int)

make uv coverage from bases. 
uv is a boolean matrix. The center is at (n_pix/2, n_pix/2)
"""
function make_uv(bases::Matrix{Int}, n_pix::Int)     
    uv = zeros(Bool, n_pix, n_pix)
    for k in eachrow(bases)
        uv[mod(k[1], n_pix) + 1,  mod(k[2], n_pix) + 1] = true
        uv[mod(n_pix - k[1], n_pix) + 1,  mod(n_pix - k[2], n_pix) + 1] = true
    end
    fftshift(uv)
end

"""
    uv2psf(uv::Matrix{Bool} )

make a psf from uv coverage
"""
uv2psf(uv::Matrix{Bool}) = real.(ifftshift(ifft(fftshift(uv))))

"""
adj(psf::Matrix{T}) where {T <: Real}

Computes the adjoint psf for cirucular convolution
"""
adj(psf::Matrix{T}) where {T <: Real} = circshift(reverse(psf), ntuple(i -> 1, ndims(psf)))

"""
    shrink(v::T, c::T) where {T<:AbstractFloat}

Soft thresholding
"""
shrink(v::T, c::T) where {T<:AbstractFloat} = v > c ? v - c : (v < -c ? v + c : zero(T))

# make dirty images

"""
    make_dirty(psf::PSF, σ²::Float64)

Computes the three dirty images.
"""
function make_dirty(psf::PSF, sky::Matrix{T}, σ²::Float64) where {T<:Real}
    
    size_ = (psf.n_pix, psf.n_pix)
    (size_ == size(sky)) || error("sky and PSF have different size")

    i = imfilter(sky, psf.full) + sqrt(σ²)*randn(size_)
    i_low = imfilter(sky, psf.low) + sqrt(σ²)*randn(size_)
    i_high = imfilter(sky, psf.high) + sqrt(σ²)*randn(size_)

    Dirty(i, i_low, i_high, psf.n_pix)
end

# deconv

"""
    imfilter(img::Matrix{Float64}, psf::Matrix{Float64})

convolves img by the psf. It computes a linear convolution using a circular 
convolution after zero padding by half the width

- the center of the psf is at size(psf)/2 
"""
function imfilter(img::Matrix{T}, psf::Matrix{T}) where {T<:Real}

    @assert size(img) == size(psf) "Image and psf must have the same size"
    n, = size(img)
    np = Int(n/2)

    img_pad = vcat(zeros(np, 2n), hcat(zeros(n, np), img, zeros(n, np)), zeros(np, 2n))
    psf_pad = vcat(zeros(np, 2n), hcat(zeros(n, np), psf, zeros(n, np)), zeros(np, 2n))

    tmp = real(ifft(fft(img_pad).*fft(ifftshift(psf_pad))))
    tmp[np+1:np+n, np+1:np+n]
end

"""
    dwt_decomp(img::Matrix{Float64}, wlts)

computes the DWT of img using all the wavelets in wlts
"""
function dwt_decomp(img::Matrix{U}, wlts::Vector{T}) where {T<:WT.OrthoWaveletClass, U<:Real}
    (nx, ny) = size(img)
    α = zeros(nx, ny, length(wlts))
    α[:,:,1] = dwt(img, wavelet(wlts[1]))
    for b in 2:length(wlts)
        α[:,:,b] = dwt(img, wavelet(wlts[b]))
    end
    return α
end

"""
    dwt_decomp_adj(α::Array{Float64, 3}, wlts)

computes the adjoint operator of dwt_decomp()
"""
function dwt_decomp_adj(α::Array{U, 3}, wlts::Vector{T}) where {T<:WT.OrthoWaveletClass, U<:Real}
    img = idwt(α[:,:,1], wavelet(wlts[1]))
    for b in 2:length(wlts)
        img += idwt(α[:,:,b], wavelet(wlts[b]))
    end
    return img
end

"""
    low_pass(ℓ::Float64, n_pix::Int)

Compute the ideal low pass filter PSF with cutoff frequency at ℓ.
"""
function  low_pass(ℓ::Float64, n_pix::Int)
    H = zeros(n_pix, n_pix)
    for k in -n_pix:n_pix, l in -n_pix:n_pix
        if norm([k, l]) < ℓ
            H[mod(k, n_pix) + 1,  mod(l, n_pix) + 1] = 1
        end
    end

    #hanning(x) = 0.5 - 0.5*cos(2*pi*x/(n_sky-1))
    #wind = ifftshift([hanning(x)*hanning(y) for x in 0:n_sky-1, y in 0:n_sky-1])

    # no need to normalize psf TF(PSF)(0,0) = 1
    real.(ifftshift(ifft(H)))

end

"""
    fista(H, i, wlts, λ, n_iter, η, G = false, i₀ = false, sky = false)

solve minₓ ||i - H⋅wlts⋅x||² + ||G⋅(i₀ - wlt⋅x)||² + λ||x||₁ 
using FISTA algorithm

- n_iter : number of iterations
- η : gradient step
- if G = false the associated term is dropped

If sky is provided returns (x, mse)
"""
function fista(H::Matrix{U} , i::Matrix{U}, λ::Float64, n_iter::Int, η::Float64; wlts::Union{Nothing, Vector{T}}=nothing,
    G::Union{Nothing, Matrix{U}}=nothing, i₀::Union{Nothing, Matrix{U}}=nothing, 
    sky::Union{Nothing, Matrix{U}}=nothing) where {T<:WT.OrthoWaveletClass, U<:Real}
    
    # init

    (wlts == nothing) && (wlts = [WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8])
    βₚ = zeros((size(i)...,(length(wlts))...))
    α = βₚ
    tₚ = 1
    H_adj = adj(H)
    
    (G == nothing) || (G_adj = adj(G))
    (sky == nothing) || (mse = Float64[])

    @showprogress 1 "Computing..." for k in 1:n_iter

        # compute gradient

        i_ = dwt_decomp_adj(α, wlts)
        u = imfilter(imfilter(i_, H) - i, H_adj)
        (G == nothing) || (u += imfilter(imfilter(i_ - i₀, G), G_adj))
        ∇f = 2*dwt_decomp(u, wlts)

        # apply prox
        β = shrink.(α - η*∇f, η*λ)

        # update
        t = (1 + sqrt(1+4*tₚ^2))/2
        α = β + ((tₚ-1)/t)*(β - βₚ)
        βₚ = β
        tₚ = t

        if sky ≠ nothing 
            push!(mse, norm(sky - dwt_decomp_adj(α, wlts)))
        end

    end
    (sky == nothing) ? (return dwt_decomp_adj(α, wlts)) : (return dwt_decomp_adj(α, wlts), mse) 
end

# compute step

"""
    compute_step(H::Matrix{Float64}, wlts; G::Union{Bool, Matrix{Float64}}=false, n_iter=20)

Compute the optimal gradient step when applying FISTA to
    minₓ ||i - H⋅wlts⋅x||² + ||G⋅(i₀ - wlt⋅x)||² + λ||x||₁ 
"""
function compute_step(H::Matrix{U};  wlts::Union{Nothing, Vector{T}}=nothing, G::Union{Nothing, Matrix{U}}=nothing, n_iter::Int=20) where {T<:WT.OrthoWaveletClass, U<:Real}

    (wlts == nothing) && (wlts = [WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8])
    α = randn(size(H)...,length(wlts)...)
    H_adj = adj(H)
    (G == nothing) || (G_adj = adj(G))

    for n in 1:n_iter

        i_ = dwt_decomp_adj(α, wlts)
        u = imfilter(H_adj, imfilter(H, i_))        
        (G == nothing) || (u += imfilter(G_adj, imfilter(G, i_)))
        α_ = dwt_decomp(u, wlts)
        α = α_/norm(α_)

    end

    i_ = dwt_decomp_adj(α, wlts)
    u = imfilter(H_adj, imfilter(H, i_))        
    (G == nothing) || (u += imfilter(G_adj, imfilter(G, i_)))
    α_ = dwt_decomp(u, wlts)

    1/(2*dot(α,  α_))
end