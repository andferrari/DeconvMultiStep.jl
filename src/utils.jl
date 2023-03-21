
"""
    snr(x::Array{T}, xₑ::Array{T}; ndigits::Int=3) where {T<:Real}

computes the snr with ndigits
"""
function snr(x::Array{T}, xₑ::Array{T}; ndigits::Int=3) where {T<:Real}
    round(20*log10(norm(x)/norm(x-xₑ)), digits = ndigits)
end

"""
    make_bases(n_ants::Int, n_pix::Int; compress::Float64=0.9)

makes the synthetic bases. The bases are grided inside the square 
    (-n\\_pix/2 ... n\\_pix/2) × (-n\\_pix/2 ... n\\_pix/2) 

- n\\_ants: number of antennas
- n\\_pix: number of pixels 
- compress coefficient shrinks all the bases by a factor compress
"""
function make_bases(n_ants::Int, n_pix::Int; compress::Float64=0.9) 

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

"""
    make_bases(filename::String, n_pix::Int; compress::Float64=0.9)

makes the bases from a file. The bases are grided inside the square 
        (-n\\_pix/2 ... n\\_pix/2) × (-n\\_pix/2 ... n\\_pix/2) 
    
    - filename: contains the bases
    - n\\_pix: number of pixels 
    - compress coefficient shrinks all the bases by a factor compress

    6 lines of header are skipped in filename
"""
function make_bases(filename::String, n_pix::Int; compress::Float64=0.9) 

    @assert compress ≤ 1.0 "Compress coefficient must be ≤ 1"
    # make bases from random gridded antennas

    bases = readdlm(filename, skipstart=6)[:, 1:2]
    round.(Int, compress*n_pix*bases/maximum(abs.(bases))/2)

end

"""
    make_psf(bases::Matrix{Int}, n_pix::Int,  ℓ::Float64, δ::Float64)

    - \\ell: baseline center frequency (in pixels)
    - \\delta: baseline threshold half width (in pixels)
"""
function make_psf(bases::Matrix{Int}, n_pix::Int,  ℓ::Float64, δ::Float64) 

    # split bases

    ind_high = findall(x -> (norm(x) > ℓ - δ/2), eachrow(bases))
    ind_low = findall(x -> (norm(x) < ℓ + δ/2), eachrow(bases))
    bases_high = bases[ind_high, :]
    bases_low = bases[ind_low, :]

    # make uv plans

    uv = UV(make_uv(bases, n_pix), make_uv(bases_low, n_pix), make_uv(bases_high, n_pix), n_pix, ℓ, δ)

    psf_full = uv2psf(uv.full)
    psf_low = uv2psf(uv.low)
    psf_high = uv2psf(uv.high)

    psf = PSF(psf_full, psf_low, psf_high, n_pix)

    return psf, uv
end

"""
    make_uv(bases::Matrix{Int}, n_pix::Int)

make uv coverage from bases. 
uv is a boolean matrix. The center is at (n_pix/2 + 1, n_pix/2 + 1)
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
make_dirty(psf::PSF, sky::Matrix{T}, σ²::Float64) where {T<:Real}

Computes the three dirty images. Add noise on dirty images.
"""
function make_dirty(psf::PSF, sky::Matrix{T}, σ²::Float64) where {T<:Real}
    
    size_ = (psf.n_pix, psf.n_pix)
    (size_ == size(sky)) || error("sky and PSF have different size")

    i = imfilter(sky, psf.full) + sqrt(σ²)*randn(size_)
    i_low = imfilter(sky, psf.low) + sqrt(σ²)*randn(size_)
    i_high = imfilter(sky, psf.high) + sqrt(σ²)*randn(size_)

    Dirty(i, i_low, i_high, psf.n_pix)
end

"""
make_dirty(uv::UV, sky::Matrix{T}, σ²::Float64) where {T<:Real}

Computes the three dirty images. Add noise on visibilities.
"""
function make_dirty(uv::UV, sky::Matrix{T}, σ²::Float64) where {T<:Real}
    
    size_ = (uv.n_pix, uv.n_pix)
    (size_ == size(sky)) || error("sky and uv plane have different size")

    vis_n = fft(sky + sqrt(σ²)*randn(size_))
    i = real.((ifft(fftshift(uv.full).*vis_n)))
    i_low = real.((ifft(fftshift(uv.low).*vis_n)))
    i_high = real.((ifft(fftshift(uv.high).*vis_n)))

    Dirty(i, i_low, i_high, uv.n_pix)
end

# deconv

"""
    imfilter(img::Matrix{Float64}, psf::Matrix{Float64})

convolves img by the psf. It computes a linear convolution using a circular 
convolution after zero padding by half the width

- the center of the psf is at size(psf)/2 + 1
"""
function imfilter_(img::Matrix{T}, psf::Matrix{T}) where {T<:Real}

    @assert size(img) == size(psf) "Image and psf must have the same size"
    n, = size(img)
    np = Int(n/2)

    img_pad = vcat(zeros(np, 2n), hcat(zeros(n, np), img, zeros(n, np)), zeros(np, 2n))
    psf_pad = vcat(zeros(np, 2n), hcat(zeros(n, np), psf, zeros(n, np)), zeros(np, 2n))

    tmp = real(ifft(fft(img_pad).*fft(ifftshift(psf_pad))))
    tmp[np+1:np+n, np+1:np+n]
end


"""
    imfilter(img::Matrix{Float64}, psf::Matrix{Float64})

convolves img by the psf using a circular 

- the center of the psf is at size(psf)/2  + 1
"""
function imfilter(img::Matrix{T}, psf::Matrix{T}) where {T<:Real}

    @assert size(img) == size(psf) "Image and psf must have the same size"
    real(ifft(fft(img).*fft(ifftshift(psf))))
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
    fista(H::Matrix{U} , id::Matrix{U}, λ::Float64, n_iter::Int; wlts::Union{Nothing, Vector{T}}=nothing, η::Union{Nothing, Float64}=nothing, 
    G::Union{Nothing, Filters}=nothing, ip::Union{Nothing, Matrix{U}}=nothing, 
    sky::Union{Nothing, Matrix{U}}=nothing, show_progress=false) where {T<:WT.OrthoWaveletClass, U<:Real}

solve minₓ ||G_high⋅(id - H⋅wlts⋅x)||² + ||G_low⋅(ip - wlt⋅x)||² + λ||x||₁ 
using FISTA algorithm

- id : high frequency dirty image
- n_iter : number of iterations
- η : gradient step. If not provided η is computed from ∇f Lipschitz constant
- ip : low frequency image
- G_low : low pass filter PSF
- G_high : high pass filter PSF

If sky is provided returns (x, mse)
"""
function fista(H::Matrix{U} , id::Matrix{U}, λ::Float64, n_iter::Int; wlts::Union{Nothing, Vector{T}}=nothing, η::Union{Nothing, Float64}=nothing, 
    G::Union{Nothing, Filters}=nothing, ip::Union{Nothing, Matrix{U}}=nothing, 
    sky::Union{Nothing, Matrix{U}}=nothing, show_progress=false) where {T<:WT.OrthoWaveletClass, U<:Real}
    
    # init

    (wlts == nothing) && (wlts = [WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8])
    βₚ = zeros((size(id)...,(length(wlts))...))
    α = βₚ
    tₚ = 1
 
    # precomputations
    H_adj = adj(H)

    if G === nothing
        H2 = imfilter(H, H_adj)
        Hid = imfilter(id, H_adj)
        
        (η == nothing) && (η =compute_step(H; wlts = wlts))
    else 
        Glow = G.LowPass
        Ghigh = G.HighPass
        Glow_adj = adj(Glow)
        Ghigh_adj = adj(Ghigh)

        HG2 = imfilter(imfilter(imfilter(H, Ghigh), Ghigh_adj), H_adj) + imfilter(Glow, Glow_adj)
        HGidip = imfilter(imfilter(imfilter(id, Ghigh), Ghigh_adj), H_adj) + imfilter(imfilter(ip, Glow), Glow_adj)

        (η == nothing) && (η =compute_step(H; wlts = wlts, G=G))
    end

    (sky === nothing) || (mse = Float64[])

    p_bar = Progress(n_iter; enabled=show_progress)
    for k in 1:n_iter

        # compute gradient

        i_ = dwt_decomp_adj(α, wlts)

        if G === nothing
            u = imfilter(i_, H2) - Hid
        else
            u = imfilter(i_, HG2) - HGidip
        end
        ∇f = 2*dwt_decomp(u, wlts)

        # apply prox
        β = shrink.(α - η*∇f, η*λ)

        # update
        t = (1 + sqrt(1+4*tₚ^2))/2
        α = β + ((tₚ-1)/t)*(β - βₚ)
        βₚ = β
        tₚ = t

        if sky ≠ nothing 
            push!(mse, norm(sky - dwt_decomp_adj(α, wlts))^2)
        end

        next!(p_bar)
    end
    (sky == nothing) ? (return dwt_decomp_adj(α, wlts)) : (return dwt_decomp_adj(α, wlts), mse) 
end

# compute step


"""
    compute_step_(H::Matrix{U};  wlts::Union{Nothing, Vector{T}}=nothing, G::Union{Nothing, Filters}=nothing; n_iter::Int=20) where {T<:WT.OrthoWaveletClass, U<:Real}

Compute the optimal gradient step when applying FISTA to
    minₓ ||G_high⋅(id - H⋅wlts⋅x)||² + ||G_low⋅(ip - wlt⋅x)||² + λ||x||₁ 
using power iterations.
"""
function compute_step_(H::Matrix{U};  wlts::Union{Nothing, Vector{T}}=nothing, G::Union{Nothing, Filters}=nothing, n_iter::Int=20) where {T<:WT.OrthoWaveletClass, U<:Real}

    (wlts === nothing) && (wlts = [WT.db1, WT.db2, WT.db3, WT.db4, WT.db5, WT.db6, WT.db7, WT.db8])
    α = randn(size(H)...,length(wlts)...)
    H_adj = adj(H)
 
    # precomputations

    if G === nothing
        H2 = imfilter(H, H_adj)
    else
        Glow = G.LowPass
        Ghigh = G.HighPass
        Glow_adj = adj(Glow)
        Ghigh_adj = adj(Ghigh)
        H2 = imfilter(imfilter(imfilter(H, Ghigh), Ghigh_adj), H_adj) + imfilter(Glow, Glow_adj)
    end

    for n in 1:n_iter

        u = imfilter(dwt_decomp_adj(α, wlts), H2)      
        α_ = dwt_decomp(u, wlts)
        α = α_/norm(α_)

    end

    u = imfilter(dwt_decomp_adj(α, wlts), H2)      
    α_ = dwt_decomp(u, wlts)

    1/(2*dot(α,  α_))
end



"""
    filt_rad(r::Float64, ℓ::Float64, δ::Float64; σ² = 1.0, η² = 1.0)

    - compute the low pass and high pass radial transfer function
"""
function filt_rad(r::Real, ℓ::Real, δ::Real; σ² = 1.0, η² = 1.0)
    
    if (0 ≤ r < ℓ - δ) 

        gl = 1.0/sqrt(η²)
        gh = 0.0

    elseif (ℓ + δ < r) 

        gl = 0.0
        gh = 1.0/sqrt(σ²)

    else

        gl = 0.5*(1 - sin(2*π*(r - ℓ)/(4δ)))
        gh = 0.5*(1 + sin(2*π*(r - ℓ)/(4δ)))
        tmp = sqrt(σ²*gh^2 + η²*gl^2)
        gl /= tmp
        gh /= tmp

    end

    return gl, gh
end


"""
    make_filters(ℓ::Real, δ::Real, n_pix::Int64; σ² = 1.0, η² = 1.0) 

    - computes the low pass and high pass PSFs
"""
function make_filters(ℓ::Real, δ::Real, n_pix::Int64; σ² = 1.0, η² = 1.0) 
    
    # compute the center in FFT plane
    xc = yc = iseven(n_pix) ? n_pix/2 + 1 : (n_pix+1)/2

    Hl = zeros(n_pix, n_pix)
    Hh = zeros(n_pix, n_pix)
    for k in 1:n_pix, l in 1:n_pix
        rep = filt_rad(sqrt((k-xc)^2 + (l-yc)^2), ℓ, δ; σ², η²)
        Hl[k,l] = rep[1]
        Hh[k,l] = rep[2]
    end
    LowPass = real.(ifftshift(ifft(fftshift(Hl))))
    HighPass = real.(ifftshift(ifft(fftshift(Hh))))
    Filters(LowPass, HighPass, n_pix)
end

"""
    compute_step(H::Matrix{U};  wlts::Union{Nothing, Vector{T}}=nothing, G::Union{Nothing, Filters}=nothing) where {T<:WT.OrthoWaveletClass, U<:Real}

Compute the optimal gradient step when applying FISTA to
    minₓ ||G_high⋅(id - H⋅wlts⋅x)||² + ||G_low⋅(ip - wlt⋅x)||² + λ||x||₁ 
for convolution operators
"""
function compute_step(H::Matrix{U};  wlts::Union{Nothing, Vector{T}}=nothing, G::Union{Nothing, Filters}=nothing) where {T<:WT.OrthoWaveletClass, U<:Real}

    (wlts === nothing) ? M = 8 : M = length(wlts)

    if G === nothing
        step = 1/(2*M*maximum(abs2.(fft(H))))
    else
        Glow = G.LowPass
        Ghigh = G.HighPass
        step = 1/(2*M*maximum(abs2.(fft(H).*fft(Ghigh)) + abs2.(fft(Glow))))
    end
    return step
end

