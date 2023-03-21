using CairoMakie
using FFTW

# Attributes(
#     font = "Chilanka",
#     backgroundcolor = :gray,
#     color = :blue,
#     linestyle = :dot,
#     linewidth = 3
# )

# theme isa Attributes && set_theme!(theme)


function plot_psf(psf::PSF, uv::UV; zoom = 1.0)

    nx, ny  =size(psf.full)
    lx, ly = LinRange(-nx/2, nx/2, nx), LinRange(-ny/2, ny/2, ny)

    fig = Figure(resolution = (600, 800))
    ax = [Axis(fig[i, j], aspect=DataAspect()) for i in 1:3, j in 1:2]

    heatmap!(ax[1,1], lx, ly, uv.full, colormap=["white", "black"])
    hidedecorations!.(ax[1,1])
    ax[1,1].title="UV full"
    hm1 = heatmap!(ax[1,2], lx, ly, psf.full)
    xlims!(ax[1,2], -zoom, zoom)
    ylims!(ax[1,2], -zoom, zoom)
    ax[1,2].title="PSF full"
    Colorbar(fig[1, 3], hm1)

    heatmap!(ax[2,1], lx, ly, uv.low, colormap=["white", "black"])
    hidedecorations!.(ax[2,1])
    ax[1,1].title="UV low"
    r_max = (uv.ℓ + uv.δ/2)*exp.(im*LinRange(0, 2*π, 256)) 
    r_min = (uv.ℓ - uv.δ/2)*exp.(im*LinRange(0, 2*π, 256)) 
    lines!(ax[2,1], real.(r_max), imag.(r_max), color = :red)
    lines!(ax[2,1], real.(r_min), imag.(r_min), color = :red)
    hm2 = heatmap!(ax[2,2], lx, ly, psf.low)
    xlims!(ax[2,2], -zoom, zoom)
    ylims!(ax[2,2], -zoom, zoom)
    ax[2,2].title="PSF low"
    Colorbar(fig[2, 3], hm2)

    heatmap!(ax[3,1], lx, ly, uv.high, colormap=["white", "black"])
    hidedecorations!.(ax[3,1])
    ax[1,1].title="UV high"
    lines!(ax[3,1], real.(r_max), imag.(r_max), color = :red)
    lines!(ax[3,1], real.(r_min), imag.(r_min), color = :red)
    hm3 = heatmap!(ax[3,2], lx, ly, psf.high)
    xlims!(ax[3,2], -zoom, zoom)
    ylims!(ax[3,2], -zoom, zoom)
    ax[3,2].title="PSF high"
    Colorbar(fig[3, 3], hm3)

    fig
end

function plot_dirty(dirty::Dirty, sky::Matrix{T}) where  {T<:Real}

    fig = Figure(resolution = (800, 700))
    ax = [Axis(fig[i, j], aspect=DataAspect()) for i in 1:2, j in 1:2:3]
    hidedecorations!.(ax)

    cmap = cgrad(:turbo, scale=:exp10 )

    hm1 = heatmap!(ax[1,1], sky, colormap=cmap)
    ax[1,1].title="Sky"
    Colorbar(fig[1, 2], hm1, minorticksvisible=true)

    hm2 = heatmap!(ax[1,2], dirty.full, colormap=cmap)
    ax[1,2].title="Dirty full"
    Colorbar(fig[1, 4], hm2, minorticksvisible=true)

    hm3 = heatmap!(ax[2,1], dirty.low, colormap=cmap)
    ax[2,1].title="Dirty low"
    Colorbar(fig[2, 2], hm3, minorticksvisible=true)

    hm4 = heatmap!(ax[2,2], dirty.high, colormap=:turbo)
    ax[2,2].title="Dirty high"
    Colorbar(fig[2, 4], hm4, minorticksvisible=true)

    rowsize!(fig.layout, 1, Aspect(1, 1))
    rowsize!(fig.layout, 2, Aspect(1, 1))

    fig
end

function plot_dirty_(dirty::Dirty, sky::Matrix{T}; δ::Float64=1e-2) where  {T<:Real}

    fig = Figure(resolution = (900, 600))
    ax = [Axis(fig[i, j], aspect=DataAspect()) for i in 1:2, j in [1, 4]]
    hidedecorations!.(ax)


    hm1l = heatmap!(ax[1,1], map(x -> ( x<δ ? NaN : x), sky); colormap =:heat)
    hm1h = heatmap!(ax[1,1], map(x -> ( x>δ ? NaN : x), sky); colormap =:grays)
    Colorbar(fig[1, 2], hm1l, minorticksvisible=true)
    Colorbar(fig[1, 3], hm1h, minorticksvisible=true)
    ax[1,1].title="Sky"

    hm2l = heatmap!(ax[1,2], map(x -> ( x<δ ? NaN : x), dirty.full); colormap =:heat)
    hm2h = heatmap!(ax[1,2], map(x -> ( x>δ ? NaN : x), dirty.full); colormap =:grays)
    Colorbar(fig[1, 5], hm2l, minorticksvisible=true)
    Colorbar(fig[1, 6], hm2h, minorticksvisible=true)
    ax[1,2].title="Dirty full"

    hm3l = heatmap!(ax[2,1], map(x -> ( x<δ ? NaN : x), dirty.low); colormap =:heat)
    hm3h = heatmap!(ax[2,1], map(x -> ( x>δ ? NaN : x), dirty.low); colormap =:grays)
    Colorbar(fig[2, 2], hm3l, minorticksvisible=true)
    Colorbar(fig[2, 3], hm3h, minorticksvisible=true)
    ax[2,1].title="Dirty low"


    hm4l = heatmap!(ax[2,2], map(x -> ( x<δ ? NaN : x), dirty.high); colormap =:heat)
    hm4h = heatmap!(ax[2,2], map(x -> ( x>δ ? NaN : x), dirty.high); colormap =:grays)
    Colorbar(fig[2, 5], hm4l, minorticksvisible=true)
    Colorbar(fig[2, 6], hm4h, minorticksvisible=true)
    ax[2,2].title="Dirty high"

    rowsize!(fig.layout, 1, Aspect(1, 1))
    rowsize!(fig.layout, 2, Aspect(1, 1))

    fig
end

function plot_deconv(dirty::Matrix{T}, sky::Matrix{T}, i_rec::Matrix{T}, 
    mse::Vector{T}, title::Vector{String}) where  {T<:Real}
    
    fig = Figure()

    ax1 = Axis(fig[1, 1], aspect=DataAspect())
    ax1.title = title[1]
    ax2 = Axis(fig[1, 3], aspect=DataAspect())
    ax2.title = title[2]
    ax3 = Axis(fig[2, 1:4], xlabel = "Iterations")
    hidedecorations!(ax1)
    hidedecorations!(ax2)

    cmap = cgrad(:turbo, scale=:exp10 )
    
    hm1 = heatmap!(ax1, dirty, colormap=cmap)
    cbar = Colorbar(fig[1, 2], hm1,  minorticksvisible=true)
    
    hm2 = heatmap!(ax2, i_rec, colormap=cmap)
    cbar = Colorbar(fig[1, 4], hm2,  minorticksvisible=true)
    
    l1 =  lines!(ax3, mse, label = "MSE")
    axislegend(ax3)
    rowsize!(fig.layout, 1, Aspect(1, 1))
    rowsize!(fig.layout, 2, Relative(1/4))
    
    fig
end

function plot_deconv_(dirty::Matrix{T}, sky::Matrix{T}, i_rec::Matrix{T}, 
    mse::Vector{T}, title::Vector{String}; δ::Float64=1e-2) where  {T<:Real}
    
    fig = Figure(resolution = (900, 500))

    ax1 = Axis(fig[1, 1], aspect=DataAspect())
    ax1.title = title[1]
    ax2 = Axis(fig[1, 4], aspect=DataAspect())
    ax2.title = title[2]
    ax3 = Axis(fig[2, 1:6], xlabel = "Iterations")
    hidedecorations!(ax1)
    hidedecorations!(ax2)

    hm1l = heatmap!(ax1, map(x -> ( x<δ ? NaN : x), dirty); colormap =:heat)
    hm1h = heatmap!(ax1, map(x -> ( x>δ ? NaN : x), dirty); colormap =:grays)
    Colorbar(fig[1, 2], hm1l, minorticksvisible=true)
    Colorbar(fig[1, 3], hm1h, minorticksvisible=true)
    

    hm2l = heatmap!(ax2, map(x -> ( x<δ ? NaN : x), i_rec); colormap =:heat)
    hm2h = heatmap!(ax2, map(x -> ( x>δ ? NaN : x), i_rec); colormap =:grays)
    Colorbar(fig[1, 5], hm2l, minorticksvisible=true)
    Colorbar(fig[1, 6], hm2h, minorticksvisible=true)

    
    l1 =  lines!(ax3, mse, label = "MSE")
    axislegend(ax3)
    rowsize!(fig.layout, 1, Aspect(1, 1))
    rowsize!(fig.layout, 2, Relative(1/4))
    
    fig
end

function plot_filters_rad(ℓ::T, δ::T, n_pix::Int64; σ² = 1.0, η² = 1.0, zoom = nothing) where {T<:Real}

    #radial freq. response

    g = filt_rad.(0:n_pix/2-1, ℓ, δ; σ²=σ², η²=η²)
    gl = [tup[1] for tup in g]
    gh = [tup[2] for tup in g]

    fig_rad = Figure(resolution = (600, 400), font = "CMU Serif")
    ax = Axis(fig_rad[1, 1],  xlabel = (L"r"), xlabelsize = 22)
    lines!(ax, gl, label = L"g_\mathcal{L}(r)")
    lines!(ax, gh, label = L"g_\mathcal{H}(r)")
    lines!(ax, σ²*gh.^2+η²*gl.^2, label = L"\sigma^2 g_\mathcal{L}(r)^2 + \eta^2 g_\mathcal{H}(r)^2")

    vlines!(ax, [ℓ+δ]; color = :black, linewidth = 1, linestyle = :dashdot)
    vlines!(ax, [ℓ-δ]; color = :black, linewidth = 1, linestyle = :dashdot, xticks = ([ℓ+δ], ["p"]))
    text!(ax, L"\ell + \delta", position = (ℓ+δ, -0.1))
    text!(ax, L"\ell - \delta", position = (ℓ-δ , -0.1))

    band!([ℓ-δ, ℓ+δ], [0, 0], [1, 1]; color = (:blue, 0.2))
    axislegend(position = :rc)
    zoom === nothing || xlims!(ax, 0, zoom)

    fig_rad

end

function plot_filters2D(G::Filters; zoom = nothing) 

    lx = LinRange(-G.n_pix/2, G.n_pix/2, G.n_pix)

    fig_psf = Figure(resolution = (800, 400), font = "CMU Serif")
    ax1 = Axis(fig_psf[1, 1], aspect=DataAspect(); title="Low pass PSF") 
    ax2 = Axis(fig_psf[1, 3], aspect=DataAspect(); title="High pass PSF") 

    hm1 = heatmap!(ax1, lx, lx, abs.(fftshift(fft(ifftshift(G.LowPass)))))
    hm2 = heatmap!(ax2, lx, lx, abs.(fftshift(fft(ifftshift(G.HighPass)))))

    if zoom !== nothing 
        xlims!(ax1, -zoom, zoom)
        ylims!(ax1, -zoom, zoom)
        xlims!(ax2, -zoom, zoom)
        ylims!(ax2, -zoom, zoom)
    end

    Colorbar(fig_psf[1, 2], hm1)
    Colorbar(fig_psf[1, 4], hm2)

    rowsize!(fig_psf.layout, 1, Aspect(1, 1))

    fig_psf

end


function plot_recon(i, uv, ℓ, δ)
    nx, ny  =size(i)
    lx, ly = LinRange(-nx/2, nx/2, nx), LinRange(-ny/2, ny/2, ny)

    fig = Figure(resolution = (600, 600), font = "CMU Serif")
    ax = Axis(fig[1, 1], aspect=DataAspect(); title = "Frequency domain")

    hm = heatmap!(ax, lx, ly, abs.(fftshift(fft(ifftshift(i)))))
    r_max = (uv.ℓ + uv.δ/2)*exp.(im*LinRange(0, 2*π, 256)) 
    r_min = (uv.ℓ - uv.δ/2)*exp.(im*LinRange(0, 2*π, 256)) 
    lines!(ax, real.(r_max), imag.(r_max), color = :red)
    lines!(ax, real.(r_min), imag.(r_min), color = :red)

    Colorbar(fig[1, 2], hm)

    fig 
end