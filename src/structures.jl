# structures

mutable struct PSF{T <: Matrix{<:Real}}
    full::T 
    low::T
    high::T
    n_pix::Int
    
    function PSF(full::T, low::T, high::T, n_pix::Int) where {T <: Matrix{<:Real}}
        
        @assert size(full) == size(low) == size(high) "PSF's must have the same size"
        (m, n) = size(full)
        @assert m == n == n_pix "PSF's must be square with size n_pix"

        new{T}(full, low, high, n)
    end

end

mutable struct Filters{T <: Matrix{<:Real}}
    LowPass::T 
    HighPass::T
    n_pix::Int
    
    function Filters(LowPass::T, HighPass::T, n_pix::Int) where {T <: Matrix{<:Real}}
        
        @assert size(LowPass) == size(HighPass) "Filters PSF's must have the same size"
        (m, n) = size(LowPass)
        @assert m == n == n_pix "Filters PSF's must be square with size n_pix"

        new{T}(LowPass, HighPass, n)
    end

end

mutable struct Dirty{T <: Matrix{<:Real}}
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

mutable struct UV
    full::Matrix{Bool} 
    low::Matrix{Bool} 
    high::Matrix{Bool} 
    n_pix::Int
    ℓ::Float64
    δ::Float64
    
    function UV(full::Matrix{Bool}, low::Matrix{Bool}, high::Matrix{Bool}, n_pix::Int, ℓ::Float64, δ::Float64)
        
        @assert size(full) == size(low) === size(high) "UV coverages must have the same size"
        (m, n) = size(full)
        @assert m == n == n_pix "PSF's must be square with size n_pix"
        
        new(full, low, high, n, ℓ, δ)
    end
end