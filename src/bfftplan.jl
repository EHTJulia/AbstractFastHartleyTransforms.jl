##############################################################################
# This is a generic FHT Plan using BFFT (unnormalized backward FFT)

struct BFFTPlan{T, P} <: DHTPlan{T}
    bfftplan::P
    function BFFTPlan(bfftplan::Plan)
        F = eltype(bfftplan)
        P = typeof(bfftplan)
        T = fieldtypes(F)[1]
        return new{T, P}(bfftplan)
    end
end
size(p::BFFTPlan) = size(p.bfftplan)
fftdims(p::BFFTPlan) = fftdims(p.bfftplan)

function LinearAlgebra.mul!(
        y::AbstractArray{T}, p::BFFTPlan{T, P}, x::AbstractArray{T}
    ) where {T, P}
    fx = p.bfftplan * x
    return y .= real.(fx) .+ imag.(fx)
end

function *(p::BFFTPlan{T, P}, x::AbstractArray{T}) where {T, P}
    z = similar(x)
    return LinearAlgebra.mul!(z, p, x)
end

function plan_inv(p::BFFTPlan)
    return ScaledDHTPlan(p, normalization(eltype(p), size(p), fftdims(p)))
end

struct BFFTPlanInplace{T, P} <: DHTPlan{T}
    bfftplan::P
    function BFFTPlanInplace(bfftplan::Plan)
        F = eltype(bfftplan)
        P = typeof(bfftplan)
        T = fieldtypes(F)[1]
        return new{T, P}(bfftplan)
    end
end
size(p::BFFTPlanInplace) = size(p.bfftplan)
fftdims(p::BFFTPlanInplace) = fftdims(p.bfftplan)

function LinearAlgebra.mul!(
        y::AbstractArray{T}, p::BFFTPlanInplace{T, P}, x::AbstractArray{T}
    ) where {T, P}
    fx = p.bfftplan * x
    return y .= real.(fx) .+ imag.(fx)
end

function *(p::BFFTPlanInplace{T, P}, x::AbstractArray{T}) where {T, P}
    return LinearAlgebra.mul!(x, p, x)
end

function plan_inv(p::BFFTPlanInplace)
    return ScaledDHTPlan(p, normalization(eltype(p), size(p), fftdims(p)))
end

is_inplace(::BFFTPlanInplace) = true
