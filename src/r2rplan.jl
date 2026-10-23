##############################################################################
# This is a generic FHT Plan using R2R transform

struct R2RPlan{T,P} <: DHTPlan{T}
    r2rplan::P
    function R2RPlan(r2rplan::Plan)
        T = eltype(r2rplan)
        P = typeof(r2rplan)
        return new{T,P}(r2rplan)
    end
end
size(p::R2RPlan) = size(p.r2rplan)
fftdims(p::R2RPlan) = fftdims(p.r2rplan)

function LinearAlgebra.mul!(
    y::AbstractArray{T}, p::R2RPlan{T}, x::AbstractArray{T}
) where {T}
    return LinearAlgebra.mul!(y, p.r2rplan, x)
end

function *(p::R2RPlan{T,P}, x::AbstractArray{T}) where {T,P}
    return p.r2rplan * x
end

function plan_inv(p::R2RPlan)
    return ScaledDHTPlan(p, normalization(eltype(p), size(p), fftdims(p)))
end

struct R2RPlanInplace{T,P} <: DHTPlan{T}
    r2r!plan::P
    function R2RPlanInplace(r2rplan::Plan)
        T = eltype(r2rplan)
        P = typeof(r2rplan)
        return new{T,P}(r2rplan)
    end
end
size(p::R2RPlanInplace) = size(p.r2r!plan)
fftdims(p::R2RPlanInplace) = fftdims(p.r2r!plan)

function LinearAlgebra.mul!(
    y::AbstractArray{T}, p::R2RPlanInplace{T,P}, x::AbstractArray{T}
) where {T,P}
    y === x || copyto!(y, x)
    return LinearAlgebra.mul!(y, p.r2r!plan, y)
end

function *(p::R2RPlanInplace{T,P}, x::AbstractArray{T}) where {T,P}
    return LinearAlgebra.mul!(x, p, x)
end

function plan_inv(p::R2RPlanInplace)
    return ScaledDHTPlan(p, normalization(eltype(p), size(p), fftdims(p)))
end

is_inplace(::R2RPlanInplace) = true
