# Implementing a backend

A backend supplies a plan subtype, shape/axes accessors, execution, and inverse construction. Generic transform calls then dispatch through `plan_fht(x, dims) * x`. The following deliberately slow, one-dimensional reference backend demonstrates the protocol without an FFT dependency. It allocates a dense matrix and is intended only for explanation and small correctness tests.

```@example reference_backend
using AbstractFastHartleyTransforms, LinearAlgebra
const AFHT = AbstractFastHartleyTransforms

# Store the dense reference operator in the plan for reuse.
struct ReferenceHartleyPlan{T} <: AFHT.DHTPlan{T}
    kernel::Matrix{T}
end

# Expose the input shape and transformed axes through the shared interface.
Base.size(p::ReferenceHartleyPlan) = size(p.kernel)[1:1]
AFHT.fftdims(::ReferenceHartleyPlan) = (1,)
Base.:*(p::ReferenceHartleyPlan{T}, x::Vector{T}) where {T} = p.kernel * x
function LinearAlgebra.mul!(
        y::Vector{T}, p::ReferenceHartleyPlan{T}, x::Vector{T}
    ) where {T}
    return mul!(y, p.kernel, x)
end
function AFHT.plan_inv(p::ReferenceHartleyPlan{T}) where {T}
    # For this symmetric n-point transform, H * H = n * I.
    return AFHT.ScaledDHTPlan(p, inv(T(size(p, 1))))
end

function reference_plan(x::Vector{T}) where {T <: Union{Float32, Float64}}
    n = length(x)
    kernel = Matrix{T}(undef, n, n)
    # Iterate down each column; subtract one for the mathematical indices.
    for j in 1:n, k in 1:n
        phase = T(2) * T(pi) * T((j - 1) * (k - 1)) / T(n)
        kernel[k, j] = cos(phase) + sin(phase)
    end
    return ReferenceHartleyPlan(kernel)
end

x = [1.0, 2.0, 3.0, 4.0]
p = reference_plan(x)
y = p * x
# Verify that the normalized inverse recovers the original input.
@assert inv(p) * y ≈ x
# The output-buffer API should agree with plan multiplication.
out = similar(x)
mul!(out, p, x)
@assert out ≈ y
nothing # hide
```

A production backend should specialize `plan_fht(x::YourArray, dims; kwargs...)` and, if supported, `plan_fht!` for an array type it owns. This example uses a separate constructor to avoid overwriting another package's methods for ordinary vectors. Return a compatible `DHTPlan` and implement `size(p)`, `fftdims(p)`, `*(p, x)`, and `plan_inv(p)`. For preallocated output, add three-argument `LinearAlgebra.mul!`.

Generic `plan_ifht` and `plan_ifht!` scale the corresponding forward plan by the product of the selected axis lengths. Do not scale by the total array size when untransformed batch axes exist. Provide explicit device execution methods instead of relying on a generic fallback that copies into a host array.

## Mutation and differentiation

The internal `is_inplace` trait defaults to false and delegates through scaled plans. The bundled mutating wrappers return true. A new plan with overwritten-input semantics must represent that behavior correctly before using the bundled Enzyme rules. This trait and the qualified plan facilities describe the backend integration protocol.

The CPU Enzyme rules assume a fixed real symmetric Hartley operator and independent derivative buffers. Subtyping `DHTPlan` does not prove that an arbitrary backend meets those assumptions. A new operator, device array, or mutation model needs its own verification and possibly separate differentiation rules. See [Differentiation contract](@ref).

## Checking a backend

Compare small transforms to a direct kernel calculation. Test Float32/Float64, odd/even/singleton lengths, inverse normalization, multidimensional and partial axes, input mutation, and independent batches. Verify `dot(p*x, g) ≈ dot(x, p*g)` for a real symmetric plan, output-buffer behavior, and analytic/finite-difference gradients. Inspect inference and allocations after warmup, and separate plan construction from execution. A dense reference matrix is useful for correctness, not a performance implementation.
