module AbstractFastHartleyTransforms
import AbstractFFTs: Plan, eltype, size, fftdims
import Base: size, ndims, length, eltype, *, inv, \, size
using LinearAlgebra: LinearAlgebra

export plan_fht, plan_fht!, plan_ifht, plan_ifht!, fht, fht!, ifht, ifht!, fftdims

include("definitions.jl")
include("bfftplan.jl")
include("r2rplan.jl")

end # module
