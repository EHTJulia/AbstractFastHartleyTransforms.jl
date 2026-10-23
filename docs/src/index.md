# AbstractFastHartleyTransforms.jl

A common Julia interface for real Hartley transforms. This package defines transform functions, plan types, inverse normalization, FFT-backed wrappers, and optional EnzymeCore rules. Concrete backends supply planning methods and execution.

For everyday transforms, install [FastHartleyTransform.jl](https://EHTJulia.github.io/FastHartleyTransform.jl/dev/), which implements the interface using FFTW and optional CUDA/Reactant paths. Loading the abstract package alone does not provide a general CPU planner.

## Read next

  - [Interface and conventions](@ref): mathematical definition, axes, normalization, and mutation.
  - [Implementing a backend](@ref): a complete small reference backend and the extension protocol.
  - [Differentiation contract](@ref): optional CPU Enzyme rules and their limitations.
  - [API reference](@ref): shared transform and planning functions.

Julia 1.10 or later is required. The Enzyme behavior described here requires AbstractFastHartleyTransforms 0.1.3 or later.
