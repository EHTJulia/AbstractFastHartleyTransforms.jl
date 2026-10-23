# Differentiation contract

Loading Enzyme activates the optional EnzymeCore extension. Enzyme is not a mandatory runtime dependency. The verified baseline is Julia 1.10.12, Enzyme 0.13.216, and EnzymeCore 0.8.26.

## Supported CPU behavior

Single-lane forward and reverse rules cover dense Float32/Float64 host arrays with fixed Hartley plans, including scaled inverse plans, mutating plans, and three-argument `LinearAlgebra.mul!`. The concrete FFTW implementation and its convenience calls are tested through FastHartleyTransform.

For a symmetric forward operator `H`, the tangent is `H*dx` and the reverse contribution is `H*g`. For `H/M`, both include the inverse scale. Nonmutating reverse calls accumulate into the input cotangent. Mutating reverse calls replace the cotangent of the overwritten value; output-buffer rules account for the overwritten output and exact input/output aliasing.

## Activity and plans

  - Build plans outside differentiated workloads and pass them as `Enzyme.Const`.
  - Keep differentiable numerical parameters in separate active arrays. Plan metadata, dimensions, flags, and stored scalar coefficients are fixed for these rules.
  - Provide independent, correctly shaped primal and derivative buffers via `Duplicated`.
  - Enable runtime activity for mixed active/constant broadcasts when required, as in the integration examples.
  - Restore direction vectors and primal inputs when a mutating operation overwrites them.

With runtime activity enabled, an array that dynamically becomes constant may have an aliased primal and shadow. Transform rules explicitly reject that situation to avoid silently producing an incorrect gradient; compute constant transforms outside AD.

Batched AD lanes, higher-order AD, views, and CUDA derivatives are not validated. Independent samples stored on untransformed axes of an ordinary array are supported and are distinct from batched AD lanes.

## Reactant

Reactant support belongs to the concrete FastHartleyTransform extension. It lowers FFT operations to a compiled representation and differentiates that compiled workload. Native EnzymeCore rules and Reactant differentiation are separate mechanisms. See the concrete package's [backend guide](https://ehtjulia.github.io/FastHartleyTransform.jl/dev/backends/) for device support.

See the concrete package's [differentiation guide](https://EHTJulia.github.io/FastHartleyTransform.jl/dev/differentiation/) for Enzyme usage and Reactant compilation examples.
