# AbstractFastHartleyTransforms
| Documentation | Code Status | Guidelines |
| :------- | :------ | :------ |
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://EHTJulia.github.io/AbstractFastHartleyTransforms.jl/stable/)[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://EHTJulia.github.io/AbstractFastHartleyTransforms.jl/dev/) | [![Build Status](https://github.com/EHTJulia/AbstractFastHartleyTransforms.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/EHTJulia/AbstractFastHartleyTransforms.jl/actions/workflows/CI.yml?query=branch%3Amain)[![Coverage](https://codecov.io/gh/EHTJulia/AbstractFastHartleyTransforms.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/EHTJulia/AbstractFastHartleyTransforms.jl) | [![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac) |


A general framework for fast Hartley transforms (FHTs) in Julia. This framework is following [AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl), providing a popular framework for the fast Fourier Transforms (FFTs) in Julia. Similar to `AbstractFFTs.jl`, **this package is not intended to be used directly**. Instead, packages that implement FHTs extend the types/functions defined in this `AbstractFastHartleyTransforms.jl` package. In this way, multiple FHT packages can use the same interface: for instance, `fht(x)` and `ifht(x)` for the forward/inverse transforms, `plan_fht(x)` and `plan_ifht(x)` for the corresponding plans. Please check the implementations listed below for FHT functions built upon this framework.



## Example implementations of FHTs using this framework
- [`FastHartleyTransform.jl`](https://github.com/EHTJulia/FastHartleyTransform.jl) for CPU and NVIDIA CUDA GPUs.

## Differentiation

Load Enzyme to enable optional forward and reverse rules for fixed Hartley plans on CPU arrays. See the [differentiation guide](docs/src/differentiation.md) for supported operations and limitations.

## Documentation

See the [interface](docs/src/interface.md), [backend implementation guide](docs/src/backends.md), and [API reference](docs/src/api.md).
