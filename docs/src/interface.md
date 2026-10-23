# Interface and conventions

## Hartley transform

With zero-based mathematical indices on the transformed axes,

```math
(Hx)_{\mathbf{k}}=\sum_{\mathbf{n}}x_{\mathbf{n}}
\left[\cos\!\left(2\pi\sum_d n_d k_d/N_d\right)
+\sin\!\left(2\pi\sum_d n_d k_d/N_d\right)\right].
```

Julia indices start at one. The `dims` argument selects transformed axes; other axes are independent batches. Let `M` be the product of lengths of the selected axes. The default transform is real and unnormalized, with `H' = H`, `H*H = M*I`, and inverse `H/M`. Therefore the generic inverse plan scales the forward plan by `1/M`.

A product of one-dimensional Hartley transforms differs from the single-kernel multidimensional transform. FFTW's explicit `use_r2r=true` option in the concrete CPU backend preserves that separable convention. Both conventions have symmetric inverses with the same scale, but their forward values generally differ.

## Calls and plans

| Interface                          | Contract                                         |
|:---------------------------------- |:------------------------------------------------ |
| `fht(x, dims)` / `ifht(x, dims)`   | Forward / normalized inverse; preserve input     |
| `fht!(x, dims)` / `ifht!(x, dims)` | Overwrite input                                  |
| `plan_fht(x, dims; kwargs...)`     | Reusable forward plan                            |
| `plan_fht!(x, dims; kwargs...)`    | Reusable forward plan with input mutation        |
| `plan_ifht` / `plan_ifht!`         | Generic scaled forward plans                     |
| `p * x`                            | Apply a plan with its declared mutation behavior |
| `inv(p)` / `p \ x`                 | Invert / apply the inverse plan                  |
| `mul!(out, p, x)`                  | Execute into an existing output when implemented |

Omitting `dims` selects all axes. The generic convenience transform methods do not accept planner keywords: create a plan to choose backend options. A plan is tied to compatible element type, shape, axes, and backend layout. Input conversion and offset-axis handling can copy data; do not infer allocation-free execution from the shared interface.

## FFT wrappers

Qualified names `AbstractFastHartleyTransforms.BFFTPlan` and `BFFTPlanInplace` wrap an unnormalized backward FFT and extract `real + imag`. `R2RPlan` and `R2RPlanInplace` wrap a real-to-real DHT plan. They delegate shape and axes to the wrapped AbstractFFTs-compatible plan and implement three-argument `mul!`.

The mutating wrappers return the transformed input for multiplication. With separate input/output buffers, `mul!` writes to the output. Exact same-buffer operation is tested; partially overlapping views are outside the verified contract. FFT workspace allocation can remain even when the real output is preallocated.

`ScaledDHTPlan` represents a scalar multiple of a plan and supplies generic inverse scaling. These qualified plan types are backend integration facilities, not all exported user-facing names. Review the implementation and tests before extending the protocol for a new device.
