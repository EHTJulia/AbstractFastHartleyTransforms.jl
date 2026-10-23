# Contributing

See [AGENTS.md](AGENTS.md) for coding standards and verification requirements.

## Tests

```sh
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

## Formatting

```sh
julia --project=tools -e 'using Pkg; Pkg.instantiate()'
julia --project=tools tools/format.jl --write
julia --project=tools tools/format.jl
```

JuliaFormatter 2.1.2 is pinned in the tooling environment. The final command checks Julia files and documentation examples without rewriting them.
