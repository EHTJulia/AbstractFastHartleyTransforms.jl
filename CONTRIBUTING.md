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

## Documentation

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
DOCS_LOCAL=true julia --project=docs docs/make.jl
```

The local build executes examples and checks references without deployment. Open `docs/build/index.html` to preview it.
