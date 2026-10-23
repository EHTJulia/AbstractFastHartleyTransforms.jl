using AbstractFastHartleyTransforms
using Documenter

DocMeta.setdocmeta!(
    AbstractFastHartleyTransforms,
    :DocTestSetup,
    :(using AbstractFastHartleyTransforms);
    recursive=true,
)

local_build = get(ENV, "DOCS_LOCAL", "false") == "true"

makedocs(;
    remotes=local_build ? nothing : Dict(),
    modules=[AbstractFastHartleyTransforms],
    authors="Kazunori Akiyama",
    sitename="AbstractFastHartleyTransforms.jl",
    format=Documenter.HTML(;
        canonical="https://EHTJulia.github.io/AbstractFastHartleyTransforms.jl",
        edit_link=local_build ? nothing : "main",
        repolink="https://github.com/EHTJulia/AbstractFastHartleyTransforms.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Interface and conventions" => "interface.md",
        "Implementing a backend" => "backends.md",
        "Differentiation" => "differentiation.md",
        "API reference" => "api.md",
    ],
)

if !local_build
    deploydocs(;
        repo="github.com/EHTJulia/AbstractFastHartleyTransforms.jl", devbranch="main"
    )
end
