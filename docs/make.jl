using Documenter, AbstractIterativeSolvers

makedocs(;
    modules=[AbstractIterativeSolvers],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/tomtrogdon/AbstractIterativeSolvers.jl/blob/{commit}{path}#L{line}",
    sitename="AbstractIterativeSolvers.jl",
    authors="Tom Trogdon",
    assets=String[],
)

deploydocs(;
    repo="github.com/tomtrogdon/AbstractIterativeSolvers.jl",
)
