using ARPESMakie
using Documenter

DocMeta.setdocmeta!(ARPESMakie, :DocTestSetup, :(using ARPESMakie); recursive = true)

makedocs(;
    modules = [ARPESMakie],
    authors = "Ryuichi Arafune",
    sitename = "ARPESMakie.jl",
    format = Documenter.HTML(;
        canonical = "https://arafune.github.io/ARPESMakie.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(;
    repo = "github.com/arafune/ARPESMakie.jl.git",
    devbranch = "main",
    push_preview = true,
)
