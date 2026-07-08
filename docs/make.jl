using Documenter
using HybridDynamics

makedocs(
    sitename = "HybridDynamics",
    format = Documenter.HTML(),
    modules = [HybridDynamics],
    # checkdocs = :none,              # Remove when...documentation is written?
    # remotes = nothing,               # Remove when...the git repo is public?
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
