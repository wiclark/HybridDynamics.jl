using Documenter
using HybridDynamics

# Automatically populate index.md with the git readme
write(joinpath(@__DIR__, "src", "index.md"),
      read(joinpath(".", "README.md"), String))

# Make the build
makedocs(
    sitename = "HybridDynamics",
    format = Documenter.HTML(mathengine = Documenter.MathJax3()),  # The mathengine should already be MathJax
    modules = [HybridDynamics],
    checkdocs = :none,              # Remove when...documentation is written?
    # remotes = nothing,               # Remove when...the git repo is public?
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Change Log" => "changelog.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
