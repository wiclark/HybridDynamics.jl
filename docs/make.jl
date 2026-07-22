using Documenter
using HybridDynamics

# Automatically populate index.md with the git readme
write(joinpath(@__DIR__, "src", "index.md"),
      read(joinpath(".", "README.md"), String))

# Make the build
makedocs(
    sitename = "HybridDynamics",
    format = Documenter.HTML(), 
    modules = [HybridDynamics],
    checkdocs = :none,              # Remove when...documentation is written?
    # remotes = nothing,               # Remove when...the git repo is public?
    pages = [
        "Home" => "index.md",
        "System Types" => Any[
            "General" => "System_Types/General.md",
            "Linear and Affine" => "System_Types/Linear_Affine.md",
            "Mechanical" => "System_Types/Mechanical.md",
            "Nonholonomic" => "System_Types/Nonholonomic.md",
            "Stochastic" => "System_Types/Stochastic.md",
            "Filippov" => "System_Types/Filippov.md"
        ],
        "Solver Algorithms" => "Solver_Algorithms/all.md",
        "Examples" => Any[
            "General" => "Examples/ex_general.md",
            "Linear and Affine" => "Examples/ex_linear_affine.md",
            "Mechanical" => "Examples/ex_mechanical.md",
            "Nonholonomic" => "Examples/ex_nonholonomic.md",
            "Stochastic" => "Examples/ex_stochastic.md",
            "Filippov" => "Examples/ex_filippov.md"
        ],
        "Change Log" => "changelog.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/wiclark/HybridDynamics.jl.git",
    devbranch = "main"
)
