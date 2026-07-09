### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 152874e3-76da-45e7-bd8b-8186b49b851e
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate("..")
    Pkg.resolve()
    Pkg.instantiate()
    Pkg.precompile()

    import HybridDynamics as HD
    import Plots as plt
    using LaTeXStrings
    using PlutoUI
end

# ╔═╡ 71a93050-39ce-4180-a20e-a6b5afe9bb70
# David's solution for an interactable plot
begin
    using PlotlyJS
    plt.plotly()
end

# ╔═╡ 9d28edb2-8291-493a-94c6-f2dac0495e42
@bind ϵ Slider(0.0:0.001:0.01, show_value=true)

# ╔═╡ c63136a5-924f-42de-84fc-3afd247434ab
# Define concrete system parameters, intital condition, and time span
begin
    λ = [1, 1, 1]

    C = [0 1 0
        0 0 1+ϵ
        1-ϵ 0 0]

    A = [0 1 0
        0 0 1
        0 0 0]

    #guard(x) =

    v = C * AbstractVector{Float64}([1.0, 1.0, -2.0])
    time_span = (0.0, 50.0)
end

# ╔═╡ c3c2d095-87c6-42a9-90f6-c5eec78b6ac4
# Construct system and compute solution
begin
    system = HD.LinearSystem(A, λ, C)
    problem = HD.prob(system, v, time_span)
    solution = HD.solve(problem, HD.RK45())

    t_list, x_list = HD.split_jumps(solution)
end

# ╔═╡ f0a61c2e-af58-4b15-8b0e-222f2f63cf54
# Plot solution
begin
    p = plt.plot(title="System Trajectory with v = $(v)",
        xlabel="x",
        ylabel="y",
        zlabel="z",
        #aspect_ratio=:equal,
        grid=true,
        #legend=false
    )

    #plt.surface!(p, guard, label="Guard", color=:blue, alpha=0.3, lw=2)

    plt.plot!(p,
        getindex.(x_list, 1),
        getindex.(x_list, 2),
        getindex.(x_list, 3),
        color=:red,
        #linestyle=:dot,
        lw=3.0,
        label="Trajectory"
    )

    p
end

# ╔═╡ 9a37c68d-2958-41da-bc1d-56dcc18c254b
begin
    diff = diff(solution.jump_times)
    plt.scatter(diff, title="Jump Times", color=:blue)
end

# ╔═╡ Cell order:
# ╠═152874e3-76da-45e7-bd8b-8186b49b851e
# ╠═c63136a5-924f-42de-84fc-3afd247434ab
# ╠═c3c2d095-87c6-42a9-90f6-c5eec78b6ac4
# ╠═71a93050-39ce-4180-a20e-a6b5afe9bb70
# ╠═9d28edb2-8291-493a-94c6-f2dac0495e42
# ╠═f0a61c2e-af58-4b15-8b0e-222f2f63cf54
# ╠═9a37c68d-2958-41da-bc1d-56dcc18c254b
