### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# ╔═╡ 7a6c05e0-74ad-11f1-a8e2-5b2375a4e38b
begin
	using Pkg
	Pkg.activate("..")
	Pkg.resolve()
	Pkg.instantiate()
	Pkg.precompile()
end

# ╔═╡ c56c396c-ad27-49b9-938e-ae946c05ca57
begin
	import HybridDynamics as HD
	import Plots as plt
	using LaTeXStrings
	using PlutoUI
end

# ╔═╡ 2dcc2a4a-4c80-46d9-a2c4-affd0c580d21
begin
	f(x, t) = [1.0]
	g(x, t) = [0.1;;]
	h(x) = x[1]-1.0
	Δ(x) = x .-1.0
end

# ╔═╡ c38e9792-3ae7-4355-b5c0-350f7389b75e
length(f(1.0, 1.0))

# ╔═╡ c584bf14-945c-4faf-8983-e32e2a6c959c
size(g(1.0, 1.0))[1]

# ╔═╡ 101e8feb-588b-4a01-8e4a-8267e7603f9a
begin
	sysST = HD.StochasticSystem(f, g, h, Δ)
	probST = HD.prob(sysST, [0.5], (0.0, 3.0))
	solST = HD.solve(probST)
end

# ╔═╡ 9f30bd50-f027-4a5f-9b55-03bd2898d44e
xd = getindex.(solST.x, 1)

# ╔═╡ 8a1237c0-e6a7-4afe-b690-ef8dbb96dc20
plt.plot(solST.t, xd)

# ╔═╡ 9aa5731d-3e45-4f5a-9f6d-64bbb87a6019
size([1.0;;])

# ╔═╡ Cell order:
# ╠═7a6c05e0-74ad-11f1-a8e2-5b2375a4e38b
# ╠═c56c396c-ad27-49b9-938e-ae946c05ca57
# ╠═2dcc2a4a-4c80-46d9-a2c4-affd0c580d21
# ╠═c38e9792-3ae7-4355-b5c0-350f7389b75e
# ╠═c584bf14-945c-4faf-8983-e32e2a6c959c
# ╠═101e8feb-588b-4a01-8e4a-8267e7603f9a
# ╠═8a1237c0-e6a7-4afe-b690-ef8dbb96dc20
# ╠═9f30bd50-f027-4a5f-9b55-03bd2898d44e
# ╠═9aa5731d-3e45-4f5a-9f6d-64bbb87a6019
