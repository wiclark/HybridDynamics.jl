# testing

using HybridDynamics
using LaTeXStrings
# import Pkg; Pkg.add("LaTeXStrings")
import Plots as plt

init = [1.,0.];
tspan = (0., 5.);

function L(x, t)
	q = x[1]
	v = x[2]

	1/2 * v^2 - q^2
end;

Lsys = HybridDynamics.LagSys(L)

probl = HybridDynamics.prob(Lsys, init, tspan)
sl = HybridDynamics.solve(probl, HybridDynamics.ForwardEuler(); dt_initial=0.01)

ql = [x[1] for x in sl[2]]
vl = [x[2] for x in sl[2]]
pl = plt.plot(ql, vl, title = "Harmonic osc. phase portrait via Euler-Lagrange", label = "", xlabel = L"q", ylabel = L"\dot{q}")
display(pl)


## Hamiltonian example phase portrait

# function H(x)
# 	q = x[1]
#     p = x[2]

#     p^2  + q^2
# end;

# 
# probh = HybridDynamics.Prob(H, init, tspan)
# sh = HybridDynamics.solve(probh, HybridDynamics.forward_euler; dt=0.01)

# qh = [x[1] for x in sh[2]]
# vh = [x[2] for x in sh[2]]
# ph = plt.plot(qh, vh, title = "Harmonic osc. phase portrait via Hamilton's", label = "", xlabel = L"q", ylabel = L"\dot{q}")
# display(ph)