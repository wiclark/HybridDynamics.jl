# testing

using HybridDynamics
using LaTeXStrings
# import Pkg; Pkg.add("LaTeXStrings")
import Plots as plt

init = [1.,0.];
tspan = (0., 5.);

function M(q)
    [1.0]
end

function V(q)
    q[1]^2
end

# Table/guard
g(q) = q[1]

Lsys = HybridDynamics.LagSys(M, V; guard = g)

probl = HybridDynamics.prob(Lsys, init, tspan)
sl = HybridDynamics.solve(probl, HybridDynamics.ForwardEuler(); dt_initial=0.01)

ql = [x[1] for x in sl.X]
vl = [x[2] for x in sl.X]
pl1 = plt.plot(ql, vl, title = "Harmonic osc. phase portrait via Euler-Lagrange", label = "", xlabel = L"q", ylabel = L"\dot{q}")
display(pl1)
pl2 = plt.plot(sl.T, ql, title = "1-D SHM", xlabel = L"t", ylabel = L"q", label = "")
display(pl2)


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