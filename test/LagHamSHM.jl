# testing

using HybridDynamics
using LaTeXStrings
import Plots as plt

init = [1.,0.];
tspan = (0., 5.);

function L(x, t)
	q = x[1]
	v = x[2]

	1/2 * v^2 - q^2
end;

Lsys = HybridDynamics.LagSys(L)
probl = HybridDynamics.Prob(Lsys, init, tspan)
sl = HybridDynamics.solve(probl, HybridDynamics.forward_euler; dt=0.01)

ql = [x[1] for x in sl[2]]
vl = [x[2] for x in sl[2]]
pl = plt.plot(ql, vl, title = "Harmonic osc. phase portrait via Euler-Lagrange", label = "", xlabel = L"q", ylabel = L"\dot{q}")
display(pl)


function H(x)
	q = x[1]
    p = x[2]

    p^2  + q^2
end;

# The Hamiltonain has a nice clean functional representation so use the automatic differentiating struct ADHamiltonian()
probh = HybridDynamics.Prob(HybridDynamics.ADHamiltonian(H), init, tspan)
sh = HybridDynamics.solve(probh, HybridDynamics.forward_euler; dt=0.01)

qh = [x[1] for x in sh[2]]
vh = [x[2] for x in sh[2]]
ph = plt.plot(qh, vh, title = "Harmonic osc. phase portrait via Hamilton's", label = "", xlabel = L"q", ylabel = L"\dot{q}")
display(ph)