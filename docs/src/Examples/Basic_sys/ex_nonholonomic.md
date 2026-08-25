## A Nonholonomic Example
The Chaplygin sleigh is a nonholonomic system on $(x,y,\theta) \in \mathbb{R}^2\times \mathbb{S}^1$ with mass matrix
```math
	M(q) = \begin{bmatrix}
		m & 0 & -ma\sin\theta \\
		0 & m & ma\cos\theta \\
		-ma\sin\theta & ma\cos\theta & I+ma^2
	\end{bmatrix}
```
potential energy $V = 0$, and constraint matrix
```math
	A(q) = \begin{bmatrix}
		-\sin\theta & \cos\theta & 0
	\end{bmatrix}.
```
Creating a 'nonholonomic billiard' with a circular table makes the event function
```math
	h(x,y) = 1 - (x^2+y^2) \geq 0.
```

```@example noholo
import HybridDynamics as HD
import Plots as plt
using LaTeXStrings

a = 0.2

M_sleigh(q) = [1.0 0.0 -a*sin(q[3]);
        0.0 1.0 a*cos(q[3]);
        -a*sin(q[3]) a*cos(q[3]) 1+a^2]
V_sleigh(q) = 0.0
A_sleigh(q) = [-sin(q[3]) cos(q[3]) 0.0]
h_sleigh(q) = 1 - (q[1]^2+q[2]^2)
∇h_sleigh(q) = [-2*q[1], -2*q[2], 0.0]

sysNH = HD.NonholonomicSystem(M_sleigh, V_sleigh; A = A_sleigh, guard=h_sleigh, normal=∇h_sleigh, e=1)
probNH = HD.prob(sysNH, [0.0, 0.0, 0.0, 1.0, 0.0, 1.5], (0.0, 10.0))
solNH = HD.solve(probNH, solver=HD.RK4())
```

```@example noholo
θ = LinRange(0, 2π, 100)
xs = getindex.(solNH.x, 1)
ys = getindex.(solNH.x, 2)
plt.plot(xs, ys, label="Sleigh Trajectory", lw=2)
plt.plot!(cos.(θ), sin.(θ), label="", lc=:black, aspect_ratio=1)
```