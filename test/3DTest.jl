using HybridDynamics
using Plots

use_affine = true #CHANGE THIS FOR LINEAR VS AFFINE TEST

A = [0.0 -4.0 0.0;
    4.0 -0.2 -1.5;
    0.0 1.5 -0.4]

λ = [1.0, 2.0, 0.0]

C = [1.30 0.20 0.00;
     0.00 0.80 0.00;
     0.00 0.00 0.90]

x0 = [1.0, .5, -.5]

tspan = (0.0, 10.0)

if use_affine
    b = [1.5, 0.0, .5]
    a = .2
    κ = [.1, -.2, .3]

    sys = AffineSystem(A, b, λ, a, C, κ)
    problem = prob(sys, x0, tspan)
    
    solver = RK45()
    sol = solve(problem, solver)

    println("You just used $(typeof(solver))")

    _, X = split_jumps(sol)

    p = plot(title = "3D Affine System Test", xlabel = "x1", ylabel = "x2", zlabel = "x3", aspect_ratio = :equal, legend = false)
    y = range(-5, 5, length=25)
    z = range(-5, 5, length=25)

    Y = repeat(y', length(z), 1)
    Z = repeat(z, 1, length(y))

    X_guard = fill(-a, size(Y))
    surface!(p,X_guard,Y,Z,alpha = 0.25,label = "Guard")
    plot!(p, getindex.(X,1), getindex.(X,2), getindex.(X,3), lw=2)
    scatter!(p,[x0[1]],[x0[2]],[x0[3]],markersize = 5,label = "Initial Condition")

    display(p)

else
    sys = LinearSystem(A, λ, C)
    problem = prob(sys, x0, tspan)

    solver = RK45()
    sol = solve(problem, solver)

    println("You just used $(typeof(solver))")

    _, X = split_jumps(sol)

    p = plot(title = "3D Linear System", xlabel = "x1",ylabel = "x2",zlabel = "x3",aspect_ratio = :equal,legend = false)
    y = range(-5, 5, length=25)
    z = range(-5, 5, length=25)

    Y = repeat(y', length(z), 1)
    Z = repeat(z, 1, length(y))

    X_guard = fill(0.0, size(Y))

    surface!(p,X_guard,Y,Z,alpha = 0.25,label = "Guard")

    plot!(p,getindex.(X,1),getindex.(X,2),getindex.(X,3),lw = 2)
    scatter!(p,[x0[1]],[x0[2]],[x0[3]],markersize = 5,label = "Initial Condition")

    display(p)
end