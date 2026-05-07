function greet_your_package_name()
	return "Hello HybridDynamics!"
end

#Forward Euler method. We solve an ODE defined by $f(u,t)$ starting at u0 and over tspan with step size dt. 

function forward_euler(f::Function,u0,tspan::Tuple{Float64,Float64},dt::Float64)
    t_start, t_end = tspan

    #Create a time vector:
    t = collect(t_start:dt:t_end) #gets the range of values for time and puts thm into an array (vector). So we can index. 
    num_steps = length(t)

    #intialize the solution array to match the type of initial condition 
    u = Vector{typeof(u0)}(undef, num_steps) #typeof so we can keep things straight. undef is supposedly faster than zeros() but I didnt fact check that. Makes sense though as we skip some values
    u[1] = u0

    for i in 1:(num_steps-1)
        #Eulers update: u_next = u_now + dt * slope
        u[i+1] = u[i] + dt* f(u[i], t[i])
    end
    return t,u
end