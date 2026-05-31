
#This file is to establish the foundation of everything we use with Multiple Dispatch. It doesnt need to be here exactly but I think it will be a good place to keep it organized. 

#Instead of writing a lot of if/else statements in our solveloop we define these abstract types.
#The methods at the bottom are used to define parts of new systems we might add so that we can still use the architecture here. 
#If we want to add a new system we can simply declare it as a subtype of AbstractHybridSystem and define the four functions at the bottom. The solver then will automatically know how to simulate it.

#---------------------------
#ABSTRACT TYPES
#These are the empty categories that hold nothing themselves. By restricting our solver function args to these types, we can make sure the solver can accept ANY system or solver that is in these families.

#Parent category for any continuous ODE integration method
abstract type AbstractODESolver end

#Parent Category for any root-finding algorithm used to pinpoint the exact time of a guard crossing.
abstract type AbstractEventLocator end

#Parent Category for any physical system that features both continuous flow and discrete jumps (unsure how Filippov will integrate but I hope its not too bad)
abstract type AbstractHybridSystem end

#Parent Category for a simulation scenario which takes the system, initial conditions, and timespan 
abstract type AbstractHybridProblem end

#-----------------
#SOLVE TAGS
#These are empty structs. Their purpose is to act as a tag for the Julia compiler. When a user passes 'ModifiedMidpoint()' to the solver, the compiler uses multiple dispatch to route to the 'take_step' function written for that tag. 

#Tag for standard Forward Euler integration 
struct ForwardEuler <: AbstractODESolver end

#Tag for the Modifed Trapezoidal method 
struct ModifiedTrap <: AbstractODESolver end

#Tag for the Modifed Midpoint method 
struct ModifiedMidpoint <: AbstractODESolver end

#Tag for exact integration of Linear systems via the matrix Exponential 
struct ExponentialSolver <: AbstractODESolver end

#Extrapolation struct tags

#Richardson Extrapolation tag
struct RichardsonExtrapolation <: AbstractODESolver end

#-----------------------
#LOCATOR TAGS
#These are similar to the solver tags. But these define how the solver finds the exact crossing time once an event is detected. 

#Tag to use Linear Interpolation. Very fast but can be innacurate for higher order methods. 
struct LinearLocator <: AbstractEventLocator end

#Tag to use a bisection method serach. Can be very accurate but also very slow with complex systems
struct BisectionLocator <: AbstractEventLocator end


#---------------------
#SYSTEM INTERFACE
#These four (later 6?) functions represent the minimum requirements to simulate a system. This is the part I think is most subject to change. I added it here to try to fully future proof things but only time will tell. 

#We set them to throw errors here as a safety. If one of us defines a new system but forgets to add these the code will give the error which is very easy to pinpoint, rather than crashing and being annoying.

#Must return the number of state variables (length of state vector 'x')
get_dimension(sys::AbstractHybridSystem) = error("get_dimension not implemented for $(typeof(sys))")

#Must return a function of the form (x,t) -> dx/dt representing the continuous dynamics.
vector_field(sys::AbstractHybridSystem) = error("vector_field not implemented for $(typeof(sys))")

#Must evaluate the guard surface at state 'x'. A sign change will trigger the jump detector
guard(sys::AbstractHybridSystem, x) = error("guard not implemented for $(typeof(sys))")

#Must apply discrete jump to the pre-jump state 'x' to return the post-jump state. 
apply_reset(sys::AbstractHybridSystem, x) = error("apply_reset not implemented for $(typeof(sys))")

#OPTIONAL: ONLY USE IF YOUR SYSTEM EXHIBITS BEATING/BLOCKING/ZENO 
#determins if the system is stuck beating or maybe Zeno?
#Default behavior for any system is to ignore this as its mainly for Linear/Affine stuffs. 
check_beating_status(sys::AbstractHybridSystem, instant_jumps, n, x_star, t_star, tol) = :continue

#This will help us figure out Zeno I think but I have yet to iron that out. 
#process_jump_event(sys::AbstractHybridSystem, x_minus, t_star, instant_jumps, tol) = :continue

#-----------------------------------------------
#This may be a little strange to implement with something like a Filippov system but here are my thoughts. Instead of a smooth function we would have the
# vector_field condition return a function that has if/else logic depending on the side of the state we are on. 
#The guard would still work I believe. And the reset does not have us 'jump' so the apply_reset would just be 'return x' or smth. 
#The caveat is the sliding issue. I am unsure how to solve that right now. We might just need to find a way to detect it and then if we do send it to an entirely new solver tag (e.g. struct SlidingModSolver <: AbstractODESolver end) to make it work. 

#This method should also be fine with regards to having no guard in our system and solving it normally like an ODE solver would. We would just return guard(sys::ContinuousSystem, x) = 1.0
#as when it takes a step it will always return 1 and the sign never changes to indicate a jump. Of course this doesnt apply to jumps where the sign doesnt change but I havent figured that out yet. 