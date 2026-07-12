# Change Log

### Version 0.3.0

## External
* Added Documentation! See our dedicated website for the full list of functionality of this package.

* Updated logic of all solvers to utilize the current, predicted and previous (from midpoint) for use in quadratic event detection.

* All systems 'solve' dispatch use prob, solver;... instead of some previously being prob; solver,.... This is done to better match DifferentialEquations.jl.

* Updated 'pluto_examples' to reflect all our changes above. Give it a test!

* Linear and Affine systems now correctly identify if initial conditions start on the guard.

* Linear and Affine systems now also have a toggle for logic order. Flow -> Event check or Event check -> Flow. To utilize this see the arg event_before_flow in the 'solve' dispatch. 

* Updated many many different small tweaks to all the solvers to allow (most of) them to work with all system types. See known issues below for more details. 

# New Functions

* 'extract_jumps'. Used to allow you to customize just the lines connecting discrete jumps within Plotting architecture. (See 'split_jumps' to get rid of them completely!)

## Internal

* Unified variables names across all systems and solvers. Mostly internal but it may have changed external calling. Those changes are below. 

* Added 'take_step_systemtype!' to all system specific solve dispatches. 


## Potential Breaks

Changing some of these things may have affected users code that utilizes the old version. Below are things to watch out for:

* 'jump_times' and 'jump_indices' in Lin/Aff and General systems are now 'event_times' and 'event_indices' to match the rest of the package. 

* As 'solve' dispatches are now prob, solver;... Thus, if you are calling solve in previous versions and write 'HybridDynamics.solve(prob, solver = ...,) you will need to get rid of the 'solver = ' portion as Julia will default to the default solver of the system type instead of use your specified one (usually 'RK45' or 'RK4'). This is a quirk of Julia and it WILL NOT notify you of this default. 

* Let us know of any more we missed!

## Known Bugs

# General Systems

* 'QuadraticLocator' the arg for 'EventLocator' does not currently work for most of the solvers. It will run but the output is not accurate. Currently 'RK45' is the only working solver.  

* 'BisectionLocator' the arg for 'EventLocator' does not currently work for most of the solvers. It will run but the output is not accurate. Currently 'RK45' is the only working solver. 

# Linear/Affine Systems

* Zeno detection is still primitive. Be wary of what it tells you and give us your findings so we can fine tune! Be sure to utilize the tolerances in the 'solve' dispatch to better fit your use case. 

# Mechanical Systems

* Some results will 'fall out' of the accepted region. This is OK! It should only do this for a step before catching itself and getting back on track. 

* 'RK23' struggles to provide accurate data. 

* All event locators (except for Linear!) will not provide accurate data. Some better than others but be careful.

# Nonholonomic Systems

* Some particularly bad solvers (such as Forward/Backwards Euler) will not behave nicely. This is expected but keep it in mind. 

* 'BDF2', 'ModifiedTrap', 'ModifiedMidpoint', 'RK23' and 'RK45' all provide inaccurate data with ALL locator options. 

# Stochastic Systems

* All is in order (for now).

# Filippov Systems

* Bug fixed to work with all adaptive solvers now!

#### Version 0.2.0
Finished implementing all system types. 

#### Version 0.1.0
Initial functionality