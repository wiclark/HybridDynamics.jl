# HybridDynamics.jl

Documentation for HybridDynamics.jl

# HybridDynamics

[![Build Status](https://github.com/wiclark/HybridDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wiclark/HybridDynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package for studying the dynamics of hybrid systems.

## Usage
This package supports the following system types: "GeneralSystem", "MechanicalSystem", "NonholonomicSystem", "StochasticSystem", "FilippovSystem", and "Linear/AffineSystem".

For a demonstration, see the Pluto notebooks: "bouncing_ball.jl" and "pluto_examples.jl".

## Systems

- **GeneralSystem** ($f, h, \Delta$)

    Required inputs:
    - $f$: Continuous dynamics ($f = dx/dt$)
    - $h$: Guard
    - $\Delta$: Reset map
<br><br>

- **MechanicalSystem** ($M, V$)

    Required inputs:
    - $M$: Mass matrix
    - $V$: Potential energy

    Optional inputs:
    - $h$: Guard
    - $\Delta$: Reset map
    - $e$: Coefficient of restitution
<br><br>

- **NonholonomicSystem** 

    Required inputs:
    - $M$: Mass matrix
    - $V$: Potential energy

    Optional inputs:
    - $A$: The Pfaffian constraint matrix
    - $h$: Guard
    - $\Delta$: Reset map
    - $e$: Coefficient of restitution

- **StochasticSystem** 

    Required inputs:
    - $f$: The drift dynamics
    - $g$: The diffusion term
    - $h$: Guard
    - $\Delta$: Reset map

- **LinearSystem** ($A, \lambda, C$)

    Required inputs:
    - $A$: State transition matrix
    - $\lambda$: Normal vector to the guard
    - $C$: Reset map matrix ($x⁺ = Cx$)
<br><br>

- **AffineSystem** ($A, b, \lambda, a, C, \kappa$)

    Required inputs:
    - $A$: State transition matrix
    - $b$: Continuous affine vector ($dx/dt = Ax + b$)
    - $\lambda$: Normal vector to the guard
    - $a$: Guard Offset const $dot(λ, x) + a = 0$
    - $C$: Reset map matrix ($x⁺ = Cx$)
    - $\kappa$: Discrete affine vector const x⁺ = Cx + κ
<br><br>

- **FilippovSystem** ($F, G, H$)

    Required inputs:
    - $F$: Continuous dynamics where $H(x) > 0$
    - $G$: Continuous dynamics where $H(x) < 0$
    - $H$: Guard

    Optional inputs:
    - $N$: Normal to the guard