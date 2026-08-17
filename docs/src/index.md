# HybridDynamics.jl

[![Build Status](https://github.com/wiclark/HybridDynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/wiclark/HybridDynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package for studying the dynamics of hybrid systems.

## Documentation
Documentation is available [here](https://wiclark.github.io/HybridDynamics.jl/dev/)

## Installation

To install this package, run:
```julia
using Pkg
Pkg.add("HybridDynamics")
```
in a Julia REPL.

The package can then be loaded using the command:
```julia
import HybridDynamics as HD
```

## Usage

This package supports the following system types: "GeneralSystem", "MechanicalSystem", "NonholonomicSystem", "StochasticSystem", "FilippovSystem", and "Linear/AffineSystem".

For a demonstration, see the Pluto notebooks: "bouncing_ball.jl" and "pluto_examples.jl".

## Systems

- **GeneralSystem** ($f, h, \Delta$)

    Required inputs:
    - ``f``: Continuous dynamics ($f = dx/dt$)
    - ``h``: Guard
    - ``reset``: Reset map

- **MechanicalSystem** ($M, V$)

    Required inputs:
    - ``M``: Mass matrix
    - ``V``: Potential energy

    Optional inputs:
    - ``h``: Guard
    - ``Δ``: Reset map
    - ``e``: Coefficient of restitution

- **LinearSystem** ($A, \lambda, C$)

    Required inputs:
    - ``A``: State transition matrix
    - ``λ``: Normal vector to the guard
    - ``C``: Reset map matrix ($x⁺ = Cx$)

- **AffineSystem** ($A, b, \lambda, a, C, \kappa$)

    Required inputs:
    - ``A``: State transition matrix
    - ``b``: Continuous affine vector ($dx/dt = Ax + b$)
    - ``λ``: Normal vector to the guard
    - ``a``: Guard Offset const $dot(λ, x) + a = 0$
    - ``C``: Reset map matrix ($x⁺ = Cx$)
    - ``κ``: Discrete affine vector const x⁺ = Cx + κ

- **FilippovSystem** ($F, G, H$)

    Required inputs:
    - ``F``: Continuous dynamics where $H(x) > 0$
    - ``G``: Continuous dynamics where $H(x) < 0$
    - ``H``: Guard

    Optional inputs:
    - ``N``: Normal to the guard
