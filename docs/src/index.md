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

For demonstrations, see the demos folder for Jupyter and Pluto notebooks:
- "examples_pluto.jl": A Pluto notebook providing an example of each problem type.
- "bouncing_ball_pluto.jl": A Pluto notebook covering the dissipative bouncing ball as different problems.
- "life_after_zeno_pluto.jl": A Pluto notebook showcasing various mechanical examples with Zeno trajectories.
- "FilippovFuller_jupyter.ipynb": A Jupyter notebook presenting Fuller's problem as an example of a chattering Filippov system.
- "FilippovSSF_jupyter.ipynb": A Jupyter notebook with modeling a system with friction as a Filippov system.
- "TestingNotebook_jupyter.ipynb": A Jupyter notebook that performs error analysis for various methods.

The demonstrations require the pacakges: Plots, LaTeXStrings, Printf (only 'TestingNotebook_jupyter.ipynb'), and PlutoUI (all Pluto notebooks).

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

- **NonholonomicSystem** 

    Required inputs:
    - ``M``: Mass matrix
    - ``V``: Potential energy

    Optional inputs:
    - ``A``: The Pfaffian constraint matrix
    - ``h``: Guard
    - ``Δ``: Reset map
    - ``e``: Coefficient of restitution

- **StochasticSystem** 

    Required inputs:
    - ``f``: The drift dynamics
    - ``g``: The diffusion term
    - ``h``: Guard
    - ``Δ``: Reset map

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

## Acknowledgments
This work was supported by AFOSR grant FA9550-32-0400.