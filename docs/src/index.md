# InitialValueProblems.jl Documentation

```@contents
```

## Introduction

InitialValueProblems.jl is a light-weight package for solving ordinary differential equations (ODEs) with single-step integrators. It provides:

- Abstract types for defining ODEs
- Single-step numerical integration schemes (explicit and implicit)
- Analytical solution propagators for select ODE types
- A small catalog of model problems
- Plots.jl integration for easy visualization

## Installation

Since this package is not yet registered, install it directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/vlc1/InitialValueProblems.jl")
```

## Quick Start

Here's a simple example solving a sinusoidal ODE and comparing numerical methods against the analytical solution:

```julia
using Plots
using InitialValueProblems

# Define a sinusoidal ODE: dy/dx = A*sin(ω*x) - λ*y
eq = Sinusoidal(lambda = 0.2, omega = 4.0, amp = [0.5])

# Initial condition and time range
y0 = [1.0]
xs = range(0.0, 10.0, length=100)

# Plot analytical solution
prop = Propagator(eq)
plot(prop, xs, y0, 1, label="Analytical", linewidth=2)

# Compare with numerical integrators
tau = step(xs)
plot!(ForwardEuler(eq), xs, copy(y0), 1, (), label="Forward Euler")
plot!(BackwardEuler(eq), xs, copy(y0), 1, ([0.0],), label="Backward Euler")
plot!(Midpoint(eq), xs, copy(y0), 1, ([0.0],), label="Midpoint")

xlabel!("x")
ylabel!("y")
```

## Core Types

```@docs
OrdinaryDifferentialEquation
Integrator
Propagator
```

## Numerical Integrators

The package provides single-step integrators for first-order ODEs:

```@docs
ForwardEuler
BackwardEuler
Midpoint
```

### Usage

All integrators are callable objects that update the state in-place:

```julia
# Explicit method (no cache needed)
scheme = ForwardEuler(eq)
scheme(x, y, tau)  # Updates y in-place

# Implicit methods (require cache for nonlinear solver)
cache = (similar(y),)
scheme = BackwardEuler(eq)
scheme(x, y, tau, cache)  # Updates y in-place
```

## Model Catalog

```@docs
Sinusoidal
```

## Analytical Solutions

For ODEs with known analytical solutions, the `Propagator` type provides exact solution evaluation:

```@docs
propagate
```

### Example

```julia
eq = Sinusoidal(lambda = 0.2, omega = 4.0, amp = [0.5])
prop = Propagator(eq)

x0, y0 = 0.0, [1.0]
tau = 1.0
i = 1

# Compute analytical solution at x0 + tau
y_exact = prop(x0, y0, tau, i)
```

## Plots Integration

The package provides Plots.jl recipes for both `Propagator` and `Integrator` types, enabling direct plotting:

```julia
using Plots

# Plot analytical solution
plot(propagator, time_range, initial_condition, component_index)

# Plot numerical solution
plot(integrator, time_range, initial_condition, component_index, cache)
```

The recipes automatically evaluate the solution over the provided time range and return plottable data.

## Implementing Custom ODEs

To implement a custom ODE, subtype `OrdinaryDifferentialEquation{N}` where `N` is the order. Your ODE type must implement two call signatures:

### 1. Explicit Integrator Signature

Required for explicit methods like `ForwardEuler`:

```julia
function (eq::MyODE)(x, y::AbstractArray, α::Number)
    # Overwrite y with f(x, y) * α + y
    # Return y
end
```

This signature computes `y ← y + α·f(x, y)` in-place.

### 2. Implicit Integrator Signature

Required for implicit methods like `BackwardEuler` and `Midpoint`:

```julia
function (eq::MyODE)(res::AbstractArray, x::Number, y::AbstractArray,
                      inc::AbstractArray, α::Number, β::Number=one(α))
    # Overwrite res with f(x, y + β * inc) * α + inc
    # Return res
end
```

This signature computes `res ← inc + α·f(x, y + β·inc)` in-place. The optional parameter `β` defaults to `one(α)`:
- For `BackwardEuler`: `β = 1` (evaluates `f(x, y + inc)`)
- For `Midpoint`: `β = 0.5` (evaluates `f(x, y + inc/2)`)

### Complete Example

```julia
struct MyODE{T} <: OrdinaryDifferentialEquation{1}
    param::T
end

# Explicit signature
function (eq::MyODE)(x, y::AbstractArray, α::Number)
    @. y += eq.param * sin(x) * α
end

# Implicit signature
function (eq::MyODE)(res::AbstractArray, x::Number, y::AbstractArray,
                      inc::AbstractArray, α::Number, β::Number=one(α))
    @. res = inc + eq.param * sin(x) * α
end

# Optional: implement analytical solution
function InitialValueProblems.propagate(eq::MyODE, x, y, tau, i)
    # Return analytical solution at x + tau
end
```

## Index

```@index
```

