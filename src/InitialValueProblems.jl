module InitialValueProblems

using NLsolve

export OrdinaryDifferentialEquation,
       Propagator,
       propagate,
       Integrator,
       ForwardEuler,
       BackwardEuler,
       Midpoint,
       Sinusoidal
#    Oscillator

"""

    OrdinaryDifferentialEquation{N} <: Function

Abstract base type for ordinary differential equations of order `N`.

The ODE is assumed to be explicitly defined in the form:

```math
y ^ {(n)} (t) = f(t, y(t), \\dot{y}(t), \\ldots, y ^ {(n-1)}(t))
```

Subtypes should implement the call signatures required by the integrators.

"""
abstract type OrdinaryDifferentialEquation{N} <: Function end

const ODE = OrdinaryDifferentialEquation

"""

    Propagator{Q<:ODE}

A callable object used to evaluate the analytic solution of an ODE.

"""
struct Propagator{Q<:ODE}
    eq::Q
end

"""

    (this::Propagator)(x, y::AbstractArray, tau, i)

Apply the propagator to compute the `i`th component of the solution at time `x + tau` based on the current state `y` at time `x`.

# Returns

The value of the `i`th component of the solution at time `x + tau`.

"""
function (this::Propagator)(x, y::AbstractArray, tau, i)
    (; eq) = this

    propagate(eq, x, y, tau, i)
end

"""

    propagate(eq::ODE, x, y, tau, i)

Evaluate the analytic solution for the `i`th component when the model provides it.

Specialized methods for specific ODE types should be implemented to compute the solution based on the parameters of the ODE and the initial conditions.

"""
function propagate end

include("odes/sinusoidal.jl")
#include("odes/oscillators.jl")

"""

    Integrator{N}

Abstract base type for numerical time-stepping schemes for solving first-order initial value problems.

# Type Parameter

`N::Int`: The number of previous steps required by the scheme.

- `N=1` for single-step methods (e.g., Forward Euler)
- `N=2` for two-step methods (e.g., Adams-Bashforth), etc.

"""
abstract type Integrator{N} end

include("integrators/single.jl")

end
