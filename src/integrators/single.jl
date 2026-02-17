# single-step integrators for explicit first-order ODEs

"""

    ForwardEuler(eq)

Representation of the forward Euler recurrence relation for explicit first-order ODEs.
```math
y ^ {n + 1} = y ^ n + \\tau f ( x ^ n, y ^ n )
```

This is a callable object that applies one step of the forward Euler scheme to update a state vector in-place.

# Fields

- `eq::ODE{1}`: The explicit first-order ODE to solve

"""
struct ForwardEuler{Q<:ODE{1}} <: Integrator{1}
    eq::Q
end

"""

    (scheme::ForwardEuler)(x, y, tau)

Apply one step of the forward Euler scheme.

# Arguments

- `x::Number`: Current time
- `y::AbstractArray`: Current state (modified in-place by accumulating the time-scaled increment)
- `tau::Number`: Time step size

# Returns

The updated time `x + tau`.

"""
function (this::ForwardEuler)(x, y::AbstractArray, tau, _=())
    (; eq) = this
    eq(x, y, tau)
end


"""

    BackwardEuler(eq, buffer)

"""
struct BackwardEuler{Q<:ODE{1}} <: Integrator{1}
    eq::Q
end

function (this::BackwardEuler)(x, y::AbstractArray, tau, (inc,))
    (; eq) = this

    x += tau

    fill!(inc, zero(eltype(inc)))
    sol = nlsolve(inc) do res, inc
        eq(res, x, y, inc, -tau)
    end

    y .+= sol.zero
end


"""

    Midpoint(eq, buffer)

"""
struct Midpoint{Q<:ODE{1}} <: Integrator{1}
    eq::Q
end

function (this::Midpoint)(x, y::AbstractArray, tau, (inc,))
    (; eq) = this

    x += tau / 2

    fill!(inc, zero(eltype(inc)))
    sol = nlsolve(inc) do res, inc
        eq(res, x, y, inc, -tau, one(tau) / 2)
    end

    y .+= sol.zero
end
