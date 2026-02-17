module PlotsExt

using InitialValueProblems
using Plots

using InitialValueProblems: Integrator

Plots.@recipe function f(flow::Propagator, xs::AbstractRange, y, i=firstindex(y))
    n, tau = length(xs), step(xs)

    yis = similar(y, n)

    x = first(xs)
    j = firstindex(yis)
    yis[j] = y[i]

    delta = zero(tau)

    while j < lastindex(yis)
        j = nextind(yis, j)
        delta += tau
        yis[j] = flow(x, y, delta, i)
    end

    xs, yis
end

Plots.@recipe function f(scheme::Integrator{1}, xs::AbstractRange, y, i=firstindex(y), cache=())
    n, tau = length(xs), step(xs)

    yis = similar(y, axes(xs))

    j = firstindex(xs)
    x, yis[j] = xs[j], y[i]

    while j < lastindex(xs)
        scheme(x, y, tau, cache)

        j = nextind(xs, j)
        x, yis[j] = xs[j], y[i]
    end

    xs, yis
end

end
