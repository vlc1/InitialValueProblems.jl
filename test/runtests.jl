using Plots
using InitialValueProblems
using Test

@testset "Sinusoidal solution" begin
    x, y, eq = 0., ones(1), Sinusoidal()
    cache = (similar(y),)

    ref = Propagator(eq)
    fwd = ForwardEuler(eq)
    bwd = BackwardEuler(eq)
    mid = Midpoint(eq)

    xs = 0.:0.1:1.

    p = plot(ref, xs, ones(1), 1)
    scatter!(p, fwd, xs, ones(1), 1)
    scatter!(p, bwd, xs, ones(1), 1, cache)
    scatter!(p, mid, xs, ones(1), 1, cache)

    @test isa(p, Plots.Plot)
end
