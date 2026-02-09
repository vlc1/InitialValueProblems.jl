using Plots
using Ene4302a
using Test

@testset "Sinusoidal solution" begin
    x, y = 0., ones(1)

    eq = Sinusoidal()
    z = eq(x, y)

    ic = InitialCondition(x, y)
    sol = Solution(eq, ic)

    y′ = sol(x)

    @test isapprox(y, y′)

    p = plot(sol, 0:1, 1)
    @test isa(p, Plots.Plot)
end
