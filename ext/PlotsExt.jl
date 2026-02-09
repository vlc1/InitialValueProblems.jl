module PlotsExt

using Ene4302a
using Plots

Plots.@recipe function f(sol::Solution, x::AbstractVector, i)
    x, sol.(x, i)
end

end
