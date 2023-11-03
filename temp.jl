#
#using Pkg
#Pkg.activate(".")
#
using Formatting
using ForwardDiff
using MINPACK
using OptimalControl
using Plots
using Plots.PlotMeasures # for leftmargin, bottommargin


println("OCP")

t0 = 0
tf = 1
x0 = -1
xf = -0.451475
α  = 0.0
@def ocp begin
    t ∈ [ t0, tf ], time
    x ∈ R, state
    u ∈ R, control
    x(t0) == x0
    x(tf) == xf
    ẋ(t) == -α*x(t) + u(t) * x(t)
    ∫( 0.5u(t)^2 + 0.5*(1-x(t))^2) → min
end


println("Shooting function")

u(x, p) = p;

φ = Flow(ocp, u); # flow from the optimal control problem ocp and maximising control u

# definition of the state projection
π(x, p) = x
π(z::Tuple{Number, Number}) = π(z...) # z = (x, p)

# shooting function
S(p0) = π( φ(t0, x0, p0, tf) ) - xf;

Jₛ(p0) = ForwardDiff.jacobian(x -> [S(x[1])], [p0])[1, 1]; # S'(p0)

println("Shooting method")

p0 = 1.5                                # initial costate for the Newton solver

#    iterates = [p0]                         # list of iterates
#    iterations = 5                          # number of iterations
#    for i ∈ 1:iterations
#        d  = - Jₛ(p0) \ S(p0)           # Newton direction
#        p0 = p0 + d                     # update
#        push!(iterates, p0)             # save
#    end

println("Shooting method: MINPACK")

# Resolution with Minpack hybrj solver
nle = (s, ξ) -> s[1] = S(ξ[1])          # auxiliary function
ξ = [ p0 ]                              # initial guess
indirect_sol = fsolve(nle, ξ)           # resolution of S(p0) = 0
p0_solution  = indirect_sol.x[1]        # costate solution

println(indirect_sol)
p0 = p0_solution
println("p0 solution = ", p0_solution)
println("JS(p0 sol) = ", Jₛ(p0))

#
println("Plot")

#
include("temp_plot.jl")

#
idx = 3

# copy the plt_flow
plt_flow_copy = deepcopy(plt_flow)
plt_shoot_copy = deepcopy(plt_shoot)

# add extremals
#    if idx ≥ 0
#        plot_extremals!(plt_flow_copy, iterates, idx)
#    end

# add iterates on the plot of shooting function
#    if idx ≥ 0
#        plot_iterates!(plt_shoot_copy, iterates, idx)
#    end


# plot all
plot(plt_flow_copy, plt_shoot_copy, layout=(1,2), size=(900, 450),
        leftmargin=5mm, bottommargin=5mm)

