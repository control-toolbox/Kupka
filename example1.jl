#
#using Pkg
#Pkg.activate(".")
#
using OptimalControl
using MINPACK # NLE solver
using NLsolve

# cellule 2
t0 = 0
tf = 1
x0 = [-1 ; 0]
xf = [0 ; 0]
@def ocp begin
    t ∈ [ t0, tf ], time
    x ∈ R^2, state
    u ∈ R, control
    x(t0) == x0
    x(tf) == xf
    #ẋ₁(t) == x₂(t)
    #ẋ₂(t) == u(t)
    ẋ(t) ==  [x₂(t) ; u(t)]
    ∫( 0.5u(t)^2 ) → min
end;

# Compute the  Flow
u(x, p) = p[2]

f = Flow(ocp, u)
p0 = [0;0]
println(f((t0, tf), x0, p0, saveat=[]).ode_sol(tf))

# Define the shooting function
S(p0) = f((t0, tf), x0, p0, saveat=[]).ode_sol(tf)[1:2] - xf;

println(S([0. ; 0.]))
println(S([12. ; 6.]))

println(S([1. ; 1.]))

# solve the shooting equation
nl_sol = nlsolve(S, [0.,0.],autodiff = :forward)

# ajouter une fonction qui renvoie une solution du pr de control de type OptimalControl_solution afin de faire le plot

println("nl_sol = ", nl_sol)
p_sol = nl_sol.zero
println("p_sol = ", p_sol)

# draw the solution
#plot(f((t0, tf), x0, p_sol, saveat=[]).ode_sol)

# compute the solution: state, costate, control...
# if we don't have the saveat  option, then we obtain too few points from the numerical solution to have  a nice plot
Δt = (tf - t0)/100

flow_sol = f((t0, tf), x0, p_sol, saveat = Δt)    

plt = plot(flow_sol)

plot(plt)

flow_p0 = f((t0, tf), x0, p0, saveat = Δt)          

plt = plot!(plt,flow_p0) # pb car pas les même ordonnées


# Conjugate points



      







