using LinearAlgebra
using ForwardDiff
using DifferentialEquations
using Plots

#using MINPACK # NLE solver
#using NLsolve
using LaTeXStrings
using Revise

# analytic solutions
# ------------------
# state
"""
    Compute the state equation 
    x_t = x(t,x0,p0)
"""
function x(t,x0,p0)
    λ = atanh(p0/x0)
    β = x0/cosh(λ)
    α = -λ*β
    return β*cosh((t-α)/β)
end
"""
    Compute the costate equation 
    p_t = p(t,x0,p0)
"""
function p(t,x0,p0)
    λ = atanh(p0/x0)
    β = x0/cosh(λ)
    α = -λ*β
    return β*sinh((t-α)/β)
end

"""
    Compute the second member of the Hamiltanian Flow
    Hvec = H_vec(z)
"""
function H_vec(z)
    x = z[1]
    p = z[2]
    den = sqrt(x^2-p^2)
    return [sign(x)*p/den ; x/den]
end

"""
    Compute the second member of the Jacobi equation
    rhs  =  jacobi_ana(δz, par, t)
"""
function jacobi_ana(δz, par, t)
    x0, p0 = par   # (x_0,p_0)
    x_t = x(t,x0,p0)
    p_t = p(t,x0,p0)
    temp = 1/(x_t^2 - p_t^2)^(3/2)
    return temp*[-x_t*p_t x_t^2 ; -p_t^2 x_t*p_t]*δz
end

"""
    Compute the flow of the Jacobian equation for the initial condition δz(0)=(0,1)
    sol = flow_jacobi_ana(t0tf,x0,p0)
    sol is return by the solve function of  the DifferentialEquations package
"""
function flow_jacobi_ana(t0tf,x0,p0)
    t0, tf = t0tf
    prob = ODEProblem(jacobi_ana,[0;1],t0tf,[x0;p0])
    sol = DifferentialEquations.solve(prob,reltol = 1e-10, abstol = 1e-10)
    return sol
end

function F_ana(τ,p0)
    tspan = (0,τ)
#    return [flow_jacobi(tspan,x₀,p0)(τ)[1]]
    return flow_jacobi_ana(tspan,x0,p0)(τ)[1] 
end
function rhs_path_ana(tau, par, p0)
    τ = tau[1]
    Ftau(p0) = F_ana(τ, p0)
    δz = flow_jacobi_ana((0,τ),x0,p0)(τ)
    derivee_τ = jacobi_ana(δz, [x0,p0], τ)[1]
    derivee_p0 = ForwardDiff.derivative(Ftau, p0)
    return [-(1/derivee_τ)*derivee_p0]
end

