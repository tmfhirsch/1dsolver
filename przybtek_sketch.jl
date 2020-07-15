#= Playing around with matching coefficients with Przybytek quintet potential,
with the aim of producing a scattering length

Description last updated 14/07/20=#

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using OrdinaryDiffEq, Plots, ComponentArrays, LinearAlgebra
using StepPotential: schrodinger_solver
using Przybytek: przybytek, Unitful, UnitfulAtomic
using Revise


"""TISE solver for IC of (ψ,ψ')=(0,1)
All inputs in units, or assumed to be in atomic units if not.
Inputs: energy, limits; V(R)=przybytek, m=He₂ mass
Output: (ψ,ψ') on RHS of limits"""
function rhs_solver(ϵ;
                    stapt=20.0u"bohr", # location of IC
                    endpt=500.0u"bohr", # RHS to be solved to
                    V=przybytek, # interatomic potential
                    m=2*4.002602u"u" # He₂ mass)
                    )
    function TISE!(du,u,p,x) # TISE gives (ψ,ψ')'=f((ψ,ψ'))
        du.w = u.d # ψ'(x)≡ψ'(x)
        du.d = 2*m*(V(x)-ϵ)*(u.w)/1u"ħ^2" # (ψ'(x))'=2m(V(x)-E)ψ(x)/ħ^2
    end
    w0 = [0.0] #ψ₀ = 0
    d0 = [1.0]u"bohr^-1" #ψ₀'= 1
    IC=ComponentArray(w=w0,d=d0)
    prob=ODEProblem(TISE!,IC,[stapt,endpt])
    sol=solve(prob,Tsit5())
    return sol
end
sol=rhs_solver(10u"hartree",stapt=1.0u"bohr",endpt=150.0u"bohr")
Rs=sol.t
ψs=[i[1] for i in sol.u[:]]
ψplot=plot(ustrip.(Rs), ustrip.(ψs), legend=false)
Vs=przybytek.(Rs)
Vplot=plot(ustrip.(Rs), ustrip.(Vs), legend=false)
plot(ψplot, Vplot, layout=(2,1))
#=
# Test of above
ϵ=1u"hartree"
limits=[10.,11.]u"bohr"
V=przybytek
m=2*4.002602u"u"
function TISE!(du,u,p,x) # TISE gives (ψ,ψ')'=f((ψ,ψ'))
    du.w = u.d # ψ'(x)≡ψ'(x)
    du.d = 2*m*(V(x)-ϵ)*(u.w)/1u"ħ^2" # (ψ'(x))'=2m(V(x)-E)ψ(x)/ħ^2
end
w0 = [0.0] #ψ₀ = 0
d0 = [1.0]u"bohr^-1" #ψ₀'= 1
IC=ComponentArray(w=w0,d=d0)
prob=ODEProblem(TISE!,IC,limits)
sol=solve(prob,Tsit5())
=#

#=
#Simple example demonstrating bug I find with second order ODEs
function deriv!(du,u,p,t)
    du.w = u.d
    du.d = (u.w)u"s^-2"
end
w0 = [0.0]
d0 = [1.0]u"s^-1"
Δt = 1u"s"
u0 = ComponentArray(w=w0, d=d0)
prob = ODEProblem(deriv!, u0, (0.0u"s", Δt))
sol = solve(prob,Vern8(), dt=1e-1u"s", adaptive=false)
=#
