#= Playing around with matching coefficients with Przybytek quintet potential,
with the aim of producing a scattering length

Description last updated 14/07/20=#

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using Revise
using OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays
using Unitful, UnitfulAtomic
using Przybytek: przybytek


"""TISE solver for IC of (ψ,ψ')=(0,1)
All inputs in units, or assumed to be in atomic units if not.
Inputs: energy, limits; V(R)=przybytek, m=He₂ mass
Output: (ψ,ψ') on RHS as fn of x (a₀)"""
function rhs_solver(ϵ; # energy [E]
                    stapt=20.0u"bohr", # location of IC [L]
                    endpt=500.0u"bohr", # RHS to be solved to [L]
                    U=przybytek, # interatomic potential [L->E]
                    m=0.5*4.002602u"u" # He₂ reduced mass [M]
                    )
    # check units are dimensionally correct
    @assert dimension(2*m*(U(stapt)-ϵ)/1u"ħ^2")==dimension(1u"m^-2") "TISE units mismatch"
    # convert and strip units before TISE
    ϵ = austrip(ϵ)
    stapt, endpt = austrip(stapt), austrip(endpt)
    V(x)=austrip(U((x)u"bohr"))
    m = austrip(m)
    ħ = austrip(1u"ħ")
    # TISE solver
    function TISE(u,p,x) # TISE gives (ψ,ψ')'=f((ψ,ψ'))
        du = u[2] # ψ'(x)≡ψ'(x)
        dd = 2*m*(V(x)-ϵ)*u[1]/ħ^2 # (ψ'(x))'=2m(V(x)-E)ψ(x)/ħ^2 {𝐋⁻²}
        SVector{2}([du,dd]) # (ψ,ψ')'
    end
    IC = SVector{2}([0.0, 1.0]) # (ψ,ψ')=(0,1)
    prob=ODEProblem(TISE,IC,(stapt,endpt))
    sol=solve(prob,Tsit5(),reltol=1e-10)
    Usol(x)=sol(ustrip(uconvert(u"bohr",x)))
    return Usol
end

plot_ϵ=1e-10u"hartree"
plot_stapt=1.0u"bohr"
plot_endpt=20000.0u"bohr"
sol=rhs_solver(plot_ϵ,stapt=plot_stapt,endpt=plot_endpt)
Rs=LinRange(plot_stapt,plot_endpt,1000)
ψs=[i[1] for i in sol.(Rs)]
plot(ustrip.(Rs), ustrip.(ψs), legend=false, title="$plot_ϵ")
#savefig("prz_wfn_ϵ-$plot_ϵ.png")

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
