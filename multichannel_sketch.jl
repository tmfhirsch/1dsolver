#=
Initial goal: solve TISE for singlet, triplet and quintet channels,  with
no coupling of channels. TISE rewritten in terms of matrix vector mult.
Description last updated 27/07/20
=#
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using Revise
using OrdinaryDiffEq, LinearAlgebra, StaticArrays, SpecialFunctions
using Unitful, UnitfulAtomic
using Potentials: Singlet, Triplet, Quintet
using Plots

"""
TISE solver for (u₁,u₂,u₃,<derivs>) with IC of (0,0,0,1,1,1)
"""

function solver(k, # wavenumber [L]⁻¹
                l::Int; # angular momentum
                stapt=3.0u"bohr", # location of IC [L]
                endpt=100.0u"bohr", # RHS to be solved to [L]
                pot1=Singlet, # interatomic potentials [L]->[E]
                pot2=Triplet,
                pot3=Quintet,
                μ=0.5*4.002602u"u", # He₂ reduced mass [M]
                reltol=1e-10, #relative tolerance for DE solver
                maxiters=1e6 #max iterations for DE solver
                )
    ϵ=auconvert(k^2*1u"ħ^2"/(2*μ)) # E=ħ²k²/2μ
    # add centrifugal potential
    V1 = R -> auconvert(pot1(R)+l*(l+1)u"ħ^2"/(2*μ*R^2))
    V2 = R -> auconvert(pot2(R)+l*(l+1)u"ħ^2"/(2*μ*R^2))
    V3 = R -> auconvert(pot3(R)+l*(l+1)u"ħ^2"/(2*μ*R^2))
    # convert and strip units before TISE
    V1⁰ = x -> austrip(V1((x)u"bohr")) # unitless -> unitless potential
    V2⁰ = x -> austrip(V2((x)u"bohr"))
    V3⁰ = x -> austrip(V3((x)u"bohr"))
    ϵ⁰ = austrip(ϵ) # strip energy
    stapt⁰, endpt⁰ = austrip(stapt), austrip(endpt) # strip start/end points
    μ⁰ = austrip(μ) # strip mass
    ħ⁰ = austrip(1u"ħ") # strip ħ
    # TISE solver
    function TISE(u,p,x) # TISE gives (u₁,...,u₃')'=f((u₁,...,u₃'))
        t1 = 2*μ⁰*(V1⁰(x)-ϵ⁰)/ħ⁰^2 # (u'(x))'=2m(V(x)-E)u(x)/ħ^2 {𝐋⁻²}
        t2 = 2*μ⁰*(V2⁰(x)-ϵ⁰)/ħ⁰^2
        t3 = 2*μ⁰*(V3⁰(x)-ϵ⁰)/ħ⁰^2
        M = @SMatrix   [0 0 0 1 0 0;
                        0 0 0 0 1 0;
                        0 0 0 0 0 1;
                       t1 0 0 0 0 0;
                       0 t2 0 0 0 0;
                       0 0 t3 0 0 0]
        M*u
    end
    # Specify problem and solve
    IC = SVector{6}([0.0,0.0,0.0, 1.0,1.0,1.0]) # (∀i uᵢ=0 ⩓ uᵢ'=1)
    prob=ODEProblem(TISE,IC,(stapt⁰,endpt⁰))
    sol_unitless=solve(prob,Tsit5(),reltol=reltol,maxiters=maxiters)
    # Add back units
    sol_input = x -> sol_unitless(austrip(x)) # add unit input
    sol = x -> sol_input(x).*SA[1u"bohr",1u"bohr",1u"bohr",1,1,1] # u~[L],u'~1
    return sol
end

# Test solver
function test_solver()
    k=1e-5u"bohr^-1"
    l=0
    plot_k=1e-5u"bohr^-1"
    plot_l=0
    plot_stapt=3.0u"bohr"
    plot_endpt=1e6u"bohr"
    sol=solver(plot_k,plot_l,stapt=plot_stapt,endpt=plot_endpt)
    Rs=LinRange(plot_stapt,plot_endpt,1000)
    sols=sol.(Rs)
    S_us=getindex.(sols,1) # singlet
    T_us=getindex.(sols,2) # triplet
    Q_us=getindex.(sols,3) # quintet
    Splt=plot(ustrip.(Rs), ustrip.(S_us), legend=false, title="Singlet, k=$plot_k")
    Tplt=plot(ustrip.(Rs), ustrip.(T_us), legend=false, title="Triplet, k=$plot_k")
    Qplt=plot(ustrip.(Rs), ustrip.(Q_us), legend=false, title="Quintet, k=$plot_k")
    plot(Splt,Tplt,Qplt,layout=(3,1))
end
