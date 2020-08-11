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
    TISE solver for (u₁,u₂,u₃,u₁',u₂',u₃') with IC of (0,0,0,1,1,1)
    Inputs: ϵ~[E], l::Int; stapt/endpt for DE solver~[L], potentials~[L]->[E],
    μ~[M], reltol and maxiters for DE solver
    Outputs: DE solution for (u₁,u₂,u₃,u₁',u₂',u₃')~([L],[L],[L],1,1,1)
"""
function solver(ϵ, # energy ~ [E]
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

# Test solver - can compare to test_rhs_solver in przybytek_sketch.jl
function test_solver()
    plot_ϵ=1e-14u"hartree"
    plot_l=0
    plot_stapt=3.0u"bohr"
    plot_endpt=1e6u"bohr"
    sol=solver(plot_ϵ,plot_l,stapt=plot_stapt,endpt=plot_endpt)
    Rs=LinRange(plot_stapt,plot_endpt,1000)
    sols=sol.(Rs)
    S_us=getindex.(sols,1) # singlet
    T_us=getindex.(sols,2) # triplet
    Q_us=getindex.(sols,3) # quintet
    Splt=plot(ustrip.(Rs), ustrip.(S_us), legend=false, title="Singlet, ϵ=$plot_ϵ")
    Tplt=plot(ustrip.(Rs), ustrip.(T_us), legend=false, title="Triplet, ϵ=$plot_ϵ")
    Qplt=plot(ustrip.(Rs), ustrip.(Q_us), legend=false, title="Quintet, ϵ=$plot_ϵ")
    plot(Splt,Tplt,Qplt,layout=(3,1))
end


"""
Spherical bessel functions
"""
j(l,x)=sphericalbesselj(l,x)
n(l,x)=sphericalbessely(l,x)


"""
    Returns (Aₗ,Bₗ)(R) that match wavefunction to spherical bessel functions
    Inputs: (u,u')(R) solution ~ ([L],1), ϵ~[E], V(R)~[L]→[E], l; μ~[M]
    Output: (Aₗ,Bₗ)(R) for k matching ħ²k(R)²/2m = ϵ-V(R), for that partial wave
"""
function matchAB(sol, # (u,u')(R) ~ [L]->([L],1)
                 ϵ, # energy ~ [E]
                 pot, # potential ~[L]->[E]
                 l; # ang mom
                 μ=0.5*4.002602u"u", # He₂ reduced mass [M]
                 )
    @assert l>=0 "l=$l is not a valid angular momentum"
    ϵₖ(R)=ϵ-pot(R) # kinetic energy
    k=auconvert(sqrt(2*μ*ϵ)/1u"ħ") # kinetic energy wavenumber
    M = R -> [R*j(l,k*R)                    R*n(l,k*R); #Rjₗ(kR), Rnₗ(kR) & derivs
              (l+1)*j(l,k*R)-k*R*j(l+1,k*R) (l+1)*n(l,k*R)-k*R*n(l+1,k*R)]
    AB = R -> ustrip.(M(R))\ustrip.(sol(R))
    return AB
end
