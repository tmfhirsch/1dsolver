#= Goal: verify my method of solving the multichannel TISE in the case of a
square barrier with coupling, identical to the situation in PHYS4304
Description last updated 30/07/2020=#

using Revise
using OrdinaryDiffEq, LinearAlgebra, StaticArrays, SpecialFunctions
using Unitful, UnitfulAtomic
using Plots

"""
Coupled square barrier
Inputs: R~[L]; lhs=1a₀,rhs=2a₀, V₀=10Eₕ, Vₓ=1Eₕ
Outputs: Potential as a matrix in the basis of channels, units of energy
"""
function coupled_square_barrier(R;lhs=1u"bohr",rhs=2u"bohr",V₀=10u"hartree",Vₓ=1u"hartree")
    if R<lhs # left of barrier
        M = SA[0 0; 0 0]
    elseif R<rhs # inside barrier
        M = SA[V₀ Vₓ; Vₓ V₀]
    else # right of barrier
        M = SA[0 0; 0 0]
    end
    return M
end

"""
    TISE solver for potential
    Input: energy~[E], other constants
    Output: matrix relating (u₁,u₂,u₁',u₂') on RHS to ICs
"""
function solver(ϵ=5u"hartree", V=coupled_square_barrier, μ=0.5*4.002602u"u", stapt=-10.0u"bohr",endpt=12.0u"bohr")
    @assert dimension(2*μ*(1u"hartree"-ϵ)/1u"ħ^2")==dimension(1u"m^-2") "TISE units mismatch"
    ϵ⁰=austrip(ϵ)
    V⁰ = R->austrip.(V((R)u"bohr"))
    μ⁰=austrip(μ)
    ħ⁰=austrip(1u"ħ")
    stapt⁰, endpt⁰ = austrip(stapt), austrip(endpt)
    # TISE solver
    function TISE(u,p,x) # TISE gives (u₁,...,u₃')'=f((u₁,...,u₃'))
        E_matrix=SA[ϵ⁰ 0; 0 ϵ⁰] # energy term in schrodinger equation #NOTE can have diff E_channel
        V_matrix=V⁰(x) # potential term
        DD_matrix=-2*μ⁰*(E_matrix-V_matrix)/ħ⁰^2 # second derivatives
        M = SMatrix{4,4}([zeros(2,2) I;
                          DD_matrix zeros(2,2)])
        M*u
    end
    # Iterating over basis vectors for ICs to get map from ICs to RHS values.
    # TODO: nicer way than iterating over each of the 2x2=4 elements
    IC1, IC2, IC3, IC4 = SA[1.,0,0,0], SA[0,1.,0,0], SA[0,0,1.,0], SA[0,0,0,1.]
    prob1=ODEProblem(TISE,IC1,(stapt⁰,endpt⁰))
    sol1_unitless=solve(prob1,Tsit5())
    sol1 = x -> sol1_unitless(austrip(x)).*SA[1u"bohr",1u"bohr",1,1]
    prob2=ODEProblem(TISE,IC2,(stapt⁰,endpt⁰))
    sol2_unitless=solve(prob2,Tsit5())
    sol2 = x -> sol2_unitless(austrip(x)).*SA[1u"bohr",1u"bohr",1,1]
    prob3=ODEProblem(TISE,IC3,(stapt⁰,endpt⁰))
    sol3_unitless=solve(prob3,Tsit5())
    sol3 = x -> sol3_unitless(austrip(x)).*SA[1u"bohr",1u"bohr",1,1]
    prob4=ODEProblem(TISE,IC4,(stapt⁰,endpt⁰))
    sol4_unitless=solve(prob4,Tsit5())
    sol4 = x -> sol4_unitless(austrip(x)).*SA[1u"bohr",1u"bohr",1,1]
    return [sol1(endpt) sol2(endpt) sol3(endpt) sol4(endpt)]
end
