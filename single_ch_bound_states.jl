# calculate bound state energies, for use with identifying FB res bound states

using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using Unitful, UnitfulAtomic
using OrdinaryDiffEq, LinearAlgebra
using Interactions, StateStructures
using Plots, Plots.PlotMeasures

const num=1000 # number of grid points to save

"""Callback function for renormalisation of wavefunction. Code by DC"""
function CreateRenormalisedCallback(maxval=1e5)
    maxvalsqr = maxval^2
    condition = (u,t,int) -> any(abs2.(u) .> maxvalsqr)
    DiscreteCallback(condition, _Renormalise!, save_positions=(true,false))
end
function _Renormalise!(int)
    maxval = sqrt.(maximum(abs2.(int.u), dims=1))
    int.u ./= maxval
    for i = eachindex(int.sol.u)
        int.sol.u[i] ./= maxval
    end
    nothing
end

""" Single channel solver, *without* H_sd, H_zee, Γ_gms"""
function sng_solver(lookup, IC, ϵ, lhs, rhs; μ=0.5*4.002602u"u")
    # Check units of ϵ, μ, lhs, rhs
    @assert dimension(ϵ)==dimension(1u"J") "ϵ not an energy"
    @assert dimension(μ)==dimension(1u"g") "μ not a mass"
    @assert dimension(lhs)==dimension(1u"m") "lhs not a length"
    @assert dimension(rhs)==dimension(1u"m") "rhs not a length"
    # Check length and units of IC
    n = length(lookup) # number of channels
    @assert size(IC)[1]==2*n "Initial condition has wrong number of channels"
    if length(size(IC))==1 # IC a vector
        for i=1:n # check wavefunction entries
            @assert dimension(IC[i])==dimension(1u"m") "IC[$i] not a length"
        end
        for i=(n+1):2*n # check derivative entries
            @assert dimension(IC[i])==dimension(1) "IC[$i] not dimensionless"
        end
    elseif length(size(IC))==2 # IC a matrix of initial condition vectors
        for i=1:n, j=1:size(IC)[2] # check wavefunction entries
            @assert dimension(IC[i,j])==dimension(1u"m") "IC[$i,$j] not a length"
        end
        for i=(n+1):2*n, j=1:size(IC)[2] # check derivative entries
            @assert dimension(IC[i,j])==dimension(1) "IC[$i,$j] not dimensionless"
        end
    else # IC not a 1D vector or 2D array
        error("IC not 1D or 2D")
    end
    # strip units from constants
    ϵ⁰, μ⁰ = austrip(ϵ), austrip(μ)
    ħ⁰ = austrip(1.0u"ħ")
    lhs⁰, rhs⁰ = austrip(lhs), austrip(rhs)
    # TISE differential equation
    function TISE(u,p,x)
        # Construct V(R) matrix
        V = zeros(Float64,n,n)u"hartree" # initialise
        for i=1:n, j=1:n
            V[i,j] = H_rot(lookup[i],lookup[j], x*1u"bohr", μ) # rotational
            V[i,j]+= H_el(lookup[i],lookup[j], x*1u"bohr") # electronic
        end
        V⁰=austrip.(V) # strip units from V
        M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matix
        D = ([0*I I
              M 0*I])
        D*u # ⃗u' = D . ⃗u
    end
    # strip units from IC
    IC⁰ = austrip.(IC)
    # solve
    prob=ODEProblem(TISE,IC⁰,(lhs⁰,rhs⁰))
    callback=CreateRenormalisedCallback()
    sol_unitless=solve(prob,Tsit5(),reltol=1e-10,dense=true, callback=callback)
    # add units back
    units = vcat(fill(1.0u"bohr",n),fill(1.0,n))
    sol = x -> sol_unitless(austrip(x)).*units
    return sol
end

BOmins = let
    grid=LinRange(5u"bohr",7u"bohr",10_000)
    S0ket = filter(x->x.S==0,SmS_lookup_generator(2))[1]
    V0(R) = H_el(S0ket,S0ket,R)
    S1ket = filter(x->x.S==1,SmS_lookup_generator(2))[1]
    V1(R) = H_el(S1ket,S1ket,R)
    S2ket = filter(x->x.S==2,SmS_lookup_generator(2))[1]
    V2(R) = H_el(S2ket,S2ket,R)
    minimum(V0.(grid)),minimum(V1.(grid)),minimum(V2.(grid))
end

""" Find lhs and rhs where H_el + H_rot ~= given energy
    Input: |Smₛ⟩ ket, ϵ~[E]
    Output: lhs~[L], rhs~[L]"""
function lhsrhs(ket::SmS_ket, ϵ; μ=0.5*4.002602u"u", δR=1e-2u"bohr", rhsbound=1e5u"bohr")
    V(R)=H_el(ket,ket,R)+H_rot(ket,ket,R,μ)
    Δ(R)=ϵ-V(R)
    lhs, rhs = 0e0u"bohr", 0e0u"bohr" # initialise
    needle = 3e0u"bohr" # initialise first search
    @assert sign(Δ(needle))==-1 "Not starting lhs search outside well"
    while sign(Δ(needle))==-1
        needle += δR
    end
    lhs=needle # save result of first search
    needle=7e0u"bohr" # initialise second search
    @assert sign(Δ(needle))==1 "Not starting rhs search inside well"
    while sign(Δ(needle))==1
        needle += δR
        needle < rhsbound || error("rhs search reached $rhsbound")
    end
    rhs=needle # save result of second search
    return lhs, rhs
end

"""Number of zero crossings calculator
    Input: S∈{0,1,2}, l::Int, Energy~[E]
    Output: number of zero crossings (on grid as defined above)"""
function numcrossings(S::Int, l::Int, ϵ::Unitful.Energy)
    #checks
    @assert ϵ<=0u"hartree" "ϵ > 0"
    @assert S in [0,1,2] "S ∉ {0,1,2}"
    S==0 && ϵ<BOmins[1] && error("ϵ below well depth")
    S==1 && ϵ<BOmins[2] && error("ϵ below well depth")
    S==2 && ϵ<BOmins[3] && error("ϵ below well depth")
    # finished checks
    lookup = [filter(x->(x.S==S && x.l == l), SmS_lookup_generator(l))[1]]
    IC = [0e0u"bohr"; 1e0]
    lhs, rhs = lhsrhs(lookup[1],ϵ)
    sol=sng_solver(lookup, IC, ϵ, lhs, rhs)
    grid=LinRange(lhs,rhs,num)
    ψvals=getindex.(sol.(grid), 1)
    crossings=count(ψvals[1:end-1].*ψvals[2:end].<0u"bohr^2")
    #plt=plot(ustrip.(grid),ustrip.(ψvals)) plot of wavefunction
    return crossings
end

####################Script work###################################3
function sketch()
    l=2; S=2
    ϵ=-7.044e-5u"hartree"; δϵ=1e-10u"hartree"
    while numcrossings(S,l,ϵ)==12
        println(ϵ)
        ϵ-=δϵ
    end
    println(ϵ)
end

# S=0, l=2 has 27 crossings at 0 energy; S=2, l=2 has 14
# manually entering values into here # 1e-10 Eₕ precision
S0l2_bstate_en=[-6.7031e-6u"hartree",-4.20576e-5u"hartree",-1.273626e-4u"hartree"]
S2l2_bstate_en=[-1.0794e-6u"hartree",-1.87529e-5u"hartree",-7.04485e-5u"hartree"]
