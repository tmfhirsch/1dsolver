#= Multichannel equations solver. Uses symmetrised |a⟩ basis states and all
of the interactions except for the Zeeman interaction.
All states up to and including the lmax parameter are used
Description last updated 12/08/20=#
using Revise
using HalfIntegers, LinearAlgebra, StaticArrays, OrdinaryDiffEq, WignerSymbols
using Unitful, UnitfulAtomic

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")
using Interactions
using StateStructures
using Plots


""" Multichannel TISE solver
    Input: lookup vector, IC~SA{{[L],...,1,...}}, pot~|a⟩×|a⟩×[L]×[M]→[E],
    energy~[E], lhs~[L], rhs~[L]
    Output: sol(R) [where R ∈ (lhs,rhs)] ~ IC"""
function solver(lookup, IC, ϵ, lhs, rhs; μ=0.5*4.002602u"u")
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
    # precalculate coupling coefficients for spin-dipole interaction
    C_sd = zeros(n,n) # initialise
    for i=1:n, j=1:n
        C_sd[i,j] = H_sd_coeffs(lookup[i], lookup[j])
    end
    # TISE differential equation
    function TISE(u,p,x)
        # Construct V(R) matrix
        V = zeros(n,n)u"hartree" # initialise
        for i=1:n, j=1:n
            V[i,j] = H_rot(lookup[i],lookup[j], x*1u"bohr", μ) # rotational
            V[i,j]+= H_el(lookup[i],lookup[j], x*1u"bohr") # electronic
            V[i,j]+= C_sd[i,j]*H_sd_radial(x*1u"bohr") # spin-dipole
            #TODO hyperfine interaction
            #TODO zeeman interaction (will also need to fix channels for this)
        end
        V⁰=austrip.(V) # strip units from V
        M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matix
        D = SMatrix{2*n,2*n}([0*I I
                              M 0*I])
        D*u # ⃗u' = D . ⃗u
    end
    # strip units from IC
    IC⁰ = austrip.(IC)
    # solve
    prob=ODEProblem(TISE,IC⁰,(lhs⁰,rhs⁰))
    sol_unitless=solve(prob,BS5()) #BS5 used to avoid Tsit5 querying x<3a₀<lhs⁰
    # add units back
    units = SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol = x -> sol_unitless(austrip(x)).*units
    return sol
end

#=
""" Old multichannel TISE solver, which took an arbitrary potential.
    Input: lookup vector, IC~SA{{[L],...,1,...}}, pot~|a⟩×|a⟩×[L]×[M]→[E],
    energy~[E], lhs~[L], rhs~[L]
    Output: sol(R) [where R ∈ (lhs,rhs)] ~ IC"""
function solver(lookup, IC, pot, ϵ, lhs, rhs; μ=0.5*4.002602u"u")
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
            @assert dimension(IC[i,j])==dimension(1u"m") "IC[$i,$j] not dimensionless"
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
        V = zeros(n,n)u"hartree" # initialise
        for i=1:n, j=1:n
            V[i,j] = pot(lookup[i], lookup[j], x*1u"bohr", μ)
        end
        V⁰=austrip.(V) # strip units from V
        M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matix
        D = SMatrix{2*n,2*n}([0*I I
                              M 0*I])
        D*u # ⃗u' = D . ⃗u
    end
    # strip units from IC
    IC⁰ = austrip.(IC)
    # solve
    prob=ODEProblem(TISE,IC⁰,(lhs⁰,rhs⁰))
    sol_unitless=solve(prob,Tsit5())
    # add units back
    units = SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol = x -> sol_unitless(austrip(x)).*units
    return sol
end
=#

################################################################################
# Test functions
################################################################################

# test function for solver - runs for zero potential
# tested successfully 12/08/2020
function test_solver(lmax=0)
    lhs, rhs = 3.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=a_lookup_generator(lmax)
    n=length(lookup)
    # construct ICs
    IC=SMatrix{2*n,n}([fill(0.0u"bohr",n,n)
                       I])
    println("Starting to solve for wavefunctions, lmax=$lmax")
    @time begin sol=solver(lookup, IC, ϵ, lhs, rhs)
    end
    println("Plotting...")
    Rs=LinRange(lhs,rhs,1000);
    vals = getindex.(sol.(Rs),1,1)
    plot(austrip.(Rs), austrip.(vals),title="It works! Plotting wavefn of first channel", legend=false)
end

#TODO in2out_solver. Inputs: defaults=(see notebook). Output: matrix of IC results
