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


################################################################################
# Test functions
################################################################################

# test function for solver - runs for zero potential
# tested successfully 12/08/2020
function test_solver()
    V(bra,ket,R,μ)=0u"hartree" # free
    lhs, rhs = 1.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=a_lookup_generator(0) # s-wave only
    n=length(lookup)
    IC=SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol=solver(lookup, IC, V, ϵ, lhs, rhs)
    Rs=LinRange(lhs,rhs,1000);
    vals = getindex.(sol.(Rs),1)
    plot(austrip.(Rs), austrip.(vals),title="If you see this, solver runs for the trivial case")
end


# test of rotational interaction
function test_H_rot()
    println("Starting test_H_rot")
    V(bra,ket,R,μ)=H_rot(bra,ket,R,μ)
    lhs, rhs = 1.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=a_lookup_generator(0)
    println("Generated lookup vector")
    n=length(lookup)
    IC=SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    println("Generated ICs, about to solve")
    sol=solver(lookup, IC, V, ϵ, lhs, rhs)
    println("Solved for wavefunctions")
    Rs=LinRange(lhs,rhs,1000)
    vals1 = getindex.(sol.(Rs),1)
    vals2 = last.(sol.(Rs))
    vals=hcat(vals1,vals2)
    plot(austrip.(Rs), austrip.(vals),title="If you see this, solver runs for H_rot")
end

# test of electronic interaction
function test_H_el()
    V(bra,ket,R,μ)=H_el(bra,ket,R)
    lhs, rhs = 3.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=a_lookup_generator(0)
    n=length(lookup)
    IC=SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol=solver(lookup, IC, V, ϵ, lhs, rhs)
    Rs=LinRange(lhs,rhs,1000)
    vals1 = getindex.(sol.(Rs),1)
    vals2 = last.(sol.(Rs))
    vals=hcat(vals1,vals2)
    plot(austrip.(Rs), austrip.(vals),title="If you see this, solver runs for H_el")
end

# test of electronic interaction
function test_H_sd()
    V(bra,ket,R,μ)=H_sd(bra,ket,R)
    lhs, rhs = 3.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=a_lookup_generator(1) # s-wave only
    n=length(lookup)
    IC=SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol=solver(lookup, IC, V, ϵ, lhs, rhs)
    Rs=LinRange(lhs,rhs,1000)
    vals1 = getindex.(sol.(Rs),1)
    vals2 = last.(sol.(Rs))
    vals=hcat(vals1,vals2)
    plot(austrip.(Rs), austrip.(vals),title="If you see this, solver runs for H_el")
end


#TODO hyperfine interaction

#TODO in2out_solver. Inputs: defaults=(see notebook). Output: matrix of IC results
