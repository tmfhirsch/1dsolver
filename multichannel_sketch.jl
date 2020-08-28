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


""" Multichannel TISE solver. Produces info needed for K matrix.
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
    sol_unitless=solve(prob,Tsit5())
    # add units back
    units = SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol = x -> sol_unitless(austrip(x)).*units
    # construct 𝐤 vector
    𝐤=zeros(n,1)*1.0u"bohr^-1"

    return sol
end

using SpecialFunctions, LinearAlgebra
""" Solver for 𝐊 matrix (following Mies eqn (3.8)), where 𝐅=𝐉-𝐍𝐊.
    Inputs:
        R~[L] the asymptotic radial distance to match K at;
        sol() solution for [𝐆; 𝐆'] where 𝐆 is matrix of wavefunctions for
    different initial conditions;
        𝐤 ~ [L]⁻¹ vector of wavenumber for each channel at rhs;
    **Note: this function assumes it is only being given open channels**
        𝐥 vector of l quantum numbers for each channel.
    **Note: sol, 𝐤, 𝐥 must share the same ordering/number of channels**
    Output:
        𝐊 ~ n×n matrix (n=number of channels considered)"""
function K_matrix(R, sol, 𝐤, 𝐥)
    # match for A, B where G=J.A-N.B
    #construct G, G' matrices to match for A, B with
    eval=sol(R)
    n=Int(size(eval,1)/2) # n = number of channels. Assumes sol in above form.
    @assert size(eval,2) == n "solution matrix not of shape 2n × n"#BUG
    G, G⁻ = austrip.(copy(eval[1:n,1:n])), copy(eval[n+1:2*n,1:n])
    # solve for A,B
    A, B = zeros(n,n), zeros(n,n) # initialise
    for i in 1:n, j in 1:n
        # construct [jᵢ nᵢ; jᵢ' nᵢ'] matrix
        k=𝐤[i]
        l=𝐥[i]
        j=austrip(sqrt(k)*R)*sphericalbesselj(l,k*R)
        j⁻=austrip(sqrt(k))*((l+1)*sphericalbesselj(l,k*R)
            -k*R*sphericalbesselj(l+1,k*R))
        n=austrip(sqrt(k)*R)*sphericalbessely(l,k*R)
        n⁻=austrip(sqrt(k))*((l+1)*sphericalbessely(l,k*R)
            -k*R*sphericalbessely(l+1,k*R))
        Gᵢⱼ, G⁻ᵢⱼ = G[i,j], G⁻ᵢⱼ = G⁻ᵢⱼ
        AB = [j -n; j⁻ -n⁻]\[Gᵢⱼ; G⁻ᵢⱼ] # AB≡[Aᵢⱼ; Bᵢⱼ]
        A[i,j], B[i,j] = AB
    end
    𝐊 = B*inv(A)
end

"""SKETCH: "the whole package".
    Produces lookup, ICs, solves and outputs K matrix"""
function K_solver(lmax, ϵ, lhs, rhs; μ=0.5*4.002602u"u")
    # geberate states
    lookup=a_lookup_generator(lmax)
    n=length(lookup)
    # construct ICs
    IC=SMatrix{2*n,n}([fill(0.0u"bohr",n,n)
                       I])
    # solver for wavefunctions
    #TODO
    # generate 𝐤 vector for K matrix solver
    #TODO
    # generate 𝐥 vector for K matrix solver
    #TODO
end

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

function test_K_matrix(lmax=0)
    lhs, rhs = 3.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=a_lookup_generator(lmax)
    n=length(lookup)
    # construct ICs
    IC=SMatrix{2*n,n}([fill(0.0u"bohr",n,n)
                       I])
