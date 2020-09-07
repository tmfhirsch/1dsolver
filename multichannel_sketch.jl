#= Multichannel equations solver. Uses symmetrised |a⟩ basis states and all
of the interactions except for the Zeeman interaction.
All states up to and including the lmax parameter are used
Description last updated 12/08/20=#
using Revise
using HalfIntegers, LinearAlgebra, StaticArrays, OrdinaryDiffEq, WignerSymbols
using Unitful, UnitfulAtomic
using Plots
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using Interactions, StateStructures

""" Multichannel TISE solver. Produces info needed for K matrix.
    Input: lookup vector, IC~SA{{[L],...,1,...}}, pot~|a⟩×|a⟩×[L]×[M]→[E],
    energy~[E], lhs~[L], rhs~[L]; B~[Tesla] magnetic field, μ~[m] He* mass
    Output: sol(R) [where R ∈ (lhs,rhs)] ~ IC"""
function solver(lookup, IC, ϵ, lhs, rhs; B=0.0u"T", μ=0.5*4.002602u"u")
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
        V = zeros(ComplexF64,n,n)u"hartree" # initialise
        for i=1:n, j=1:n
            V[i,j] = H_rot(lookup[i],lookup[j], x*1u"bohr", μ) # rotational
            V[i,j]+= H_el(lookup[i],lookup[j], x*1u"bohr") # electronic
            V[i,j]+= C_sd[i,j]*H_sd_radial(x*1u"bohr") # spin-dipole
            #imaginary ionization potential
            Γ(i,j,x) = (i==j && lookup[i].S∈[0,1]) ? 0.3*exp(-x/1.086) : 0.0
            V[i,j]-= im*Γ(i,j,x)*1u"hartree"
            #TODO hyperfine interaction
            V[i,j] += H_zee(lookup[i],lookup[j],B) # zeeman
        end
        V⁰=austrip.(V) # strip units from V
        M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matix
        D = SMatrix{2*n,2*n}([0*I I
                              M 0*I])
        D*u # ⃗u' = D . ⃗u
    end
    # strip units from IC
    IC⁰ = austrip.(complex.(IC))
    # solve
    prob=ODEProblem(TISE,IC⁰,(lhs⁰,rhs⁰))
    sol_unitless=solve(prob,Tsit5(),reltol=1e-12)
    # add units back
    units = SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol = x -> sol_unitless(austrip(x)).*units
    return sol
end

using SpecialFunctions, LinearAlgebra
""" Solver for 𝐊 matrix (following Mies eqn (3.8)), where 𝐅=𝐉-𝐍𝐊.
    Inputs:
        R~[L] the asymptotic radial distance to match K at;
        𝐅 = [𝐆; 𝐆'] where 𝐆 is matrix of wavefunctions for
    different initial conditions evaluated at R;
        𝐤 ~ [L]⁻¹ vector of wavenumber for each channel at rhs;
    **Note: this function assumes it is only being given open channels**
        𝐥 vector of l quantum numbers for each channel.
    **Note: eval, 𝐤, 𝐥 must share the same ordering/number of channels**
    Output:
        𝐊 ~ n×n matrix (n=number of channels considered)"""
function K_matrix(R, 𝐅, 𝐤, 𝐥)
    # match for A, B where G=J.A-N.B
    #construct G, G' matrices to match for A, B with
    n=Int(size(𝐅,1)/2) # n = number of channels. Assumes sol in above form.
    @assert size(𝐅,2) == n "solution matrix not of shape 2n × n"
    G, G⁻ = austrip.(copy(𝐅[1:n,1:n])), copy(𝐅[n+1:2*n,1:n])
    # solve for A,B
    A, B = zeros(ComplexF64,n,n), zeros(ComplexF64,n,n) # initialise
    for i in 1:n, j in 1:n
        # construct [jᵢ nᵢ; jᵢ' nᵢ'] matrix, here called [bj -bn; bj⁻ -bn⁻]
        # expressions for derivatives (⁻) calculated using Mathematica
        # function form from Mies (A2) which ≡ Cocks et al (59)
        k=𝐤[i]
        l=𝐥[i]
        bj=austrip(sqrt(k)*R)*sphericalbesselj(l,k*R)
        bj⁻=austrip(sqrt(k))*((l+1)*sphericalbesselj(l,k*R)
            -k*R*sphericalbesselj(l+1,k*R))
        bn=austrip(sqrt(k)*R)*sphericalbessely(l,k*R)
        bn⁻=austrip(sqrt(k))*((l+1)*sphericalbessely(l,k*R)
            -k*R*sphericalbessely(l+1,k*R))
        Gᵢⱼ, G⁻ᵢⱼ = G[i,j], G⁻[i,j]
        AB = [bj -bn; bj⁻ -bn⁻]\[Gᵢⱼ; G⁻ᵢⱼ] # AB≡[Aᵢⱼ; Bᵢⱼ], solve J;-N*AB=G;G⁻
        A[i,j], B[i,j] = AB
    end
    𝐊 = B*inv(A)
    return 𝐊
end

"""Boundary condition matching via SVD
    Inputs: AL, BCs on LHS;
            AR, wavefunction solution to AL at matching location;
            BL, wavefunction solution to BR at matching location;
            BR, BCs on RHS;
            isOpen::[Bool], vector isOpen[i]⟺"channel i is open"
            tol_ratio=1e-10 is the ratio used to define a zero singular value
    ***Assumes same channel ordering in AL,AR,BL,BR,isOpen***
    Output: F, 2Nₒ×Nₒ matrix of valid (i.e. open channel) wavefunction solutions
    evaluated at the RHS"""
function F_matrix(AL,AR,BL,BR,isOpen; tol_ratio=1e-10)
    # check on dimensions of inputs
    @assert size(AL)==size(AR) "AL and AR have unlike dimensions"
    @assert size(BL)==size(BR) "BL and BR have unlike dimensions"
    @assert size(AL,1)==size(BL,1) "AL/AR and BL/BR have unlike numbers of rows"
    @assert size(AL,1)==2*length(isOpen) "Number of rows of AL/AR/BL/BR not equal
        to 2* number of rows of isOpen vector"
    # numbers of channel, for reference
    N = length(isOpen) # N channels
    Nₒ = size(BL,2)-N # Nₒ open channels
    # take SVD
    x = svd(austrip.([AR -BL]), full=true) # the SVD object
    Σ, V = x.S, x.V # extract singular values and V matrix
    # find [C; D] where [A -B]*[C;D]=0 by extracting cols of V corresponding
    # to singular values of zero
    #=
    zero_cols=[] # stores indices of singular values of zero
    zero_scale=maximum(abs.(Σ))*tol_ratio # reference scale for checking ≈0
    for i in 1:length(Σ)
        if isapprox(Σ[i], 0, atol=zero_scale) # singular value ≊ zero
            push!(zero_cols, i)
        end
    end
    CD = V[:,zero_cols] # σ≈0 cols of V are the cols of [C;D]
    =#
    CD = V[:,(end-Nₒ+1):end] # 4/09/20 cols of V matching to the zero part of Σ
    # sanity check for linear combinations
    @assert size(CD,1)==2*N+Nₒ "[C; D] doesn't have 2*N+N₀ rows"
    @assert size(CD,2)==Nₒ "[C; D] doesn't have Nₒ columns"
    C = CD[1:N,:]
    D = CD[(N+1):end,:]
    # forming F
    F = BR*D
    F=F[vcat(isOpen,isOpen),:] # taking only open wavefunctions and derivatives
end



################################################################################
# Test functions
################################################################################

# test function for solver - runs and plots first channel wavefunction
function test_solver(;lmax=0,B=0u"T")
    lhs, rhs = 3.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=SmS_lookup_generator(lmax)
    n=length(lookup)
    # construct ICs
    IC=SMatrix{2*n,n}([fill(0.0u"bohr",n,n)
                       I])
    println("Starting to solve for wavefunctions, lmax=$lmax")
    @time begin sol=solver(lookup, IC, ϵ, lhs, rhs,B=B)
    end
    println("Plotting...")
    Rs=LinRange(lhs,rhs,1000)
    vals = getindex.(sol.(Rs),n,n)
    plot(austrip.(Rs), austrip.(vals),title="It works! Plotting wavefn of last channel", legend=false)
end

# unit test for K_matrix. Should produce a scattering length of 7.54 nm
# in agreement with Przybytek
function test_K_matrix(;lmax=0, ϵ=1e-9u"hartree", μ=0.5*4.002602u"u",
    lhs=3.0u"bohr", rhs=1000u"bohr")
    println("Starting test_K_matrix")
    lookup=SmS_lookup_generator(lmax)
    n=length(lookup)
    # construct ICs
    IC=SMatrix{2*n,n}([fill(0.0u"bohr",n,n)
                       I])
    # solver for wavefunctions
    println("solving for wavefunctions")
    sol=solver(lookup,IC,ϵ,lhs,rhs,μ=μ)
    eval=sol(rhs)
    # solve 𝐤 vector for K matrix solver
    println("Producing k vector")
    𝐤=fill(0.0u"bohr^-1",n)
    for i in 1:n
        γ = lookup[i] # channel
        R∞ = Inf*1u"bohr"
        V∞ = H_el(γ,γ,R∞) + H_sd_coeffs(γ,γ)*H_sd_radial(R∞) + H_rot(γ,γ,R∞,μ)
        𝐤[i] = sqrt(2*μ*(ϵ-V∞))/1u"ħ"
    end
    # 𝐥 vector for K matrix solver
    println("Producing 𝐥 vector")
    𝐥=fill(0,n)
    for i in 1:n
        𝐥[i]=lookup[i].l
    end
    println("Passing to K_matrix function")
    𝐊=K_matrix(rhs, eval, 𝐤, 𝐥)
    𝐒=(I+im*𝐊)*inv(I-im*𝐊)
    #return uconvert(u"nm", sqrt(pi*abs(1-𝐒[n,n])^2/𝐤[n]^2/(4*pi)))
    return 𝐒
end

# unit test for F_matrix. Should produce a matrix of wavefunction solutions
function test_F_matrix(;lmax=0, ϵ=1e-12u"hartree", μ=0.5*4.002602u"u",
    lhs=3.0u"bohr", mid=100.0u"bohr", rhs=1000.0u"bohr")
    println("Running test_F_matrix")
    println("Initialising AL, BR, isOpen")
    # arbitrary sample AL, AR
    lookup=SmS_lookup_generator(lmax)
    N=length(lookup)
    # construct isOpen
    isOpen=fill(true,N)
    for i in 1:N
        γ = lookup[i] # channel
        R∞ = Inf*1u"bohr"
        V∞ = H_el(γ,γ,R∞) + H_sd_coeffs(γ,γ)*H_sd_radial(R∞) + H_rot(γ,γ,R∞,μ)
        ksq = 2*μ*(ϵ-V∞)/1u"ħ^2"
        isOpen[i] = austrip(ksq) >= 0 ? true : false # k² > 0 for open channels
    end
    # construct BCs
    AL=SMatrix{2*N,N}([fill(0.0u"bohr",N,N)
                       I])
    BR = let
        Nₒ=count(isOpen)
        BFL = SMatrix{2*N,N}([fill(0.0u"bohr",N,N);I])
        BFR = SMatrix{2*N,Nₒ}([Matrix(Diagonal(ones(N))[:,isOpen]u"bohr");zeros(N,Nₒ)])
        [BFL BFR]
    end
    # solve for solutions
    println("Solving for AR and BL")
    AR = solver(lookup, AL, ϵ, lhs, mid)(mid)
    BL = solver(lookup, BR, ϵ, rhs, mid)(mid)
    #=#Bug fixing 4/09 normlaisation
    AR = AR./maximum(abs.(austrip.(AR)),dims=1)
    BL = BL./maximum(abs.(austrip.(BL)),dims=1)=#
    # see if F_matrix runs
    println("Passing to F_matrix")
    𝐅=F_matrix(AL,AR,BL,BR,isOpen)
    println("Finished test_F_matrix")
    return 𝐅
end

# combined tests for F and K functions. Should produce a Quintet scattering
# length of 7.54nm in agreement with Przybytek
function test_K_and_F(;lmax=0, ϵ=1e-12u"hartree", μ=0.5*4.002602u"u",
    lhs=3.0u"bohr", mid=100.0u"bohr", rhs=1000.0u"bohr")
    # calculate F
    println("Calculating F")
    𝐅=test_F_matrix(lmax=lmax,ϵ=ϵ,μ=μ,lhs=lhs,mid=mid,rhs=rhs)
    # calculate 𝐤 vector for K_matrix()
    lookup=SmS_lookup_generator(lmax)
    println("Calculating k vector")
    n=length(lookup)
    𝐤=fill(0.0u"bohr^-1",n)
    for i in 1:n
        γ = lookup[i] # channel
        R∞ = Inf*1u"bohr"
        V∞ = H_el(γ,γ,R∞) + H_sd_coeffs(γ,γ)*H_sd_radial(R∞) + H_rot(γ,γ,R∞,μ)
        𝐤[i] = sqrt(2*μ*(ϵ-V∞))/1u"ħ"
    end
    # 𝐥 vector for K matrix solver
    println("Constructing l vector")
    𝐥=fill(0,n)
    for i in 1:n
        𝐥[i]=lookup[i].l
    end
    println("Calculating K")
    𝐊=K_matrix(rhs,𝐅,𝐤,𝐥)
    𝐒=(I+im*𝐊)*inv(I-im*𝐊)
    return uconvert(u"nm", sqrt(pi*abs(1-𝐒[n,n])^2/𝐤[n]^2/(4*pi)))
end
