#= Multichannel equations solver. Uses symmetrised |a⟩ basis states and all
of the interactions except for the Zeeman interaction.
All states up to and including the lmax parameter are used
F matrix matches boundary conditions, K matrix matches to bessel functions
σ_matrix uses the other functions to produce cross sections from high-level inpt
Description last updated 9/09/20=#

module CrossSections
export σ_output, σ_matrix

using HalfIntegers, LinearAlgebra, OrdinaryDiffEq, WignerSymbols
using Unitful, UnitfulAtomic
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using Interactions, StateStructures

""" Multichannel TISE solver.
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
            #imaginary ionization potential width from Garrison et al 1973
            Γ(i,j,x) = (i==j && lookup[i].S∈[0,1]) ? 0.3*exp(-x/1.086) : 0.0
            V[i,j]-= im/2*Γ(i,j,x)*1u"hartree" # Cocks et al (2019)
            #TODO hyperfine interaction
            V[i,j] += H_zee(lookup[i],lookup[j],B) # zeeman
        end
        V⁰=austrip.(V) # strip units from V
        M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matix
        D = ([0*I I
              M 0*I])
        D*u # ⃗u' = D . ⃗u
    end
    # strip units from IC
    IC⁰ = austrip.(complex.(IC))
    # solve
    prob=ODEProblem(TISE,IC⁰,(lhs⁰,rhs⁰))
    sol_unitless=solve(prob,Tsit5(),reltol=1e-10,saveat=[rhs⁰])
    # add units back
    units = vcat(fill(1.0u"bohr",n),fill(1.0,n))
    sol = x -> sol_unitless(austrip(x)).*units
    return sol
end

using SpecialFunctions
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

# structure for holding cross sections and the inputs that produced them
struct σ_output
    σ # the matrix of cross sections
    γ_lookup :: Array{γ_ket,1}
    ϵ :: Unitful.Energy # energy input
    B :: Unitful.BField # B field strength input
    lmax :: Int # lmax input
end

""" Calculates matrix of cross sections between γ, only summing over l and mₗ
    Input: ϵ~[E], B~[Tesla], lmax
    Output: σ_output containing: σ where σ[i,j]=σ(γ_lookup[j]→γ_lookup[i]),
    γ_lookup describing the γ_kets involved, ϵ input, B input, lmax input"""
function σ_matrix(ϵ::Unitful.Energy,B::Unitful.BField,lmax::Int;
    lhs::Unitful.Length=3.0u"bohr", mid::Unitful.Length=100.0u"bohr",
    rhs::Unitful.Length=1000.0u"bohr",μ::Unitful.Mass=0.5*4.002602u"u")
    @info "Starting σ_matrix"
    # lookup vector of |SmS⟩ states to be considered
    lookup=SmS_lookup_generator(lmax)
    N=length(lookup)
    # construct isOpen. simultaneously construct 𝐤 and 𝐥 for K calculation
    isOpen=fill(true,N)
    𝐤Open=Array{typeof(1.0u"bohr^-1")}([])
    𝐥Open=Array{Int,1}([])
    for i in 1:N
        ϕ = lookup[i] # channel
        R∞ = Inf*1u"bohr"
        V∞ = H_el(ϕ,ϕ,R∞) + H_sd_coeffs(ϕ,ϕ)*H_sd_radial(R∞) + H_rot(ϕ,ϕ,R∞,μ) + H_zee(ϕ,ϕ,B)
        ksq = 2*μ*(ϵ-V∞)/1u"ħ^2"
        if austrip(ksq) >= 0
            isOpen[i] = true
            push!(𝐤Open,uconvert(u"bohr^-1",sqrt(ksq)))
            push!(𝐥Open,ϕ.l)
        else # ksq < 0 ⟺ closed channel
            isOpen[i] = false
        end
    end
    Nₒ=count(isOpen); Nₒ==0 && return("No open channels!")
    @assert length(findall(isOpen))==length(𝐤Open)==Nₒ "number of
    open channels disagrees between isOpen, 𝐤Open and 𝐥Open" # sanity check
    # construct BCs
    @info "constructing BCs"
    AL=[fill(0.0u"bohr",N,N); I]
    BR = let
        Nₒ=count(isOpen)
        BFL = [fill(0.0u"bohr",N,N); I]
        BFR = [Matrix(Diagonal(ones(N))[:,isOpen]u"bohr"); zeros(N,Nₒ)]
        [BFL BFR]
    end
    # solve for inividual BCs
    @info "solving for AR, BL"
    AR = solver(lookup, AL, ϵ, lhs, mid,B=B,μ=μ)(mid)
    BL = solver(lookup, BR, ϵ, rhs, mid,B=B,μ=μ)(mid)
    # find wavefunction satisfying both BCs only including open channels
    @info "solving for F"
    𝐅 = F_matrix(AL,AR,BL,BR,isOpen)
    # match to bessel functions to find K matrix
    @info "solving for K"
    𝐊 = K_matrix(rhs,𝐅,𝐤Open,𝐥Open)
    @assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
    𝐒 = (I+im*𝐊)*inv(I-im*𝐊) # Scattering matrix
    𝐓 = I-𝐒 # transition matrix
    𝐓sq= abs2.(𝐓) # |Tᵢⱼ|² for use in calculating cross sections
    # initialise γ states used for cross sections
    @info "producing γ states"
    γ_lookup=unique(γ_ket_convert.(lookup[isOpen]))
    nᵧ=length(γ_lookup)
    # create k²ᵧ vector used to calculate cross sections
    𝐤²ᵧ=let
        k²(γ::γ_ket)=uconvert(u"bohr^-2",2*μ*(ϵ-H_zee(γ,γ,B))/1u"ħ^2") # only Zeeman at R=∞
        k².(γ_lookup)
    end
    # initialise σ array
    @info "constructing σ array"
    𝛔=zeros(nᵧ,nᵧ)u"bohr^2"
    # fill in entries σᵢⱼ. row i = output = γ_, column j = input = γ
    for i in 1:nᵧ, j in 1:nᵧ
        γ_, γ = γ_lookup[i], γ_lookup[j]
        prefac = π/𝐤²ᵧ[j] # k²ᵧ from incoming state
        sum = 0 #initialise sum
        for m in 1:Nₒ, n in 1:Nₒ # row m output, col n input
             γ_ket_convert(lookup[isOpen][m])==γ_ || continue # output states match
             γ_ket_convert(lookup[isOpen][n])==γ || continue # input states match
             sum += 𝐓sq[m,n]
         end
         𝛔[i,j]=prefac*sum
     end
     return σ_output(𝛔,γ_lookup,ϵ,B,lmax)
end

end # module
