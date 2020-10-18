#Tested working 8/10/20, git committed in working state
# sketch for producing artificial Feshbach resonance

using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots, Plots.PlotMeasures
using BSON, Dates

using CrossSections: F_matrix, K_matrix
using OrdinaryDiffEq

const G = 1e-4u"T" # Gauss unit of magnetic flux density

const Smat_dir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\resultsA\Smat"
const gampwcs_dir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\resultsA\gampwcs"
const Ics_dir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\resultsA\Ics"

const L=10.0u"bohr"
const ϵ=1e-8u"hartree"
const μ=0.5*4.002602u"u"

# parameters for ICs/matching
const lhs=1e-5*L; const mid=0.75*L; const rhs=(1+1e-3)*L; const rrhs=100*L

# artificial ket structure
struct art_ket
    x :: Int
    l :: Int
    art_ket(x) = x in [0,1] ? new(x,0) : error("x not in [0,1]")
    art_ket(x,l) = x in [0,1] && l==0 ? new(x,l) : error("x not in [0,1] or l != 0")
end
struct art_output
    S # Scattering matrix
    σ # cross section
    lookup :: Array{art_ket,1}
    ϵ :: Unitful.Energy
    C :: Unitful.Energy
end

# artificial potentials
function Vcoup(bra::art_ket, ket::art_ket, R::Unitful.Length)
    bra != ket || return 0.0u"hartree" # unlike kets
    R < L/2 || return 0.0u"hartree" # within range
    return 1e-9u"hartree"#1e-7u"hartree"
end
function Vbase(bra::art_ket, ket::art_ket, R::Unitful.Length, C::Unitful.Energy)
    bra==ket || return 0.0u"hartree" # no coupling
    ket.x==0 && return 0.0u"hartree" # zero potential for x=0
    if R < 0u"bohr"
        return 1.0u"hartree"
    elseif R < L # well region
        return C
    else
        return 1.0u"hartree"
    end
end

# Necessary C for a Feshbach resonance (in theory)
necessary_C(ϵ::Unitful.Energy;n=1)=auconvert(ϵ - n^2*π^2*1u"ħ^2"/(2*μ*L^2))


###############################DE Solver #######################################

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

""" Artificial resonance multichannel TISE solver.
    Input: lookup vector, IC~{[L],...,1,...}, pot~|a⟩×|a⟩×[L]×[M]→[E],
    energy~[E], lhs~[L], rhs~[L]; B~[Tesla] magnetic field, μ~[m] He* mass
    Output: sol(R) [where R ∈ (lhs,rhs)] ~ IC"""
function art_solver(lookup, IC, ϵ, lhs, rhs; C=0u"hartree", μ=0.5*4.002602u"u")
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
        V = zeros(ComplexF64,n,n)u"hartree" # initialise
        for i=1:n, j=1:n
            V[i,j] = Vcoup(lookup[i],lookup[j],x*1u"bohr")
            V[i,j]+= Vbase(lookup[i],lookup[j],x*1u"bohr",C)
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
    callback=CreateRenormalisedCallback()
    sol_unitless=solve(prob,Tsit5(),reltol=1e-10,save_start=true,save_end=true,save_everystep=false,dense=false,
    callback=callback) #TODO added save_start 22/9/20, need to test for lmax>0
    # add units back
    units = vcat(fill(1.0u"bohr",n),fill(1.0,n))
    sol = x -> sol_unitless(austrip(x)).*units
    return sol
end

######################################S Matrix##################################

""" Calculates S and T matrices of artificial kets
    Input: ϵ~[E], C~[E]
    Output: art_output containing: S where S[i,j]=S(lookup[j]→lookup[i]),
    lookup describing the kets involved, ϵ input, C input"""
function art_ST_matrix(ϵ::Unitful.Energy,C::Unitful.Energy,
    lhs::Unitful.Length, mid::Unitful.Length,
    rhs::Unitful.Length, rrhs::Unitful.Length)
    #println("Starting art_S_matrix")
    #println("Generating isOpen, kOpen, lOpen")
    # lookup vector of |SmS⟩ states to be considered
    lookup=[art_ket(0),art_ket(1)]
    N=length(lookup)
    # construct isOpen. simultaneously construct 𝐤 and 𝐥 for K calculation
    isOpen=fill(true,N)
    𝐤Open=Array{typeof(1.0u"bohr^-1")}([])
    𝐥Open=Array{Int}([])
    for i in 1:N
        ϕ = lookup[i] # channel
        R∞ = Inf*1u"bohr"
        V∞ = Vcoup(ϕ,ϕ,R∞) + Vbase(ϕ,ϕ,R∞,C)
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
    #println("Constructing BCs")
    AL=[fill(0.0u"bohr",N,N); I]
    BR = let
        BFL = [fill(0.0u"bohr",N,N); I]
        BFR = [Matrix(Diagonal(ones(N))[:,isOpen]u"bohr"); zeros(N,Nₒ)]
        [BFL BFR]
    end
    # solve for inividual BCs
    #println("Solving AR")
    Asol = art_solver(lookup, AL, ϵ, lhs, mid,C=C,μ=μ)
    AL, AR = Asol(lhs), Asol(mid)
    #println("Solving BL")
    Bsol = art_solver(lookup, BR, ϵ, rhs, mid,C=C,μ=μ)
    BL, BR = Bsol(mid), Bsol(rhs)
    # find wavefunction satisfying both BCs only including open channels
    #println("Matching for F")
    𝐅 = F_matrix(AL,AR,BL,BR,isOpen)
    # solve matched wavefunction out to rrhs TODO added 22/9/20
    #println("Solving F")
    𝐅 = art_solver(lookup[isOpen],𝐅,ϵ,rhs,rrhs,C=C,μ=μ)(rrhs)
    # match to bessel functions to find K matrix
    #println("Matching for K")
    𝐊 = K_matrix(rrhs,𝐅,𝐤Open,𝐥Open)
    @assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
    #println("Producing S, T")
    𝐒 = (I+im*𝐊)*inv(I-im*𝐊) # Scattering matrix
    @assert length(𝐒)==1 "length(𝐒) is not 1"
    k = sqrt(2*μ*ϵ)/1u"ħ"
    𝐓 = auconvert.(π/k^2 .* abs2.(𝐒))
    return art_output(𝐒[1], 𝐓[1], lookup[isOpen], ϵ, C)
end

using ProgressMeter
# generate different art_outputs for constant ϵ and varying C
function diffE_gen_art(Cmin::Unitful.Energy, Cmax::Unitful.Energy, n::Integer)
    data=[]
    Cs=LinRange(Cmin,Cmax,n)
    @showprogress for C in Cs
        push!(data,art_ST_matrix(ϵ,C,lhs,mid,rhs,rrhs))
    end
    data
end

nec=necessary_C(ϵ)
#Cmin, Cmax, n = (1-1e-2)*nec, (1+1e-2)*nec, 1000
n=10_000
Ctarget = -1.3495204361066107e-5u"hartree"
Cmin,Cmax = Ctarget-1e-12u"hartree",Ctarget+1e-12u"hartree"

data=diffE_gen_art(Cmin, Cmax, n)
@assert all(x->length(x.lookup)==1,data) "not all data have 1 open channel"
Ss=(x->x.S).(data)
Cs=LinRange(Cmin, Cmax, n)
Ps=angle.(Ss)./2 # phases #S = exp(2iδₗ)
maxindex=findmax(Ps)[2]
plt=plot(austrip.(Cs.-Ctarget),Ps,yticks=-π:π:π,
    grid=false,legend=false,
    bottom_margin=5mm,left_margin=5mm,top_margin=5mm,
    xlabel="C - C₀ (Eh)", ylabel="Phaseshift (rad)",
    title="C₀=$Ctarget")
#vline!([austrip(nec)])

bound_n=austrip(sqrt(2*μ*(ϵ-Ctarget))*L/(π*1u"ħ"))
