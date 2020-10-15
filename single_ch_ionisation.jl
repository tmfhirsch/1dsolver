# Tested working 12/10, with the sign change to th eionisation potnetial in solver()
# Test ionisation diverges as k→0 for single channel S=0

using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots, Plots.PlotMeasures

using CrossSections: F_matrix, K_matrix, solver
using OrdinaryDiffEq
using Interactions, StateStructures

const G = 1e-4u"T" # Gauss unit of magnetic flux density
const μ=0.5*4.002602u"u"

# parameters for ICs/matching
const lhs=3e0u"bohr"; const mid=5e1u"bohr"; const rhs=2e2u"bohr"; const rrhs=1e3u"bohr"

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
######################################S Matrix##################################
struct sin_output
    S # S matrix
    T # transition matrix
    I # ionisation cross section
    lookup :: Array{SmS_ket,1}
    ϵ :: Unitful.Energy
    B :: Unitful.BField
end

function sin_S_matrix(ϵ::Unitful.Energy,B::Unitful.BField,lmax::Int,
    lhs::Unitful.Length, mid::Unitful.Length,
    rhs::Unitful.Length, rrhs::Unitful.Length;
    μ::Unitful.Mass=0.5*4.002602u"u")
    # lookup vector of |SmS⟩ states to be considered
    lookup=filter(x->x.S==0,SmS_lookup_generator(lmax))
    N=length(lookup)
    # construct isOpen. simultaneously construct 𝐤 and 𝐥 for K calculation
    isOpen=fill(true,N)
    𝐤Open=Array{typeof(1.0u"bohr^-1")}([])
    𝐥Open=Array{Int,1}([])
    for i in 1:N
        ϕ = lookup[i] # channel
        if ϕ.l>0 #TODO added 25/9/20, artificially zeroes out l>0 channels
            isOpen[i] = false
            continue
        end
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
    AL=[fill(0.0u"bohr",N,N); I]
    BR = let
        BFL = [fill(0.0u"bohr",N,N); I]
        BFR = [Matrix(Diagonal(ones(N))[:,isOpen]u"bohr"); zeros(N,Nₒ)]
        [BFL BFR]
    end
    # solve for inividual BCs
    Asol = solver(lookup, AL, ϵ, lhs, mid,B=B,μ=μ)
    AL, AR = Asol(lhs), Asol(mid)
    Bsol = solver(lookup, BR, ϵ, rhs, mid,B=B,μ=μ)
    BL, BR = Bsol(mid), Bsol(rhs)
    # find wavefunction satisfying both BCs only including open channels
    𝐅 = F_matrix(AL,AR,BL,BR,isOpen)
    # solve matched wavefunction out to rrhs TODO added 22/9/20
    𝐅 = solver(lookup[isOpen],𝐅,ϵ,rhs,rrhs,B=B,μ=μ)(rrhs)
    # match to bessel functions to find K matrix
    𝐊 = K_matrix(rrhs,𝐅,𝐤Open,𝐥Open)
    @assert size(𝐊)==(Nₒ,Nₒ) "𝐊 is not Nₒ×Nₒ"  # want sq matrix of Nₒ channels
    𝐒 = (I+im*𝐊)*inv(I-im*𝐊) # Scattering matrix
    # calculate ionisation cross sectioins
    𝐤²=let
        k²(a::SmS_ket)=uconvert(u"bohr^-2",2*μ*(ϵ-H_zee(a,a,B))/1u"ħ^2") # only Zeeman at R=∞
        k².(lookup)
    end
    σT = [π/𝐤²[i]*abs2(𝐒[i]) for i in 1:length(lookup)]
    σPI = [π/𝐤²[i]*(1-abs2(𝐒[i])) for i in 1:length(lookup)]
    return sin_output(𝐒, σT, σPI, lookup[isOpen], ϵ, B)
end

using ProgressMeter
# generate different art_outputs for constant ϵ and varying C
function sin_diffE_data(Emin_exp,Emax_exp,n::Integer,lmax::Integer;B=0u"T",
    lhs=lhs,mid=mid,rhs=rhs,rrhs=rrhs)
    Es=exp10.(LinRange(Emin_exp,Emax_exp,n))u"hartree" # energies
    data=[]
    println("Emin=1e$(Emin_exp) Eₕ, Emax=1e$(Emax_exp) Eₕ, n=$n, lmax=$lmax, B=$B")
    @showprogress for E in Es
        output=sin_S_matrix(E,B,lmax,lhs,mid,rhs,rrhs)
        push!(data,output)
    end
    data
end


Emin, Emax = -12, -8; n=50; B=0u"T"

data=sin_diffE_data(Emin, Emax, n, 0)
@assert all(x->length(x.lookup)==1,data) "not all data have 1 open channel"
Es, Ts, Is = (x->x.ϵ).(data), (x->x.T[1]).(data), (x->x.I[1]).(data)
ks = auconvert.((x->sqrt(2*μ*x)/1u"ħ").(Es))
vs = uconvert.(u"bohr/s",(1u"ħ"/μ).*ks)

kIplt=plot(ustrip.(ks),ustrip.(Is),yscale=:log10,xscale=:log10,
    xlabel="k (a₀⁻¹)", ylabel="σ (a₀²)",
    linewidth=2,legend=false,grid=false,
    left_margin=5mm,bottom_margin=5mm)
kRplt=plot(ustrip.(ks),ustrip.(uconvert.(u"bohr^3/s",vs.*Is))./1e15,xscale=:log10,#yscale=:log10,
    xlabel="k (a₀⁻¹)", ylabel="σv (10¹⁵ a₀³ s⁻¹)",
    linewidth=2,legend=false,grid=false,
    left_margin=5mm,bottom_margin=5mm)
