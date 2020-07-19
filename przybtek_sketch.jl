#= Playing around with matching coefficients with Przybytek quintet potential,
with the aim of producing a scattering length

Description last updated 14/07/20=#

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using Revise
using OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays, SpecialFunctions, Statistics
using Unitful, UnitfulAtomic
using Przybytek: przybytek


"""
TISE solver for IC of (u,u')=(0,1)
All inputs in units, or assumed to be in atomic units if not.
Inputs: energy, limits; V(R)=przybytek, m=He₂ mass
Output: (u,u')~(1,a₀⁻¹) as fn of R~a₀"""
function rhs_solver(k, # wavenumber [L]⁻¹
                    l::Int; # angular momentum
                    stapt=1.0u"bohr", # location of IC [L]
                    endpt=100.0u"bohr", # RHS to be solved to [L]
                    U=przybytek, # interatomic potential [L]->[E]
                    μ=0.5*4.002602u"u", # He₂ reduced mass [M]
                    reltol=1e-10 #relative tolerance for DE solver
                    )
    ϵ=auconvert(k^2*1u"ħ^2"/(2*μ)) # E=ħ²k²/2μ
    # add centrifugal potential
    UL = R -> auconvert(U(R)+l*(l+1)u"ħ^2"/(2*μ*R^2))
    # check units are dimensionally correct
    @assert dimension(2*μ*(UL(stapt)-ϵ)/1u"ħ^2")==dimension(1u"m^-2") "TISE units mismatch"
    # convert and strip units before TISE
    V = x -> austrip(UL((x)u"bohr")) # unitless -> unitless potential
    ϵ⁰ = austrip(ϵ) # strip energy
    stapt⁰, endpt⁰ = austrip(stapt), austrip(endpt) # strip start/end points
    μ⁰ = austrip(μ) # strip mass
    ħ⁰ = austrip(1u"ħ") # strip ħ
    # TISE solver
    function TISE(u,p,x) # TISE gives (u,u')'=f((u,u'))
        du = u[2] # u'(x)≡u'(x)
        dd = 2*μ⁰*(V(x)-ϵ⁰)*u[1]/ħ⁰^2 # (u'(x))'=2m(V(x)-E)u(x)/ħ^2 {𝐋⁻²}
        SVector{2}([du,dd]) # (u,u')'
    end
    IC = SVector{2}([0.0, 1.0]) # (u,u')=(0,1)
    prob=ODEProblem(TISE,IC,(stapt⁰,endpt⁰))
    sol_unitless=solve(prob,Tsit5(),reltol=1e-10) # comes out of DE sans units
    sol_input = x -> sol_unitless(austrip(x)) # add unit input
    sol = x -> sol_input(x).*SA[1,1u"bohr^-1"] # add unit output
    return sol
end
#=
#Test plot of rhs_solver
plot_k=1e-5u"bohr^-1"
plot_stapt=1.0u"bohr"
plot_endpt=20000.0u"bohr"
sol=rhs_solver(plot_k,stapt=plot_stapt,endpt=plot_endpt)
Rs=LinRange(plot_stapt,plot_endpt,1000)
us=[i[1] for i in sol.(Rs)]
plot(ustrip.(Rs), ustrip.(us), legend=false, title="$plot_ϵ")
#savefig("prz_wfn_ϵ-$plot_ϵ.png")
=#

"""
Spherical bessel functions
"""
j(l,x)=sphericalbesselj(l,x)
n(l,x)=sphericalbessely(l,x)

"""
Returns (Aₗ,Bₗ)(R) that match wavefunction to spherical bessel functions
Inputs: (u,u')(R) solution ~ (1,[L]⁻¹), k ~ [L⁻¹], l
Outputs: (Aₗ,Bₗ)(R)
"""
function matchAB(sol, #(u,u')(R) ~ [L]->(1,[L]⁻¹)
                 k, # [L]⁻¹
                 l::Int # ang. momentum
                 )
    @assert l>=0 "Not a valid angular momentum"
    # define matrix of matching functions M(R) ~ [L] -> [1 1; [L]⁻¹ [L]⁻¹]
    if l==0 #𝑠-wave
        M = R -> [sin(k*R)   cos(k*R); # sin(kR), cos(kR) & derivs
                  k*cos(k*R) -k*sin(k*R)]
    else
        M = R -> [j(l,k*R)                    n(l,k*R); #jₗ(kR), nₗ(kR) & derivs
                  (l/R)*j(l,k*R)-k*j(l+1,k*R) (l/R)*n(l,k*R)-k*n(l+1,k*R)]
    end
    AB = R -> ustrip.(M(R))\ustrip.(sol(R)) # AB matched at R ~ [L]->[1,1]
    return AB
end

#=
# Investigating B(R)/A(R)
lhs=1.0u"bohr"
rhs=1e8u"bohr"
k=1e-5u"bohr^-1"
l=1
no_pts=1000
zero_pot(R)=0.0u"hartree"
wavefn=rhs_solver(k,l,stapt=lhs,endpt=rhs,reltol=1e-10)
ABfn=matchAB(wavefn,k,l,lhs,rhs)
Rs=LinRange(lhs,rhs,no_pts)
ABs=ABfn.(Rs)
As=getindex.(ABs,1)
Bs=getindex.(ABs,2)
Aplt=plot(austrip.(Rs),As,label="A")
Bplt=plot(austrip.(Rs),Bs,label="B")
#plot(Aplt,Bplt,layout=(2,1),title="l=$l, k=$k")
plot(austrip.(Rs),Bs./As,xlabel="R (a₀)",ylabel="B(R)/A(R)",title="l=$l,k=$k",legend=false)
#tan_δ₀ = mean((Bs./As)[Int(floor(no_pts/2)):end])
#σ=std((Bs./As)[Int(floor(no_pts/2)):end])
#plot(austrip.(Rs),austrip.(getindex.(wavefn.(Rs),1)))
=#


"""
Finds convergence of B(R)/A(R)
Input: AB(R) function; no. grid pts, median bubble, convergence tolerance,
warning 10^index, break10^index
Output: tan(δ₀) within tolerance, or an error if 10^break reached
"""
function BoA_lim(AB; no_pts=1000::Int, bub=1.0, tol=1e-8, warn=6::Int, stop=8::Int)
    BoA(R) = AB(R)[2]/AB(R)[1] #B/A -> tanδₗ
    i=1 # Index for R∈(10ⁱ,10ⁱ⁺¹)
    while true
        Rs=LinRange(10^i,10^(i+1),no_pts)u"bohr" # linear grid of R values
        BoAs=BoA.(Rs) # AB values on the grid
        med = median(BoAs)
        fil = x -> abs(x-med)<bub # filter within <bub> of the median
        BoAs=filter(fil, BoAs)
        if std(BoAs) < tol
            return mean(BoAs)
            break
        end
        if i == warn
            @warn "B(R)/A(R) hasn't converged by 10^$warn a₀"
        elseif i == stop
            @error "B(R)/A(R) did not converge by 10^$break a₀"
        end
        i += 1
    end
end

# Investigating B(R)/A(R)
lhs=1.0u"bohr"
rhs=1e8u"bohr"
k=1e-5u"bohr^-1"
l=2
wavefn=rhs_solver(k,l,stapt=lhs,endpt=rhs,reltol=1e-10)
ABfn=matchAB(wavefn,k,l)
str="l=$l, tan(δₗ) = "*string(BoA_lim(ABfn))
@info str
