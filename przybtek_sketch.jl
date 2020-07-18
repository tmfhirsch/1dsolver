#= Playing around with matching coefficients with Przybytek quintet potential,
with the aim of producing a scattering length

Description last updated 14/07/20=#

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using Revise
using OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays, SpecialFunctions
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
                    μ=0.5*4.002602u"u" # He₂ reduced mass [M]
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
    sol_input = x -> sol_unitless(ustrip(uconvert(u"bohr",x))) # add unit input
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
j(l,x)=sqrt(pi/(2*x))*besselj(l+1/2, x) # spherical bessel function jₗ(x)
n(l,x)=sqrt(pi/(2*x))*bessely(l+1/2, x) # spherical bessel function nₗ(x)

"""
Tries to find large-R limit of Aₗ, Bₗ, matching to jₗ(kR) and nₗ(kR)
Inputs: (u,u')(R), k, l, stapt, max_endpt
Outputs: (Aₗ,Bₗ)
#########OLD##########
Returns -tan(δₗ(k))/k from u(R)=Aₗjₗ(kR)+Bₗnₗ(kR)
Inputs: k (a₀⁻¹), l::ℕ; μ, U:[L]->[E]=przybytek, start pt [L] = 1a₀
Outputs: -tan(δₗ)/k
#########OLD##########
"""
function AB_limit(sol, #(u,u')(R) ~ [L]->(1,[L]⁻¹)
                  k, # [L]⁻¹
                  l::Int, # ang. momentum
                  stapt, # [L]
                  max_endpt, # [L]
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
    #= Old code that tried to do everything in one function
    ϵ = (1u"ħ^2")*k^2/(2*μ) #E=ħ²k²/2m
    endpt = noλ*(2*pi/k) # solve out to multiple of wavelength
    if l==0 #𝑠-wave case
        V=U #V₀(R) simply equals BO potential
        sol=rhs_solver(ϵ,stapt=stapt,endpt=endpt,U=V,μ=μ) # [u,u'](R)
        M = R -> [sin(k*R)   cos(k*R); # sin(kR), cos(kR) & derivs
                  k*cos(k*R) -k*sin(k*R)]
    else
        V(R)=auconvert(U(R)+l*(l+1)u"ħ^2"/(2*μ*R^2)) # add centrifugal potential
        sol=rhs_solver(ϵ,stapt=stapt,endpt=endpt,U=V,μ=μ) # [u,u'](R)
        M = R -> [j(l,k*R)                    n(l,k*R); #jₗ(kR), nₗ(kR) & derivs
                  (l/R)*j(l,k*R)-k*j(l+1,k*R) (l/R)*n(l,k*R)-k*n(l+1,k*R)]
    end
    AB = ustrip.(M(endpt))\ustrip.(sol(endpt)) # matching, stripping to do maths
    tanδₗ = AB[2]/AB[1] # tan(δₗ)=Bₗ/Aₗ
    return -tanδₗ/k
    =#
end


lhs=1.0u"bohr"
rhs=500.0u"bohr"
k=1e-6u"bohr^-1"
l=0
wavefn=rhs_solver(k,l,stapt=lhs,endpt=rhs)
ABfn=AB_limit(wavefn,k,l,lhs,rhs)
Rs=LinRange(lhs,rhs,1000)
ABs=ABfn.(Rs)
As=getindex.(ABs,1)
Bs=getindex.(ABs,2)
Aplt=plot(austrip.(Rs),As,label="A")
Bplt=plot(austrip.(Rs),Bs,label="B")
plot(Aplt,Bplt,layout=(2,1),title="l=$l, k=$k")
