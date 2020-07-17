#= Playing around with matching coefficients with Przybytek quintet potential,
with the aim of producing a scattering length

Description last updated 14/07/20=#

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using Revise
using OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays, SpecialFunctions
using Unitful, UnitfulAtomic
using Przybytek: przybytek


"""TISE solver for IC of (u,u')=(0,1)
All inputs in units, or assumed to be in atomic units if not.
Inputs: energy, limits; V(R)=przybytek, m=He₂ mass
Output: (u,u') on RHS as fn of x (a₀)"""
function rhs_solver(ϵ; # energy [E]
                    stapt=20.0u"bohr", # location of IC [L]
                    endpt=500.0u"bohr", # RHS to be solved to [L]
                    U=przybytek, # interatomic potential [L->E]
                    μ=0.5*4.002602u"u" # He₂ reduced mass [M]
                    )
    # check units are dimensionally correct
    @assert dimension(2*μ*(U(stapt)-ϵ)/1u"ħ^2")==dimension(1u"m^-2") "TISE units mismatch"
    # convert and strip units before TISE
    ϵ = austrip(ϵ)
    stapt, endpt = austrip(stapt), austrip(endpt)
    V(x)=austrip(U((x)u"bohr"))
    μ = austrip(μ)
    ħ = austrip(1u"ħ")
    # TISE solver
    function TISE(u,p,x) # TISE gives (u,u')'=f((u,u'))
        du = u[2] # u'(x)≡u'(x)
        dd = 2*μ*(V(x)-ϵ)*u[1]/ħ^2 # (u'(x))'=2m(V(x)-E)u(x)/ħ^2 {𝐋⁻²}
        SVector{2}([du,dd]) # (u,u')'
    end
    IC = SVector{2}([0.0, 1.0]) # (u,u')=(0,1)
    prob=ODEProblem(TISE,IC,(stapt,endpt))
    sol_unitless=solve(prob,Tsit5(),reltol=1e-10)
    sol_input(x)=sol_unitless(ustrip(uconvert(u"bohr",x)))
    sol(x)=sol_input(x).*SA[1,1u"bohr^-1"]
    return sol
end
#=
#Test plot of rhs_solver
plot_ϵ=1e-10u"hartree"
plot_stapt=1.0u"bohr"
plot_endpt=20000.0u"bohr"
sol=rhs_solver(plot_ϵ,stapt=plot_stapt,endpt=plot_endpt)
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
Returns tan(δₗ) from u(R)=Aₗjₗ(kR)+Bₗnₗ(kR) using przybytek potential
Inputs: k (a₀⁻¹), l::ℕ
Outputs: tan(δₗ)
"""
function tan_shift(l::Int, #ang mom
                   k; #[L]⁻¹
                   μ=0.5*4.002602u"u", #[M]
                   U=przybytek, #[L]->[E]
                   stapt=1u"bohr" #[L]
                   )
    @assert l>=0 "Not a valid angular momentum"
    ϵ = (1u"ħ^2")*k^2/(2*μ) #E=ħ²k²/2m
    endpt = 2*pi/k #Assuming u(R) stabilises by one λ (phenomonological)
    if l==0 #𝑠-wave case
        V=U #V₀(R) simply equals BO potential
        sol=rhs_solver(ϵ,stapt=stapt,endpt=endpt,U=V,μ=μ) # [u,u'](R)
        M = R-> [sin(k*R)      cos(k*R); # sin(kR), cos(kR) & derivs
                 k*cos(k*R)    -k*sin(k*R)]
    else
        V(R)=auconvert(U(R)+l*(l+1)u"ħ^2"/(2*μ*R^2)) # add centrifugal potential
        sol=rhs_solver(ϵ,stapt=stapt,endpt=endpt,U=V,μ=μ) # [u,u'](R)
        M = R -> [j(l,k*R)                        n(l,k*R);
                  (l/R)*j(l,k*R)-k*j(l+1,k*R)     (l/R)*n(l,k*R)-k*n(l+1,k*R)]
        #M(R)=[j(l,k*R)                        n(l,k*R);
        #      (l/R)*j(l,k*R)-k*j(l+1,k*R)     (l/R)*n(l,k*R)-k*n(l+1,k*R)]
    end
    AB = ustrip.(M(endpt))\ustrip.(sol(endpt)) # matching, stripping to do maths
    return AB[2]/AB[1] # tan(δₗ)=Bₗ/Aₗ
end


"""
#TODO scattering length finder
"""
#Attempt to find scattering length
#=
as=[] # storage for -tan(δ₀)/k for different ks. Hopefully →a as k→0
for k in exp10.(LinRange(-5,-7,10))u"bohr^-1" #k ~ 10⁻¹⁰Eₕ to 10⁻¹² Eₕ
    push!(as, -tan_shift(k)/k)
end
as=uconvert.(u"nm", as)
=#
