 #= Playing around with matching coefficients with Przybytek quintet potential,
with the aim of producing a scattering length

Description last updated 14/07/20=#

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using Revise
using OrdinaryDiffEq, Plots, LinearAlgebra, StaticArrays, SpecialFunctions, Statistics
using Unitful, UnitfulAtomic
using Potentials: Quintet, Singlet, Triplet


"""
TISE solver for IC of (u,u')=(0,1)
All inputs in units, or assumed to be in atomic units if not.
Inputs: k, l; lhs~[L], rhs~[R], U(R)=Quintet, μ=He₂ mass, reltol for DE solver
Output: (u,u')~(a₀,1) as fn of R~a₀"""
function rhs_solver(k, # wavenumber [L]⁻¹
                    l::Int; # angular momentum
                    stapt=1.0u"bohr", # location of IC [L]
                    endpt=100.0u"bohr", # RHS to be solved to [L]
                    pot=Quintet, # interatomic potential [L]->[E]
                    μ=0.5*4.002602u"u", # He₂ reduced mass [M]
                    reltol=1e-10, #relative tolerance for DE solver
                    maxiters=1e6 #max iterations for DE solver
                    )
    ϵ=auconvert(k^2*1u"ħ^2"/(2*μ)) # E=ħ²k²/2μ
    # add centrifugal potential
    V_L = R -> auconvert(pot(R)+l*(l+1)u"ħ^2"/(2*μ*R^2))
    # check units are dimensionally correct
    @assert dimension(2*μ*(V_L(stapt)-ϵ)/1u"ħ^2")==dimension(1u"m^-2") "TISE units mismatch"
    # convert and strip units before TISE
    V = x -> austrip(V_L((x)u"bohr")) # unitless -> unitless potential
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
    # Add units back
    sol_unitless=solve(prob,Tsit5(),reltol=reltol,maxiters=maxiters) # comes out of DE sans units
    sol_input = x -> sol_unitless(austrip(x)) # add unit input
    sol = x -> sol_input(x).*SA[1u"bohr",1] # add unit output
    return sol
end

#Testing rhs_solver
function test_rhs_solver()
    plot_k=1e-5u"bohr^-1"
    plot_l=0
    plot_stapt=3.0u"bohr"
    plot_endpt=1e6u"bohr"
    sol=rhs_solver(plot_k,plot_l,stapt=plot_stapt,endpt=plot_endpt)
    Rs=LinRange(plot_stapt,plot_endpt,1000)
    us=[i[1] for i in sol.(Rs)]
    plot(ustrip.(Rs), ustrip.(us), legend=false, title="k=$plot_k")
    #savefig("prz_wfn_ϵ-$plot_ϵ.png")
end

"""
Spherical bessel functions
"""
j(l,x)=sphericalbesselj(l,x)
n(l,x)=sphericalbessely(l,x)

"""Zero potential, for testing functions"""
zero_pot(R)=0.0u"hartree"

"""Hard sphere potential, for testing functions"""
hardsphere(R; rad=1u"bohr") = R > rad ? 0u"hartree" : 1e3u"hartree"

"""
Returns (Aₗ,Bₗ)(R) that match wavefunction to spherical bessel functions
Inputs: (u,u')(R) solution ~ ([L],1), k ~ [L⁻¹], l
Outputs: (Aₗ,Bₗ)(R)
"""
function matchAB(sol, #(u,u')(R) ~ [L]->([L],1)
                 k, # [L]⁻¹
                 l::Int # ang. momentum
                 )
    @assert l>=0 "Not a valid angular momentum"
    # define matrix of matching functions M(R) ~ [L] -> [1 1; [L]⁻¹ [L]⁻¹]
    M = R -> [R*j(l,k*R)                    R*n(l,k*R); #Rjₗ(kR), Rnₗ(kR) & derivs
              (l+1)*j(l,k*R)-k*R*j(l+1,k*R) (l+1)*n(l,k*R)-k*R*n(l+1,k*R)]
    AB = R -> ustrip.(M(R))\ustrip.(sol(R)) # AB matched at R ~ [L]->[1,1]
    return AB
end

# Testing that A, B match the analytic soln for hard sphere scattering
# 24/07/20 agreement shown, see Slack message to danny from this date
function AB_analytic_check(rad=1.0u"bohr", evalpt=100.0u"bohr", l=1::Int)
    Avals, Bvals = [], [] # store for AB(evalpt)
    Asols, Bsols = [], [] # store for analytic soln for AB(endpt)
    kRange=LinRange(1e-2,1e-12,50)u"bohr^-1"
    pot = R -> hardsphere(R, rad=rad)
    for k in kRange
        sol=rhs_solver(k, l, pot=pot, stapt=rad, endpt=evalpt)
        AB = matchAB(sol, k, l)
        ABeval = AB(evalpt)
        append!(Avals, ABeval[1]) # store A evaluation
        append!(Bvals, ABeval[2]) # store B evaluation
        analA = -rad*k*n(l,rad*k) # analytic A soln
        analB = +rad*k*j(l,rad*k) # analytic B soln
        append!(Asols, analA) # store A analytic solution
        append!(Bsols, analB) # store B analytic solution
    end

    Aplt = plot(austrip.(kRange), Avals, labels=["calculated"]); plot!(austrip.(kRange), austrip.(Asols),
    xlabel="k (a₀⁻¹)", ylabel="A", labels=["analytic"], xscale=:log10)
    Bplt = plot(austrip.(kRange), Bvals, labels=["calculated"]); plot!(austrip.(kRange), austrip.(Bsols),
    xlabel="k (a₀⁻¹)", ylabel="B", labels=["analytic"], xscale=:log10)
    return plot(Aplt, Bplt, layout=(2,1), title=["Hard sphere; A at 100a₀, calc vs. analytic" "Hard sphere; B at 100a₀, calc vs. analytic"])
end

 # Testing matchAB
#=lhs=1.0u"bohr"
rhs=1e2u"bohr"
k=1e-12u"bohr^-1"
l=1
no_pts=1000
wavefn=rhs_solver(k,l,pot=zero_pot,stapt=lhs,endpt=rhs,reltol=1e-10)
ABfn=matchAB(wavefn,k,l)
Rs=LinRange(lhs,rhs,no_pts)
ABs=ABfn.(Rs)
As=getindex.(ABs,1)
Bs=getindex.(ABs,2)
plot(austrip.(Rs),austrip.(getindex.(wavefn.(Rs),1)))
Aplt=plot(austrip.(Rs),As,label="A")
Bplt=plot(austrip.(Rs),Bs,label="B")
plot(Aplt,Bplt,layout=(2,1),title="l=$l, k=$k")
plot(austrip.(Rs),Bs./As,xlabel="R (a₀)",ylabel="B(R)/A(R)",title="l=$l,k=$k",legend=false)
tan_δ₀ = mean((Bs./As)[Int(floor(no_pts/2)):end])
σ=std((Bs./As)[Int(floor(no_pts/2)):end])=#


"""
Finds convergence of B(R)/A(R)=tan(δ)
Input: AB(R) ~ [L]->1; no. grid pts, median bubble, convergence tolerance,
warning 10^index, break10^index
Output: lim(r→∞) B(R)/A(R) within tolerance, or an error if 100*stop reached
"""
function BoA_lim(AB; no_pts=1000::Int, bub=1.0, tol=1e0, warn=50::Int, stop=100::Int)
    BoA(R) = AB(R)[2]/AB(R)[1] #B/A -> tanδₗ
    i=1 # Index for R∈(100i,100(i+1))
    while true
        Rs=LinRange(100i,100*(i+1),no_pts)u"bohr" # linear grid of R values
        BoAs=BoA.(Rs) # AB values on the grid
        med = median(BoAs)
        fil = x -> abs(x-med)<bub # filter within <bub> of the median
        BoAs=filter(fil, BoAs)
        if std(BoAs) < tol
            return mean(BoAs)
            break
        end
        if i == warn
            @warn "B(R)/A(R) hasn't converged by $(warn)00 a₀"
        elseif i == stop
            error("B(R)/A(R) did not converge by $(stop)00 a₀")
            break
        end
        i += 1
    end
end

# Testing BoA_lim
#=lhs=0.615124u"bohr"
rhs=1e6u"bohr"
k=1e-4u"bohr^-1"
l=1
wavefn=rhs_solver(k,l,pot=zero_pot,stapt=lhs,endpt=rhs)
ABfn=matchAB(wavefn,k,l)
Rs=LinRange(lhs,rhs,1000)
As, Bs = getindex.(ABfn.(Rs),1), getindex.(ABfn.(Rs),2)
plot(ustrip.(Rs), Bs./As,xlabel="R (a₀)",ylabel="B/A",legend=false,title="k=$k,l=$l")
str="lhs=$lhs, rhs=$rhs, k=$k, l=$l, tan(δₗ) = "*string(BoA_lim(ABfn,tol=1e-4))
@info str =#

"""
Finds scattering length of a given potential
Input: Potential ~ [L]->[E]; lhs~[L], rhs~[R], μ~[M], convergence tolerance,
warning index, max index
Output: a in nanometres, or an error if max index reached
"""
function scatlen(pot; # potential
                 lhs=3.0u"bohr", # start point for solving DE
                 rhs=1e4u"bohr", # end pt for solving DE
                 μ=0.5*4.002602u"u", # reduced mass
                 strt=5.0, # start i for k=10^-i
                 warn=10::Int,
                 stop=15::Int,
                 tol=1e-2u"nm" #tolerance in nm on the scattering length
                 )
    i=strt
    while true
        # wavenumbers to test for convergence
        k1 = 10^(-i)*1u"bohr^-1"
        k2 = 10^(-i-1)*1u"bohr^-1"
        k3 = 10^(-i-2)*1u"bohr^-1"
        sol1 = rhs_solver(k1, 0, pot=pot, stapt=lhs, endpt=rhs)
        sol2 = rhs_solver(k2, 0, pot=pot, stapt=lhs, endpt=rhs)
        sol3 = rhs_solver(k3, 0, pot=pot, stapt=lhs, endpt=rhs)
        AB1 = matchAB(sol1, k1, 0) # scattering length: l=0
        AB2 = matchAB(sol2, k2, 0)
        AB3 = matchAB(sol3, k3, 0)
        tanδ1 = BoA_lim(AB1)
        tanδ2 = BoA_lim(AB2)
        tanδ3 = BoA_lim(AB3)
        ai = uconvert.(u"nm", -[tanδ1/k1, tanδ2/k2, tanδ3/k3]) #scat lgths (nm)
        if maximum(ai)-minimum(ai) < tol # variation in scat lengths withinin tol
            return ai[end] # taking last one/smallest k value
            break
        end
        if i == warn
            @warn "tan(δ(k))/k hasn't converged by k=$k1"
        elseif i == stop
            error("tan(δ(k))/k didn't converge by k=$k1")
        end
        i += 1
    end
end

# Testing scatlen
#@info scatlen(Quintet)


"""
Finds partial cross section σₗ
Input: k~[L]⁻¹, l; potential~[L]->[E], lhs~[L], rhs~[R], μ~[M]
Outputs: σₗ~[L]^2 (cm^2 by default)
"""
function partialσ(k,
                  l::Int;
                  pot=Quintet,
                  lhs=1.0u"bohr", #lhs for DE solver
                  rhs=1e4u"bohr", #rhs for DE solver
                  μ=0.5*4.002602u"u" # reduced mass
                  )
    sol = rhs_solver(k,l,pot=pot,stapt=lhs,endpt=rhs, μ=μ)
    AB = matchAB(sol,k,l) # (A,B)(R)
    BoA = BoA_lim(AB) # lim(R→∞) B(R)/A(R) ≈ tan(δ₀)
    δₗ = atan(BoA) # partial phase shift
    σₗ = 4*pi*(2*l+1)*sin(δₗ)^2/k^2 # partial cross section
    return uconvert(u"cm^2", σₗ)
end

"""Testing hard sphere partial cross section vs energy
24/7/20: produces a good plot"""
function σ_energy_check(rad=1.0u"bohr")
    ls=[0,1,2,3,4]
    kRange = LinRange(1e1,1e-12,50)u"bohr^-1"
    pot = R -> hardsphere(R, rad=rad)
    σls = []
    labs= []
    for l in ls
        push!(labs, "l=$l")
        σs = Float64[]
        for k in kRange
            push!(σs, ustrip.(partialσ(k,l,pot=pot,lhs=rad)))
        end
        push!(σls,σs)
    end
    σls=reshape(hcat(σls...), (length(kRange),length(ls)))
    labs=reshape(convert(Array{String,1},labs), 1, length(labs))
    σt=zeros(length(kRange),1)
    for i in 1:length(kRange)
        σt[i]=sum(σls[i,:])
    end
    labs=hcat(labs, ["total"])
    return plot(austrip.(kRange), hcat(σls,σt), xlabel="k (a₀⁻¹)", ylabel="σₗ (cm²)",
    title="Hard sphere; partial cross sections vs. wavenumber",labels=labs)
end

"""Testing hard sphere cross sections; l vs k array and cumsum array
25/07/20: produces the correct convergence with l↑, more l needed for k↑"""
function σ_matrix_check(rad=1.0u"bohr")
    kRange=exp10.(1:-1:-12)u"bohr^-1"
    lRange=0:1:20
    pot = R -> hardsphere(R, rad=rad)
    nok, nol = length(kRange), length(lRange) #ls form rows, ks form columns
    lk_matrix=fill(0.0u"cm^2",(nol,nok)) # initialise matrix
    for i=1:nol, j=1:nok
        lk_matrix[i,j]=partialσ(kRange[j],lRange[i],pot=pot)
    end
    cs_matrix=cumsum(lk_matrix,dims=1)
    return lk_matrix, cs_matrix
end
