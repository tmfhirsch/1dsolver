#=
Improved version of the Step Potential Code.
Input: Initial conditions, in form of [a,b] ~ [ψ,ψ']
Optional inputs: Potential (defaults to asymmetric double step potential);
    start and end points (defaults to 0,10);
    m, ħ, ϵ (=Energy) (defaults to 1,1,1);
    functions of (k,x) to match coefficients to (defaults to exp(±ikx)).
Output: named tuple of matching coefficients (cf1, cf2) and wavefunct (wav)
Output: [(A,B),(A,B)] such that ψ=A{fn1} + B{fm2}, at each end of xrange.
Description last updated 9/07/20
=#
cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")
using Plots, OrdinaryDiffEq, LinearAlgebra
using UnPack
using Revise

#Asymmetric potential (after code by Danny)
function asym_double_well(r; start=5, width=1, sep=1.2, h1=6, h2=10)::Float64
    for (len,val) in [(start,0),
                      (width, h1),
                      (sep, 0),
                      (width, h2)]
        r < len && return val
        r -= len
    end
    return 0
end

#=TISE solver, does the differential equations work
Inputs: energy; potential, other constants, limits
Output: matrix relating (ψ,ψ') on RHS to ICs=#
function schrodinger_solver(V, #potential
                            ϵ, #energy
                            m, #mass
                            ħ, #planck's constant
                            limits #limits
                            )
    k=sqrt(2*m*ϵ)/ħ
    function TISE!(du,u,p,x) # TISE gives (ψ,ψ')'=f((ψ,ψ'))
        du[1]=u[2] # ψ'(x)≡ψ'(x)
        du[2]=2*m*(V(x)-ϵ)*u[1]/(ħ^2) # (ψ'(x))'=2m(V(x)-E)ψ(x)/ħ^2
    end

    IC1, IC2 = [1.,0.], [0.,1.] # Basis of initial conditions
    prob1, prob2 = ODEProblem(TISE!,IC1,limits), ODEProblem(TISE!,IC2,limits)
    sol1, sol2 = solve(prob1,Tsit5()), solve(prob2,Tsit5())

    return [sol1(limits[2]) sol2(limits[2])]
end

#=Finds ICs that produce desired (ψ,ψ') on RHS or LHS.
Inputs: (ψ,ψ')', potential, energy, m, ħ, limits
Output: ICs in (ψ,ψ') basis that produce desired (ψ,ψ')'
=#
function ICfinder(rvals,
                  V,
                  ϵ, m, ħ,
                  limits
                  )
    k=sqrt(2*m*ϵ)/ħ
    M = schrodinger_solver(V, ϵ, m, ħ, limits) #Matrix mapping IC to RHS
    return M\rvals #linear comb of ICs that produce desired vals
end

#=Finds Transmission coefficient (not the norm sq. version)
Inputs: Energy, direction ∈ {"L2R", "R2L"};
basis functions (fn1(x) ->, fn2(x) <-), potential, energy, m, ħ, limits
Output: [Transmission coefficient, Reflection coefficient]=#
function transmission(ϵ, rlflag;
                      fn1 = "plane wave",
                      fn2 = "plane wave",
                      fD1 = "plane wave",
                      fD2 = "plane wave",
                      V = asym_double_well,
                      m = 1., ħ = 1.,
                      limits = (-10.,15.)
                      )
    @assert rlflag in ["L2R", "R2L"] "Please choose LHS or RHS as a string"

    k=sqrt(2*m*ϵ)/ħ

    if fn1 == "plane wave"
        func1(t) = exp(im*k*t) # defaults to exp(+i k x)
    else
        func1 = fn1
    end
    if fn2 == "plane wave"
        func2(t) = exp(-im*k*t) # defaults to exp(-i k x)
    else
        func2 = fn2
    end
    if fD1 == "plane wave"
        funD1(t) = im*k*exp(im*k*t) # defaults to +i k exp (+i k x)
    else
        funD1 = fD1
    end
    if fD2 == "plane wave"
        funD2(t) = -im*k*exp(-im*k*t) # defaults to -i k exp (-i k x)
    else
        funD2 = fD2
    end
    func_M(A,B,x) = [A*func1(x) B*func2(x);
                     A*funD1(x) B*funD2(x)]
    func_V(A,B,x) = [A*func1(x)+B*func2(x),
                     A*funD1(x)+B*funD2(x)]

    if rlflag == "L2R" # B'=0 case
        rAB = (1.,0.)
        rvals = func_V(rAB[1],rAB[2],limits[2]) # ψ = 1*f₁ + 0*f₂ on RHS
        lvals = ICfinder(rvals, V, ϵ, m, ħ, limits) # (ψ,ψ') on left
        lAB = func_M(1,1,limits[1])\lvals # A,B giving left vals
        T = rAB[1]/lAB[1] # transmission coefficient = A'/A
        R = lAB[2]/lAB[1] # reflection coefficient = B/A
    else # R2L transmission, A = 0 case
        lAB = (0., 1.)
        lvals = func_V(lAB[1],lAB[2],limits[1])
        M = schrodinger_solver(V, ϵ, m, ħ, limits) #TISE matrix
        rvals = M*lvals
        rAB = func_M(1,1,limits[2])\rvals # A,B giving right vals
        T = lAB[2]/rAB[2] # transmission coefficient = B/B'
        R = rAB[1]/rAB[2] # reflection coefficient = A'/B'
    end
    return T,R
end

ϵrange = LinRange(1,50,100)
Ts = [i[1] for i in transmission.(ϵrange, "L2R")]
Rs = [i[2] for i in transmission.(ϵrange, "L2R")]
plot(ϵrange, abs2.(Ts))
plot!(ϵrange, abs2.(Rs))
