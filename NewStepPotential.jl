#=
Improved version of the Step Potential Code.
Input: Initial conditions, in form of [a,b] ~ [ψ,ψ']
Optional inputs: Potential (defaults to asymmetric double step potential);
    start and end points (defaults to 0,10);
    location of ICs (defaults to 0);
    m, ħ, ϵ (=Energy) (defaults to 1,1,1);
    functions of (k,x) to match coefficients to (defaults to exp(±ikx)).
Output: named tuple of matching coefficients (cf1, cf2) and wavefunct (wav)
Output: [(A,B),(A,B)] such that ψ=A{fn1} + B{fm2}, at each end of xrange.
Description last updated 8/07/20
=#
cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")
using Plots, OrdinaryDiffEq
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

#Main function
function coeffs(IC;
                V = asym_double_well,
                ϵ=1.0, m=1.0, ħ=1.0,
                limits=(-10.,20.0),
                IC_loc=0.0,
                fn1="plane wave",
                fn2="plane wave",
                fD1="plane wave",
                fD2="plane wave",
                )
    k=sqrt(2*m*ϵ/ħ)
    if fn1 == "plane wave"
        func1(k,t) = exp(im*k*t) # defaults to exp(+i k x)
    else
        func1 = fn1
    end
    if fn2 == "plane wave"
        func2(k,t) = exp(-im*k*t) # defaults to exp(-i k x)
    else
        func2 = fn2
    end
    if fD1 == "plane wave"
        funD1(k,t) = im*k*exp(im*k*t) # defaults to +i k exp (+i k x)
    else
        funD1 = fD1
    end
    if fD2 == "plane wave"
        funD2(k,t) = -im*k*exp(-im*k*t) # defaults to -i k exp (-i k x)
    else
        funD2 = fD2
    end

    function TISE!(du,u,p,x) # TISE gives (ψ,ψ')'=f((ψ,ψ'))
        du[1]=u[2] # ψ'(x)≡ψ'(x)
        du[2]=2*m*(V(x)-ϵ)*u[1]/(ħ^2) # (ψ'(x))'=2m(V(x)-E)ψ(x)/ħ^2
    end

    prob=ODEProblem(TISE!,IC,limits)
    sol=solve(prob,Tsit5()) # using Tsit5() algorithm

    ABMatrix(t) = [func1(k,t) func2(k,t)
                   funD1(k,t) funD2(k,t)]
    AB(t)=ABMatrix(t)\sol(t)

    xspan=LinRange(limits[1],limits[2],1000)

    return (xspan = xspan,
            wav = getindex.(sol.(xspan),1),
            cf1 = AB(limits[1]),
            cf2 = AB(limits[2])
            )
end

@unpack xspan, wav = coeffs([0,1.])
plot(xspan,abs2.(wav))
