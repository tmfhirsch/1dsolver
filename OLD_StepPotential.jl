#= This code is my first attempt at solving the TISE numerically, for the simple
rectangular step potential V(x)={3,  1<x<2
                                {0,  else
Updated header 12-6-2020
=#

cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")
using DifferentialEquations
using Plots

tFine = LinRange(-10,15,1001) # spatial range

function Vstep(x) # step potential
    if 1<x<2
        return 1 # V₀ = 1 so plot is over ϵ/V₀
    else
        return 0
    end
end
function Vsink(x)
    if 1<x<2
        return -1 # V₀ = 1 so plot is over ϵ/V₀
    else
        return 0
    end
end
function Vpara(x) # parabolic potential
    if 1<x<2
        return -4*(x-3/2)^2+1
    else
         return 0
    end
end

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

function transmitted(en)
    m = 1.0 # mass
    ħ = 1.0 # reduced Planck's constant
    ϵ = en # energy
    k = sqrt(2*m*ϵ)/ħ # wavenumber

    # Schroedinger equation
    V(x)=DoubleWell(x)
    function TISE!(du,u,p,x)
        du[1] = u[2]  # ψ'(x)≡ψ'(x)
        du[2] = 2*m*(V(x)-ϵ)*u[1]/(ħ^2)       # (ψ'(x))'=2m(V(x)-E)ψ(x)/ħ^2
    end

    # Initial conditions
    ψ₀=exp(k*im) # Initial wavefunction; A=0, B=1
    ψD₀=-k*im*exp(k*im) # Initial wavefunction gradient
    u₀=[ψ₀, ψD₀] # Initial state vector
    xspan = (-10.0,15.0) #Solve over range

    prob=ODEProblem(TISE!,u₀,xspan)
    sol=solve(prob)

    # Solve right and left wave
    ABMatrix(t) = [exp(im*k*t)      exp(-im*k*t);
                   im*k*exp(im*k*t) -im*k*exp(-im*k*t)]
    AB(t)=ABMatrix(t)\sol(t)
    ABs=AB.(tFine)
    As=[i[1] for i in ABs]
    Bs=[i[2] for i in ABs]

    # Wavefunction plot
    #plot(plot(tFine, abs2.(getindex.(sol.(tFine),1))), plot(sol.t, V.(sol.t)), layout=(2,1))
    # Transmission coefficients plot
    #plot(plot(tFine, abs2.(As)), plot(tFine, abs2.(Bs)), layout=(2,1))
    #plot(tFine, [abs2.(As) abs2.(Bs)])

    T=abs2(Bs[1])/abs2(Bs[end])
    R=abs2(As[end])/abs2(Bs[end])

    return T # return transmission coefficient
end

ϵRange = LinRange(0.01,20,1000)
Ts = transmitted.(ϵRange)
plot(ϵRange,Ts,
     xlabel="Energy", ylabel="Transmission",
     legend=false,
     width=2,
     color=:royalblue)
#savefig("Double_Barrier_Transmission.svg")

plot(LinRange(2.5,10,100),DoubleWell.(LinRange(2.5,10,100)),
     xlabel="x",
     ylabel="V(x)",
     width=2,
     legend=false,
     color=:royalblue)
#savefig("Double_Barrier_Potential.svg")
#plot(tFine, Vstep.(tFine))
#plot!(tFine, Vsink.(tFine))
#plot!(tFine, Vpara.(tFine))
