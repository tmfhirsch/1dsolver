#= This code is my first attempt at solving the TISE numerically, for the simple
rectangular step potential V(x)={3,  1<x<2
                                {0,  else
Updated header 12-6-2020
=#

using DifferentialEquations
using Plots
using Statistics

tFine = LinRange(-10,15,1001) # spatial range

function transmitted(en)
    m = 1.0 # mass
    ħ = 1.0 # reduced Planck's constant
    ϵ = en # energy
    k = sqrt(2*m*ϵ)/ħ # wavenumber

    ψ₀=exp(20*im) # Initial wavefunction; A=0, B=1
    ψD₀=-2*im*exp(20*im) # Initial wavefunction gradient
    u₀=[ψ₀, ψD₀] # Initial state vector
    xspan = (-10.0,15.0) #Solve over range (0,4)
    #=
    function V(x) # step potential
        @assert flag=="step" "Wrong flag"
        if 1<x<2
            return 1 # V₀ = 1 so plot is over ϵ/V₀
        else
            return 0
        end
    end
    =#
    #=
    function V(x)
        @assert flag=="sink" "Wrong flag"
        if 1<x<2
            return -1 # V₀ = 1 so plot is over ϵ/V₀
        else
            return 0
        end
    end
    =#

    function V(x) # parabolic potential
        @assert flag=="para" "Wrong flag"
        if 1<x<2
            return -4*(x-3/2)^2+1
        else
            return 0
        end
    end


    function TISE!(du,u,p,x)
        du[1] = u[2]  # ψ'(x)≡ψ'(x)
        du[2] = 2*m*(V(x)-ϵ)*u[1]/(ħ^2)       # (ψ'(x))'=2m(V(x)-E)ψ(x)/ħ^2
    end

    prob=ODEProblem(TISE!,u₀,xspan)
    sol=solve(prob)

    ABMatrix(t) = [exp(im*k*t)      exp(-im*k*t);
                   im*k*exp(im*k*t) -im*k*exp(-im*k*t)]
    # [A,B] = InvABMatrix.[ψ,ψ']
    # A -> ; B <-
    #InvABMatrix(t) = 0.5*[exp(-im*k*t) -im*exp(-im*k*t)/k;
    #                    exp(+im*k*t) +im*exp(+im*k*t)/k]
    #AB(t)=InvABMatrix(t)*sol(t)
    AB(t)=ABMatrix(t)\sol(t)
    ABs=AB.(tFine)
    As=[i[1] for i in ABs]
    Bs=[i[2] for i in ABs]

    # Wavefunction plot
    #lot(plot(tFine, abs2.(getindex.(sol.(tFine),1))), plot(sol.t, V.(sol.t)), layout=(2,1))
    # Transmission coefficients plot
    #plot(plot(tFine, abs2.(As)), plot(tFine, abs2.(Bs)), layout=(2,1))
    #plot(tFine, [abs2.(As) abs2.(Bs)])

    T=mean(abs2.(Bs[1:50]))/mean(abs2.(Bs[end-50:end]))
    R=mean(abs2.(As[end-50:end]))/mean(abs2.(Bs[end-50:end]))

    return T # return transmission coefficient
end

flag="para" # type of potential

ϵRange = LinRange(0.01,2,1000)
Ts = transmitted.(ϵRange)
gui()
plot(ϵRange,Ts,
     title="Transmission coefficient, $flag potential",
     xlabel="ϵ/V₀", ylabel="T",
     fmt=:pdf)


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
plot(tFine, Vstep.(tFine))
plot!(tFine, Vsink.(tFine))
plot!(tFine, Vpara.(tFine))
