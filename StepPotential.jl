#= This code is my first attempt at solving the TISE numerically, for the simple
rectangular step potential V(x)={3,  1<x<2
                                {0,  else
Updated header 12-6-2020
=#

using DifferentialEquations
using Plots

m = 1.0 # mass
ħ = 1.0 # reduced Planck's constant
ϵ = 2 # energy

function V(x) # potential
    if 1<x<2
        return 3
    else
        return 0
    end
end

function TISE!(du,u,p,x)
    du[1] = u[2]  # ψ'(x)≡ψ'(x)
    du[2] = 2*m*(V(x)-ϵ)*u[1]/(ħ^2)       # (ψ'(x))'=2m(V(x)-E)ψ(x)/ħ^2
end

ψ₀=1 # Initial wavefunction
ψD₀=0.2 # Initial wavefunction gradient
u₀=[ψ₀, ψD₀] # Initial state vector
xspan = (0.0,4.0) #Solve over range (0,4)

prob=ODEProblem(TISE!,u₀,xspan)
sol=solve(prob)
ψsq=(x->abs(x)^2).([i[1] for i in sol.u]) #Probability density

plot(sol.t,ψsq)
