#Demonstrated my Euler's method solver and the DifferentialEquations.jl package

push!(LOAD_PATH, pwd())
using nDES: nDEulerSolver
using DifferentialEquations
using Plots


# Exponential decay example
expMatrix=zeros(Float64,(1,1))
expMatrix[1,1]=-5 #exp. decay factor of 5 as given
expInitial=[1]
expAnalytic=(x->exp(-5-5*x))
expDensity,expLeft,expRight=100,-1,10

function expDecayExample()
    grid,res,err=nDEulerSolver('d',expDensity,expLeft,expRight,expMatrix,expInitial,expAnalytic)
    print("Error = $err\n")
    plot(grid,res)
    plot!(grid,expAnalytic.(grid))
end

# Conservation example
consMatrix=zeros(Float64,(4,4))
consMatrix[1,3],consMatrix[2,3]=1,-1
consMatrix[3,3],consMatrix[4,3]=-1,1
consInitial=[0,0,1,1]
consAnalytic=(x->1-exp(-x))
consDensity,consLeft,consRight=Int(1e6),0,10

function consExample()
    grid,res,err=nDEulerSolver('d',consDensity,consLeft,consRight,consMatrix,consInitial,consAnalytic)
    print("Error = $err\n")
    plot(grid,res)
    plot!(grid,consAnalytic.(grid))
end

# Wave example
waveMatrix=zeros(Float64,(2,2))
waveMatrix[1,2],waveMatrix[2,1]=1,-1
waveInitial=[2,-1]
waveAnalytic=(x->-sin(x)+2*cos(x))
waveDensity,waveLeft,waveRight=Int(1e6),0,2*π

function waveExample()
    grid,res,err=nDEulerSolver('d',waveDensity,waveLeft,waveRight,waveMatrix,waveInitial,waveAnalytic)
    print("Error = $err\n")
    plot(grid,res)
    plot!(grid,waveAnalytic.(grid))
end

function iterator(key::String,lΔ,uΔ,nΔ)
    ΔGrid=range(lΔ,uΔ,length=Int(nΔ))
    @assert key in ["exp","wave","cons"] "Invalid key: $key"
    errors=[]
    if key == "exp" #Exponential example
        for Δ in ΔGrid
            grid,res,err=nDEulerSolver('s',Δ,expLeft,expRight,expMatrix,expInitial,expAnalytic)
            push!(errors,err)
        end
    elseif key == "wave" #Wave equation example
        for Δ in ΔGrid
            grid,res,err=nDEulerSolver('s',Δ,waveLeft,waveRight,waveMatrix,waveInitial,waveAnalytic)
            push!(errors,err)
        end
    else #conservation example
        for Δ in ΔGrid
            grid,res,err=nDEulerSolver('s',Δ,consLeft,consRight,consMatrix,consInitial,consAnalytic)
            push!(errors,err)
        end
    end
    plot(ΔGrid,errors,xaxis=:log,yaxis=:log,
    xlabel="Step size", ylabel="Average error",
    title="$key n-dim Euler's Method solver, n=$nΔ")
end

################################################################################
#Now: repeat previous plots with the DifferentialEquations.jl package
function expDecay!(du,u,p,t)
    du[1]=-5u[1]
end
exp0 = [1.0]
expSpan=(-1.,10.)
expProb=ODEProblem(expDecay!,exp0,expSpan)
expSol=solve(expProb)
expPlot=plot(expSol,title="Exponential Decay solved using DE.jl", fmt = :pdf)
#savefig(expPlot, "expPlot_DE-package.pdf")

function transfer!(du,u,p,t)
    #u=[f,g]
    du[1] = -u[1] #f'=-f
    du[2] = u[1]   #g'= f
end
tr0 = [1.0;0.0] #f(0)=1, g(0)=0
trSpan = (0.,10.)
trProb=ODEProblem(transfer!,tr0,trSpan)
trSol=solve(trProb)
trPlot=plot(trSol,title="Transfer problem solved using DE.jl", fmt = :pdf)
#savefig(trPlot, "trPlot_DE-package.pdf")

function wave!(du,u,p,t)
    #u=[f,g]
    du[1] = u[2]
    du[2] = -u[1]
end
wave0 = [2.,-1.]
waveSpan = (0,2π)
waveProb=ODEProblem(wave!,wave0,waveSpan)
waveSol=solve(waveProb)
wavePlot=plot(waveSol,title="Wave equatoin solved using De.jl", fmt = :pdf)
#savefig(wavePlot, "wavePlot_DE-package.pdf")
