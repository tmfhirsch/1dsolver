push!(LOAD_PATH, pwd())
using nDES: nDEulerSolver
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
