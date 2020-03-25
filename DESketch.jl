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
grid,res,err=nDEulerSolver(expDensity,expLeft,expRight,expMatrix,expInitial,expAnalytic)
print("Error = $err\n")
plot(grid,res)
plot!(grid,expAnalytic.(grid))
end

# Conservation of population example
consMatrix=zeros(Float64,(2,2))
consMatrix[1,2],consMatrix[2,1]=1,-1
consInitial=[2,-1]
consAnalytic=(x->-sin(x)+2*cos(x))
consDensity,consLeft,consRight=100,0,2*π

function consDecayExample()
grid,res,err=nDEulerSolver(consDensity,consLeft,consRight,consMatrix,consInitial,consAnalytic)
print("Error = $err\n")
plot(grid,res)
plot!(grid,consAnalytic.(grid))
end

function iterator(key::string,lEnd,rEnd,lΔ,uΔ)
    
