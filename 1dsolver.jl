using Plots

#PARAMS
const gridSize=Int64(1e3) #number of data points
const lEnd=Int64(-1) #location of left end pt
const rEnd=Int64(1) #location of right end pt
const ψ0=1 #Boundary Condition

#derivative function
ψD(ψprev)=-5 * ψprev #A set here
soln=(x->exp(-5-5*x)).(collect(copy(xGrid)))

xGrid=range(lEnd,rEnd,length=gridSize)
Δx=(rEnd-lEnd)/gridSize
ψ=collect(copy(xGrid))
ψ[1]=ψ0
for i in 2:size(ψ)[1]
    ψ[i]=ψ[i-1] + Δx*ψD(ψ[i-1])
end


plot(xGrid,soln)
plot!(xGrid,ψ)
