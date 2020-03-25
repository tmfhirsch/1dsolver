module nDES

using Plots, LinearAlgebra

function nDEulerSolver(density::Float64,
                       lEnd::Float64, rEnd::Float64,
                       dMatrix, #the matrix of 1st order ODEs
                       initialValues::Array{Float64,1},
                       soln = false #analytic solution (a function)
                       )
@assert rEnd>lEnd "Right domain end not greater than left"
Δx = (rEnd-lEnd)/density
xGrid = range(lEnd, rEnd, step=Δx)
gridLength=size(xGrid)[1]
@assert size(dMatrix)[1]==size(dMatrix)[2] "DE Matrix is not square"
n = size(dMatrix)[1] #order of our ODE
f = initialValues #vector of 0 to (n-1) defivs. 1st index is 0th deriv
results=Array{Float64,1}(undef,gridLength) #store approximates here
results[1]=f[1] #initial value given
for i=2:gridLength #Propagate along the grid
  f += Δx * dMatrix * f #Euler's method
  results[i]=f[1] #Store result
end #for loop
end #function

end #module
