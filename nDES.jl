#n Dimensional Euler method differential equation solver
module nDES

using LinearAlgebra

function nDEulerSolver(flag, #flag density or step size
                       dOrS, #value of density or step size
                       lEnd,
                       rEnd,
                       dMatrix, #the matrix of 1st order ODEs
                       initialValues,
                       soln = false #analytic solution (a function)
                       )
    @assert rEnd>lEnd "Right domain end not greater than left"
    @assert flag in ['d','s'] "Invalid density/step size flag"
    if flag=='d' #density given
        Δx = (rEnd-lEnd)/dOrS
    else
        Δx = dOrS #step size given
    end
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
    ϵ=false #Returns false for error if no analytic comparison
    if soln != false #if analytic function argument present
        aGrid=soln.(xGrid)
        ϵ=sum([abs(results[i]-aGrid[i]) for i=1:gridLength])/gridLength
    end #if
    return xGrid, results, ϵ
end #function

end #module
