# Unit tests. Tested 9/9/20
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using Revise
using CrossSections: solver, K_matrix, F_matrix, σ_matrix
using StateStructures, Interactions
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots

# test function for solver - runs and plots first channel wavefunction
function test_solver(;lmax=0,B=0u"T")
    lhs, rhs = 3.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=SmS_lookup_generator(lmax)
    n=length(lookup)
    # construct ICs
    IC=[fill(0.0u"bohr",n,n); I]
    println("Starting to solve for wavefunctions, lmax=$lmax")
    sol=solver(lookup, IC, ϵ, lhs, rhs,B=B)
    #= plot code only works if save_everywhere=false is changed in solver()
    plt = let
        Rs=LinRange(lhs,rhs,1000)
        vals=getindex.(sol.(Rs),n,n) # first channel, first IC
        plot(austrip.(Rs),real.(austrip.(vals)))
    end
    plt=#
    println("Solved!")
end

# unit test for K_matrix. Should produce a scattering length of 7.54 nm
# in agreement with Przybytek
function test_K_matrix(;lmax=0, ϵ=1e-12u"hartree", μ=0.5*4.002602u"u",
    lhs=3.0u"bohr", rhs=1000u"bohr")
    println("Starting test_K_matrix")
    lookup=SmS_lookup_generator(lmax)
    n=length(lookup)
    # construct ICs
    IC=[fill(0.0u"bohr",n,n); I]
    # solver for wavefunctions
    println("Solving for wavefunctions")
    sol=solver(lookup,IC,ϵ,lhs,rhs,μ=μ)
    eval=sol(rhs)
    # solve 𝐤 vector for K matrix solver
    # ***Warning: the following code assumes all channels are open***
    println("Producing k vector")
    𝐤=fill(0.0u"bohr^-1",n)
    for i in 1:n
        γ = lookup[i] # channel
        R∞ = Inf*1u"bohr"
        V∞ = H_el(γ,γ,R∞) + H_sd_coeffs(γ,γ)*H_sd_radial(R∞) + H_rot(γ,γ,R∞,μ)
        𝐤[i] = sqrt(2*μ*(ϵ-V∞))/1u"ħ"
    end
    # 𝐥 vector for K matrix solver
    println("Producing 𝐥 vector")
    𝐥=fill(0,n)
    for i in 1:n
        𝐥[i]=lookup[i].l
    end
    println("Passing to K_matrix function")
    𝐊=K_matrix(rhs, eval, 𝐤, 𝐥)
    𝐒=(I+im*𝐊)*inv(I-im*𝐊)
    return uconvert(u"nm", sqrt(pi*abs(1-𝐒[n,n])^2/𝐤[n]^2/(4*pi)))
end

# unit test for F_matrix. Should produce a matrix of wavefunction solutions
function test_F_matrix(;lmax=0, ϵ=1e-12u"hartree", μ=0.5*4.002602u"u",
    lhs=3.0u"bohr", mid=100.0u"bohr", rhs=1000.0u"bohr")
    println("Running test_F_matrix")
    println("Initialising AL, BR, isOpen")
    # arbitrary sample AL, AR
    lookup=SmS_lookup_generator(lmax)
    N=length(lookup)
    # construct isOpen
    isOpen=fill(true,N)
    for i in 1:N
        γ = lookup[i] # channel
        R∞ = Inf*1u"bohr"
        V∞ = H_el(γ,γ,R∞) + H_sd_coeffs(γ,γ)*H_sd_radial(R∞) + H_rot(γ,γ,R∞,μ)
        ksq = 2*μ*(ϵ-V∞)/1u"ħ^2"
        isOpen[i] = austrip(ksq) >= 0 ? true : false # k² > 0 for open channels
    end
    # construct BCs
    AL=[fill(0.0u"bohr",N,N); I]
    BR = let
        Nₒ=count(isOpen)
        BFL = [fill(0.0u"bohr",N,N); I]
        BFR = [Matrix(Diagonal(ones(N))[:,isOpen]u"bohr"); zeros(N,Nₒ)]
        [BFL BFR]
    end
    # solve for solutions
    println("Solving for AR and BL")
    AR = solver(lookup, AL, ϵ, lhs, mid)(mid)
    BL = solver(lookup, BR, ϵ, rhs, mid)(mid)
    # see if F_matrix runs
    println("Passing to F_matrix")
    𝐅=F_matrix(AL,AR,BL,BR,isOpen)
    println("Finished test_F_matrix")
    𝐅
end

# combined tests for F and K functions. Should produce a Quintet scattering
# length of 7.54nm in agreement with Przybytek
function test_K_and_F(;lmax=0, ϵ=1e-12u"hartree", μ=0.5*4.002602u"u",
    lhs=3.0u"bohr", mid=100.0u"bohr", rhs=1000.0u"bohr")
    # calculate F
    println("Calculating F")
    𝐅=test_F_matrix(lmax=lmax,ϵ=ϵ,μ=μ,lhs=lhs,mid=mid,rhs=rhs)
    # calculate 𝐤 vector for K_matrix()
    lookup=SmS_lookup_generator(lmax)
    println("Calculating k vector")
    n=length(lookup)
    𝐤=fill(0.0u"bohr^-1",n)
    for i in 1:n
        γ = lookup[i] # channel
        R∞ = Inf*1u"bohr"
        V∞ = H_el(γ,γ,R∞) + H_sd_coeffs(γ,γ)*H_sd_radial(R∞) + H_rot(γ,γ,R∞,μ)
        𝐤[i] = sqrt(2*μ*(ϵ-V∞))/1u"ħ"
    end
    # 𝐥 vector for K matrix solver
    println("Constructing l vector")
    𝐥=fill(0,n)
    for i in 1:n
        𝐥[i]=lookup[i].l
    end
    println("Calculating K")
    𝐊=K_matrix(rhs,𝐅,𝐤,𝐥)
    𝐒=(I+im*𝐊)*inv(I-im*𝐊)
    return uconvert(u"nm", sqrt(pi*abs(1-𝐒[n,n])^2/𝐤[n]^2/(4*pi)))
end
