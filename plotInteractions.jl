# produce plots of the different Interactions in my model, for investigation
# and for thesis writing
# header last updates 23/09/20

using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using StateStructures, Interactions
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots, Plots.PlotMeasures
using BSON, Dates

const lmaxx=20 # more than enough states to be relevant to my calculations
# states up to lmax of 20 (more than enough)
lookup=SmS_lookup_generator(lmaxx)

"""Saves plots of interactions as .svg files in /PHYS4110/Plots/Interactions/"""
function save_plt(plt::Plots.Plot,title::String)
    prevdir=pwd()
    cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Plots")
    savefig(plt, title*".svg")
    cd(prevdir)
end

# Rotational Interaction
function plt_H_rot(;lmax=6,Rmin=3e0u"bohr",Rmax=1e6u"bohr",no_R_pts=10000,μ=0.5*4.002602u"u")
    firstl(i::Int)=lookup[findall(x->x.l==i,lookup)[1]] # grabs first state with correct l
    ls = [firstl(i) for i in 1:lmax] # prototypical l states
    lab=permutedims((x->"l=$(x.l)").(ls)) # labels for plot
    Rs=LinRange(Rmin,Rmax,no_R_pts)
    vals=zeros(length(Rs),length(ls))u"hartree" #initialse
    for i in 1:no_R_pts, j in 1:lmax
        vals[i,j]=H_rot(ls[j],ls[j],Rs[i],μ)
    end
    plt=plot(austrip.(Rs),austrip.(vals),xscale=:log10,yscale=:log10,
    xlabel="R (a₀)", ylabel="Hᵣₒₜ(l,R) (Eₕ)",label=lab,legend=:outertopright,
    bottom_margin=5mm,left_margin=5mm, grid=false,
    yticks=exp10.(-15:3:-3))
end

# Electronic Interaction plotting
function plt_H_el(;Rmin=3.0u"bohr",Rmax=20u"bohr",no_R_pts=10000)
    firstS(i::Int)=lookup[findall(x->x.S==i,lookup)[1]] # grabs first state with correct S
    Ss = [firstS(S) for S in 0:2] # prototypical S states
    lab=permutedims((x->"S=$(x.S)").(Ss))
    Rs=LinRange(Rmin,Rmax,no_R_pts)
    vals=zeros(length(Rs),length(Ss))u"hartree" #initialse
    for i in 1:no_R_pts, j in 1:length(0:2)
        vals[i,j]=H_el(Ss[j],Ss[j],Rs[i])
    end
    plt=plot(austrip.(Rs),austrip.(vals),
    xlabel="R (a₀)", ylabel="Hₑₗ(R) (Eₕ)",label=lab,legend=:outertopright)
end

# Electronic and rotational overlaid
function plt_H_el_and_rot(;lmax=6,Rmin=180.0u"bohr",Rmax=220.0u"bohr",no_R_pts=10000,μ=0.5*4.002602u"u")
    firstl(i::Int)=lookup[findall(x->x.l==i,lookup)[1]] # grabs first state with correct l
    ls = [firstl(i) for i in 0:lmax] # prototypical l states
    lab=permutedims((x->"l=$(x.l)").(ls)) # labels
    firstS(i::Int)=lookup[findall(x->x.S==i,lookup)[1]] # grabs first state with correct S
    Ss = [firstS(S) for S in 0:2] # prototypical S states
    Rs=LinRange(Rmin,Rmax,no_R_pts)
    # S=0 plot
    vals0=zeros(length(Rs),length(ls))u"hartree"
    for i in 1:no_R_pts, j in 1:(lmax+1)
        vals0[i,j]=H_el(Ss[1],Ss[1],Rs[i])+H_rot(ls[j],ls[j],Rs[i],μ)
    end
    lab0=permutedims((x->"l=$(x.l)").(ls))
    plt0=plot(austrip.(Rs),austrip.(vals0),label=lab,legend=:outertopright,
    title="S=0")
    # S=1 plot
    vals1=zeros(length(Rs),length(ls))u"hartree"
    for i in 1:no_R_pts, j in 1:(lmax+1)
        vals1[i,j]=H_el(Ss[2],Ss[2],Rs[i])+H_rot(ls[j],ls[j],Rs[i],μ)
    end
    lab1=permutedims((x->"l=$(x.l)").(ls))
    plt1=plot(austrip.(Rs),austrip.(vals1),label=lab,legend=:outertopright,
    title="S=1")
    # S=2 plot
    vals2=zeros(length(Rs),length(ls))u"hartree"
    for i in 1:no_R_pts, j in 1:(lmax+1)
        vals2[i,j]=H_el(Ss[3],Ss[3],Rs[i])+H_rot(ls[j],ls[j],Rs[i],μ)
    end
    lab2=permutedims((x->"l=$(x.l)").(ls))
    plt2=plot(austrip.(Rs),austrip.(vals2),label=lab,legend=:outertopright,
    title="S=2")
    return plot(plt0,plt1,plt2), plt0, plt1, plt2
end

using ProgressMeter
"""Generate plot of H_sd couplings for |lmₗ⟩ states up to lmax"""
function plt_H_sd_coupling(lmax=20)
    lookup=SmS_lookup_generator(lmax)
    n=length(lookup)
    arr=zeros(Float64,n,n)
    for i=1:n, j=1:n
        arr[i,j]=H_sd_coeffs(lookup[i],lookup[j])
    end
    unq=unique(γ_ket_convert.(lookup))
    nᵧ=length(unq)
    arrᵧ=zeros(Float64,nᵧ+2,nᵧ+2)
    for i=1:(nᵧ+2) # first row and col nonzero (delete this in editing software)
        arrᵧ[1,i]=1.0; arrᵧ[i,1]=1.0; arrᵧ[nᵧ,i]=1.0; arrᵧ[i,nᵧ]=1.0
    end
    numpts=Int(nᵧ^2*n^2); prog = Progress(numpts, 1)
    for i=1:nᵧ, j=1:nᵧ
        γ, γ_ = unq[i], unq[j]
        sum=0e0
        for p=1:n, q=1:n # p row, q column
            next!(prog)
            γ_ket_convert(lookup[p])==γ || continue
            γ_ket_convert(lookup[q])==γ_ || continue
            sum += H_sd_coeffs(lookup[p],lookup[q])
        end
        arrᵧ[i+1,j+1] = sum
    end
    spy(arrᵧ,markersize=15,framestyle=:none,legend=nothing,color=:bwr), arrᵧ, unq
end

include("./Modules/interactions/H_sd.jl")
# compare H_sd radial dependence to BO potentials
function plt_H_sd_radial(;Rmin=3.0u"bohr",Rmax=20u"bohr",no_R_pts=10000)
    Rs=LinRange(Rmin,Rmax,no_R_pts)
    # H_el values
    firstS(i::Int)=lookup[findall(x->x.S==i,lookup)[1]] # grabs first state with correct S
    Ss = [firstS(S) for S in 0:2] # prototypical S states
    ellab=permutedims((x->"S=$(x.S)").(Ss))
    elvals=zeros(length(Rs),length(Ss))u"hartree" #initialse
    for i in 1:no_R_pts, j in 1:length(0:2)
        elvals[i,j]=H_el(Ss[j],Ss[j],Rs[i])
    end
    # H_sd values
    SDvals = H_sd_radial.(Rs)
    lab=hcat(["SD"],ellab)
    vals=hcat(SDvals,elvals)
    plot(austrip.(Rs),abs.(austrip.(vals)),yscale=:log10,
    xlabel="R (a₀)",ylabel="Interaction strength (Eₕ)",
    labels=lab)
end


"""Plot Maxwell Boltzmann distribution vs k, given temp"""
function plt_MBD(T::Unitful.Temperature; n=1000, μ=0.5*4.002602u"u",
    kmin=1e-6u"bohr^-1", kmax=1e-3u"bohr^-1")
    ks=LinRange(kmin,kmax,n)
    f(k)=4π*(1u"ħ"/μ)^3*k^2*(μ/(2π*1u"k_au"*T))^(3/2)*exp(-1u"ħ^2"*k^2/(2u"k_au"*μ*T))
    plot(1e3.*austrip.(ks),austrip.(f.(ks)),title="T=$T",
    xlabel="k (10⁻³ a₀⁻¹)", ylabel="f(k) (arbitrary units)",legend=false)
end
