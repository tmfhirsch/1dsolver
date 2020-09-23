# produce plots of the different Interactions in my model, for investigation
# and for thesis writing
# header last updates 23/09/20

using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using StateStructures, Interactions
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots
using BSON, Dates

const lmaxx=20 # more than enough states to be relevant to my calculations
# states up to lmax of 20 (more than enough)
lookup=SmS_lookup_generator(lmaxx)

# Rotational Interaction
function plt_H_rot(;lmax=6,Rmin=1e-10u"bohr",Rmax=1e6u"bohr",no_R_pts=10000,μ=0.5*4.002602u"u")
    firstl(i::Int)=lookup[findall(x->x.l==i,lookup)[1]] # grabs first state with correct l
    ls = [firstl(i) for i in 1:lmax] # prototypical l states
    lab=permutedims((x->"l=$(x.l)").(ls)) # labels for plot
    Rs=LinRange(Rmin,Rmax,no_R_pts)
    vals=zeros(length(Rs),length(ls))u"hartree" #initialse
    for i in 1:no_R_pts, j in 1:lmax
        vals[i,j]=H_rot(ls[j],ls[j],Rs[i],μ)
    end
    plt=plot(austrip.(Rs),austrip.(vals),xscale=:log10,yscale=:log10,
    xlabel="R (a₀)", ylabel="Hᵣₒₜ(l,R) (Eₕ)",label=lab,legend=:outertopright)
end

# Electronic Interaction plotting
function plt_H_el(;Rmin=3.0u"bohr",Rmax=1e2u"bohr",no_R_pts=10000)
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

"""Saves plots of interactions as .svg files in /PHYS4110/Plots/Interactions/"""
function save_interaction_plt(plt::Plots.Plot,title::String)
    prevdir=pwd()
    cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Plots")
    savefig(plt, title*".svg")
    cd(prevdir)
end
