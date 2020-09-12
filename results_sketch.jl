# investigating cross sections
using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using CrossSections, StateStructures, Interactions
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots


using BSON, Dates

# saves pairwise cross section data. STRIPS UNITS FROM FILENAME, ASSUMES
# E in Eh and B field in Tesla
function save_pairwiseCS(data::σ_output)
    wd=pwd() # current directory, to move back into at the end
    cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\PairwiseCS")
    date_str=string(Dates.today())
    params_str="_ϵ"*string(round(ustrip(data.ϵ),digits=3))*"_B-"*string(round(ustrip(data.B),digits=3))*"_lmax"*string(data.lmax)
    save_str="pairwiseCS_"*date_str*params_str
    bson(save_str, σ=data.σ, γ_lookup=data.γ_lookup, ϵ=data.ϵ, B=data.B, lmax=data.lmax)
    cd(wd)
end

#checks that save_pairwiseCS runs for 1e-12 Eh, 0.01T, lmax=0
function test_save_pairwiseCS()
    output=σ_matrix(1e-12u"hartree",0.01u"T",0)
    save_pairwiseCS(output)
end

""" generate pairwiseCS for different energies, logarithmically spaced
    Inputs: log₁₀(Emin), log₁₀(Emax), n = number of different energies, lmax;
    B=0T magnetic field strength
    Outputs: / (saves files using save_pairwiseCS)"""
function diff_E_data(Emin,Emax,n::Integer,lmax::Integer;B=0u"T")
    Es=exp10.(LinRange(Emin,Emax,n))u"hartree" # energies
    for E in Es
        output=σ_matrix(E,B,lmax)
        save_pairwiseCS(output)
    end
end
