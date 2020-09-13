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
    params_str="_E"*string(ustrip(data.ϵ))*"_B"*string(ustrip(data.B))*"_lmax"*string(data.lmax)
    save_str=date_str*params_str*".pairwiseCS"
    bson(save_str, σ=data.σ, γ_lookup=data.γ_lookup, ϵ=data.ϵ, B=data.B, lmax=data.lmax)
    cd(wd)
end

#checks that save_pairwiseCS runs for 1e-12 Eh, 0.01T, lmax=0
function test_save_pairwiseCS()
    output=σ_matrix(1e-12u"hartree",0.01u"T",0)
    save_pairwiseCS(output)
end

""" generate pairwiseCS for different energies (in Eh), logarithmically spaced
    Inputs: log₁₀(austrip(Emin)), log₁₀(austrip(Emax)), n = number of different energies, lmax;
    B=0T magnetic field strength
    Outputs: / (saves files using save_pairwiseCS)"""
function diff_E_data(Emin_exp,Emax_exp,n::Integer,lmax::Integer;B=0u"T")
    Es=exp10.(LinRange(Emin_exp,Emax_exp,n))u"hartree" # energies
    for E in Es
        output=σ_matrix(E,B,lmax)
        save_pairwiseCS(output)
    end
end

""" loads pairwiseCS datas with filenames that fit
    Inputs: Emin/max ~[E], Bmin/max ~ [Tesla], lmax
    Output: array of σ_outputs found from filename filtering"""
function load_pairwiseCS(Emin::Unitful.Energy,Emax::Unitful.Energy,
    Bmin::Unitful.BField,Bmax::Unitful.BField,lmax::Integer)
    # directory change
    wd=pwd() # current working directory
    cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\PairwiseCS")
    # initialise output
    output=[]
    for f in readdir()
        occursin("pairwiseCS",f) || continue
        fE=let
            SIndex=findfirst(isequal('E'),f)+1
            EIndex=findall(isequal('_'),f)[2]-1
            str=f[SIndex:EIndex]
            parse(Float64,str)u"hartree"
        end
        fB=let
            SIndex=findfirst(isequal('B'),f)+1
            EIndex=findall(isequal('_'),f)[3]-1
            str=f[SIndex:EIndex]
            parse(Float64,str)u"T"
        end
        flmax=let
            Index=findfirst("lmax",f)[end]+1
            str=f[Index]
            parse(Int,str)
        end
        (Emin<=fE<=Emax)&&(Bmin<=fB<=Bmax)&&(flmax==lmax) || continue
        # within parameters, loading:
        load=BSON.load(f)
        data=σ_output(load[:σ],load[:γ_lookup],load[:ϵ],load[:B],load[:lmax])
        push!(output,data)
    end
    cd(wd)
    output
end

# test if load_pairwiseCS runs
test_load_pairwiseCS()=load_pairwiseCS(1e-12u"hartree",1e-8u"hartree",0.0u"T",1.0u"T",0)

# loads, then saves ages all files in pairwiseCSfunction
function pairwiseCS_resaver()
    wd=pwd()
    cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\PairwiseCS")
    for f in readdir()
        load=BSON.load(f)
        data=σ_output(load[:σ],load[:γ_lookup],load[:ϵ],load[:B],load[:lmax])
        rm(f) # delete old file
        save_pairwiseCS(data) # save new file
    end
    cd(wd)
end

# create labels for plot legend, from γ_lookup given, keeping order of the lookup
function label_from_γ_lookup(lookup::Array{γ_ket,1})
    n=length(lookup)
    lab=fill("",1,length(lookup))
    for i in 1:n
        γ=lookup[i]
        lab[i]="S=$(γ.S), mS=$(γ.mS)"
    end
    lab
end

# plot diagonal elements of pairwise σ, i.e. elastic cross sections
function diff_E_plot(Emin,Emax,B,lmax)
    datas=load_pairwiseCS(Emin,Emax,B,B,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.ϵ)) # sort by increasing energy
    γ_lookup=datas[1].γ_lookup
    # check that all data have the same γ_lookup
    @assert all((x->x.γ_lookup==γ_lookup).(datas)) "not all data has same γ_lookup"
    no_γ, no_datas = length(γ_lookup), length(datas)
    σs=zeros(no_datas,no_γ)u"bohr^2" # row = diff. energy, col = diff. state
    for γ in 1:no_γ, d in 1:no_datas
        σs[d,γ]=datas[d].σ[γ,γ] # diag elements ↔ elastic scattering
    end
    ϵs=(x->x.ϵ).(datas)
    plot(austrip.(ϵs),austrip.(σs),xlabel="Energy (Eh)", xscale=:log10,
    ylabel="σ (a₀²)",yscale=:log10, minorticks=true,
    label=label_from_γ_lookup(γ_lookup),legend=:outertopright)
    hline!([4*pi*austrip((7.54u"nm")^2)])
end
