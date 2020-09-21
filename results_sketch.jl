# investigating cross sections
using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using CrossSections, StateStructures, Interactions
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots
using BSON, Dates

const SmSpwcs_dir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\21-9-20-tests\c"
const gampwcs_dir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\21-9-20-tests\c"
# saves pairwise cross section output in ./SmSpwcs_dir/, E in Eh and B in T
function save_SmSpwcs(data::σ_output)
    wd=pwd() # current directory, to move back into at the end
    cd(SmSpwcs_dir)
    date_str=string(Dates.  today())
    params_str="_E"*string(ustrip(uconvert(u"hartree",data.ϵ)))*"_B"*
    string(ustrip(uconvert(u"T",data.B)))*"_lmax"*string(data.lmax)
    save_str=date_str*params_str*".SmSpwcs"
    bson(save_str, σ=data.σ, lookup=data.lookup, ϵ=data.ϵ, B=data.B, lmax=data.lmax)
    cd(wd)
end

#checks that save_pairwiseCS runs for 1e-12 Eh, 0.01T, lmax=0
function test_save_SmSpwcs(;ϵ=1e-12u"hartree",B=0u"T",lmax=0)
    output=σ_matrix(ϵ,B,lmax)
    save_SmSpwcs(output)
end

""" generate and save .SmSpwcs for different energies (in Eh), logarithmically spaced
    Inputs: log₁₀(austrip(Emin)), log₁₀(austrip(Emax)), n = number of different energies, lmax;
    B=0T magnetic field strength
    Outputs: / (saves files using save_pairwiseCS)"""
function gen_diffE_data(Emin_exp,Emax_exp,n::Integer,lmax::Integer;B=0u"T",
    lhs=3.0u"bohr",mid=50.0u"bohr",rhs=1000.0u"bohr")
    Es=exp10.(LinRange(Emin_exp,Emax_exp,n))u"hartree" # energies
    println("lmax=$lmax. Generating σ_output for E/Eh= ")
    for E in Es
        println("$(austrip(E)), ")
        output=σ_matrix(E,B,lmax)
        save_SmSpwcs(output)
    end
end

# loads, then resaves all files in pairwiseCSfunction; for changing name schemes
function SmSpwcs_resaver()
    wd=pwd()
    cd(SmSpwcs_dir)
    for f in readdir()
        load=BSON.load(f)
        data=σ_output(load[:σ],load[:lookup],load[:ϵ],load[:B],load[:lmax])
        #rm(f) # delete old file
        save_SmSpwcs(data) # save new file
    end
    cd(wd)
end

#################################.gampwcs functions#############################
# saves γ_output as a .gampwcs file in ./gampwcs_dir/, E in Eh and B in T
function save_gampwcs(data::γ_output)
    wd=pwd()
    cd(gampwcs_dir)
    date_str=string(Dates.today())
    params_str="_E"*string(ustrip(uconvert(u"hartree",data.ϵ)))*"_B"*
    string(ustrip(uconvert(u"T",data.B)))*"_lmax"*string(data.lmax)
    save_str=date_str*params_str*".gampwcs"
    bson(save_str, σ=data.σ, lookup=data.γ_lookup, ϵ=data.ϵ, B=data.B, lmax=data.lmax)
    cd(wd)
end

#checks that save_pairwiseCS runs for 1e-12 Eh, 0.01T, lmax=0
function test_save_gampwcs(;ϵ=1e-11u"hartree",B=0.01u"T",lmax=0)
    output=σ2γ_output(σ_matrix(ϵ,B,lmax))
    save_gampwcs(output)
end

""" loads .SmSpwcs file/s, converts to γ_output and saves as .gampwcs file
    Input: target::String, complete filename of target file in /SmS-PairwiseCS.
    If an exact match is given to a target file, converts just that file.
    Otherwise, converts all filenames where occursin(target, filename).
        This function defaults to matching/converting every file"""
function σ2γ_data(target="")
    wd=pwd()
    cd(SmSpwcs_dir)
    SmS_files=filter((x->occursin(".SmSpwcs",x)), readdir())
    for f in SmS_files # first check for exact matches; break at first match
        if target==f
            load=BSON.load(f)
            SmSdata=σ_output(load[:σ],load[:lookup],load[:ϵ],load[:B],load[:lmax])
            gamdata=σ2γ_output(SmSdata)
            save_gampwcs(gamdata)
            return
        end
    end
    for f in SmS_files # no exact match found; converting all applicable files
        if occursin(target,f)
            load=BSON.load(f)
            SmSdata=σ_output(load[:σ],load[:lookup],load[:ϵ],load[:B],load[:lmax])
            gamdata=σ2γ_output(SmSdata)
            save_gampwcs(gamdata)
        end
    end
    cd(wd)
end

###################################Load data####################################
""" load .SmSpwcs data with filenames fitting bounds on E, B, lmax
    Inputs: flag∈{"SmS","gam"}, Emin/max ~[E], Bmin/max ~ [Tesla], lmax
    Output: array of outputs found from filename filtering"""
function load_data(flag::String,Emin::Unitful.Energy,Emax::Unitful.Energy,
    Bmin::Unitful.BField,Bmax::Unitful.BField,lmax::Integer)
    @assert flag in ["SmS","gam"] "flag is neither 'SmS' or 'gam'"
    # directory change
    wd=pwd() # current working directory
    dir = flag=="SmS" ? SmSpwcs_dir : gampwcs_dir
    cd(dir)
    # initialise output
    output=[]
    filetype = flag=="SmS" ? ".SmSpwcs" : ".gampwcs"
    for f in readdir()
        occursin(filetype,f) || continue
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
        if flag=="SmS"
            data=σ_output(load[:σ],load[:lookup],load[:ϵ],load[:B],load[:lmax])
        else
            data=γ_output(load[:σ],load[:lookup],load[:ϵ],load[:B],load[:lmax])
        end
        push!(output,data)
    end
    cd(wd)
    output
end

# test if load_pairwiseCS runs
test_load_SmSpwcs()=load_data("SmS",1e-12u"hartree",1e-8u"hartree",0.0u"T",1.0u"T",0)
test_load_gampwcs()=load_data("gam",1e-12u"hartree",1e-8u"hartree",0.0u"T",1.0u"T",0)

###################################Plotting functions###########################

# create labels for plot legend, from γ_lookup given, keeping order of the lookup
function label_from_lookup(lookup::Union{Array{γ_ket,1}, Array{SmS_ket,1}})
    n=length(lookup)
    lab=fill("",1,length(lookup))
    if typeof(lookup)==Array{γ_ket,1} # converting γ_kets: S, mS relevant
        for i in 1:n
            ket=lookup[i]
            lab[i]="S=$(ket.S), mS=$(ket.mS)"
        end
    else # converting SmS_kets; S, mS, l, ml relevant
        for i in 1:n
            ket=lookup[i]
            lab[i]="S=$(ket.S), mS=$(ket.mS), l=$(ket.l), ml=$(ket.ml)"
        end
    end
    lab
end

# plot diagonal elements of pairwise σ, i.e. elastic cross sections
function diffE_gam_plot(Emin::Unitful.Energy,Emax::Unitful.Energy,B::Unitful.BField,lmax::Int)
    datas=load_data("gam",Emin,Emax,B,B,lmax)
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
    label=label_from_lookup(γ_lookup),legend=:outertopright)
    hline!([4*pi*austrip((7.54u"nm")^2)])
end
