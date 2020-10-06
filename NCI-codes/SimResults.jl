# functions to be include()'d into a julia script for submitting jobs


# investigating cross sections
using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using CrossSections, StateStructures, Interactions
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots
using BSON, Dates

using Distributed # testing Distributed today

const G = 1e-4u"T" # Gauss unit of magnetic flux density

# parameters for ICs/matching
const lhs=3e0u"bohr"; const mid=5e1u"bohr"; const rhs=2e2u"bohr"; const rrhs=1e4u"bohr"

# saves pairwise cross section output in ./Smat_dir/, E in Eh and B in T
function save_output(output::Union{S_output,γ_output,I_output})
    wd=pwd() # current directory, to move back into at the end
    if typeof(output)==S_output
        savedir, filetype = Smat_dir, ".Smat"
    elseif typeof(output)==γ_output
        savedir, filetype = gampwcs_dir, ".gampwcs"
    else # I_output
        savedir, filetype = Ics_dir, ".Ics"
    end
    cd(savedir)
    date_str=string(Dates.today())
    params_str="_E"*string(ustrip(uconvert(u"hartree",output.ϵ)))*"_B"*
    string(ustrip(uconvert(u"T",output.B)))*"_lmax"*string(output.lmax)
    save_str=date_str*params_str*filetype
    bson(save_str, data=output)
    cd(wd)
end

""" generate and save .Smat for different energies (in Eh), logarithmically spaced
    Inputs: log₁₀(austrip(Emin)), log₁₀(austrip(Emax)), n = number of different energies, lmax;
    B=0T magnetic field strength
    Outputs: / (saves files using save_Smat)"""
function gen_diffE_data(Emin_exp,Emax_exp,n::Integer,lmax::Integer;B=0u"T",
    lhs=lhs,mid=mid,rhs=rhs,rrhs=rrhs)
    Es=exp10.(LinRange(Emin_exp,Emax_exp,n))u"hartree" # energies
    println("lmax=$lmax, B=$B. Generating S_output for E/Eh= ")
    for E in Es
        println("$(austrip(E)), ")
        output=S_matrix(E,B,lmax,lhs,mid,rhs,rrhs)
        save_output(output)
    end
end

""" generate and save .Smat for different magnetic field strengths (in T), linearly spaced,
    with constant E.
    Inputs: Bmin~[BField], Bmax~[BField], n = number of different Bfields, E~[E], lmax;
    lhs~[L]=3aₒ, mid~[L]=50aₒ, rhs~[L]=200aₒ, rrhs~[L]=1e4aₒ
    Outputs: / (saves files using save_pairwiseCS)"""
function gen_diffB_constE_data(Bmin::Unitful.BField,Bmax::Unitful.BField,n::Integer,E::Unitful.Energy,lmax::Integer;
    lhs=lhs,mid=mid,rhs=rhs,rrhs=rrhs)
    Bs=LinRange(Bmin,Bmax,n) # energies
    println("lmax=$lmax, E=$E. Generating σ_output for B/1T= ")
    for B in Bs
        println("$(ustrip(uconvert(u"T",B))), ")
        output=S_matrix(E,B,lmax,lhs,mid,rhs,rrhs)
        save_output(output)
    end
end

# asymptotic wavenumber as function of state, energy, Bfield
function k∞(ket::Union{SmS_ket,γ_ket},ϵ::Unitful.Energy,B::Unitful.BField;
    μ::Unitful.Mass=0.5*4.002602u"u")
    V∞ = H_zee(ket,ket,B)
    k = auconvert(sqrt(Complex(2)*μ*(ϵ-V∞))/1u"ħ")
    k
end

# asymptotic energy as function of state, wavenumber, Bfield
function E∞(ket::Union{SmS_ket,γ_ket},k::typeof(0e0u"bohr^-1"),B::Unitful.BField;
    μ::Unitful.Mass=0.5*4.002602u"u")
    uconvert(u"hartree", 1u"ħ^2"*k^2/(2*μ) + H_zee(ket,ket,B))
end


"""Generates .Smat data for different B fields for the same k(γ), iterating over all γ
    Inputs: Bmin, Bmax, n, k, lmax
    Output: / saves output, iterating over all B, for all l=0⟺S∈{0,2} γ kets for given k"""
function gen_diffB_constk_data(Bmin::Unitful.BField,Bmax::Unitful.BField,n::Integer,
    k::typeof(0e0u"bohr^-1"), lmax::Int; μ::Unitful.Mass=0.5*4.002602u"u",
    lhs=lhs,mid=mid,rhs=rhs,rrhs=rrhs)
    Bs=LinRange(Bmin,Bmax,n)
    l0γs = let
        allγ = γ_lookup_generator() # all S=0,1,2 γ states
        allγ[findall(x->x.S!=1,allγ)] # no S=1 because zeroing out l>0 states
    end
    println("Generating outputs for k=$k, lmax=$lmax")
    @sync @distributed for B in Bs
        for γ in l0γs
            E = E∞(γ,k,B)
            println("B=$B, γ=$γ")
            preexisting=load_data("S",E,E,B,B,lmax)
            if length(preexisting)>0 # skip data that has already been generated
                continue
            end
            output=S_matrix(E,B,lmax,lhs,mid,rhs,rrhs)
            save_output(output)
        end
    end
end

# loads, then resaves all files in Smat_dir, gampwcs_dir, Ics_dir (for changing name schemes)
function resaver()
    wd=pwd()
    cd(Smat_dir)
    for f in readdir()
        data=BSON.load(f)[:data]
        #rm(f) # delete old file
        save_output(data) # save new file
    end
    cd(gampwcs_dir)
    for f in readdir()
        data=BSON.load(f)[:data]
        #rm(f) # delete old file
        save_output(data) # save new file
    end
    cd(Ics_dir)
    for f in readdir()
        data=BSON.load(f)[:data]
        #rm(f) # delete old file
        save_output(data) # save new file
    end
    cd(wd)
end

""" loads .Smat file/s, converts to γ_output and saves as .gampwcs file
    Input: target::String, complete filename of target file in /SmS-PairwiseCS.
    If an exact match is given to a target file, converts just that file.
    Otherwise, converts all filenames where occursin(target, filename).
        This function defaults to matching/converting every file"""
function S2γ_data(target="")
    wd=pwd()
    cd(Smat_dir)
    Sfiles=filter((x->occursin(".Smat",x)), readdir())
    for f in Sfiles # first check for exact matches; break at first match
        if target==f
            datum=BSON.load(f)[:data]
            @assert typeof(datum)==S_output "$f is not a S_output"
            gamdatum=S2γ_output(datum)
            save_output(gamdatum)
            return nothing # stop after finding this exactly matching file
        end
    end
    for f in Sfiles # no exact match found; converting all applicable files
        if occursin(target,f)
            datum=BSON.load(f)[:data]
            @assert typeof(datum)==S_output "$f is not a S_output"
            gamdatum=S2γ_output(datum)
            save_output(gamdatum)
        end
    end
    cd(wd)
end

""" loads .Smat file/s, converts to I_output and saves as .Ics file
    Input: target::String, complete filename of target file in /Smat_dir.
    If an exact match is given to a target file, converts just that file.
    Otherwise, converts all filenames where occursin(target, filename).
        This function defaults to matching/converting every file"""
function S2I_data(target="")
    wd=pwd()
    cd(Smat_dir)
    Sfiles=filter((x->occursin(".Smat",x)), readdir())
    for f in Sfiles # first check for exact matches; break at first match
        if target==f
            datum=BSON.load(f)[:data]
            @assert typeof(datum)==S_output "$f is not a S_output"
            Idatum=S2I_output(datum)
            save_output(Idatum)
            return nothing # stop after finding this exactly matching file
        end
    end
    for f in Sfiles # no exact match found; converting all applicable files
        if occursin(target,f)
            datum=BSON.load(f)[:data]
            @assert typeof(datum)==S_output "$f is not a S_output"
            Idatum=S2I_output(datum)
            save_output(Idatum)
        end
    end
    cd(wd)
end

###################################Load data####################################
""" load .Smat data with filenames fitting bounds on E, B, lmax
    Inputs: flag∈{"S","gam"}, Emin/max ~[E], Bmin/max ~ [Tesla], lmax
    Output: array of outputs found from filename filtering"""
function load_data(flag::String,Emin::Unitful.Energy,Emax::Unitful.Energy,
    Bmin::Unitful.BField,Bmax::Unitful.BField,lmax::Integer)
    @assert flag in ["S","gam","I"] "flag is neither 'S' or 'gam' or 'I'"
    # directory change
    wd=pwd() # current working directory
    if flag=="S"
        filetype, dir = ".Smat", Smat_dir
    elseif flag=="gam"
        filetype, dir = ".gampwcs", gampwcs_dir
    else
        filetype, dir = ".Ics", Ics_dir
    end
    cd(dir)
    # initialise output
    output=[]
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
            SIndex=findfirst("lmax",f)[end]+1
            EIndex=findfirst(filetype,f)[1]-1
            str=f[SIndex:EIndex]
            parse(Int,str)
        end
        (Emin<=fE<=Emax)&&(Bmin<=fB<=Bmax)&&(flmax==lmax) || continue
        # within parameters, loading:
        datum=BSON.load(f)[:data]
        push!(output,datum)
    end
    cd(wd)
    output
end
