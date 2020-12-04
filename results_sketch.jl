# investigating cross sections
using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using CrossSections, StateStructures, Interactions
using Unitful, UnitfulAtomic, LinearAlgebra
using Plots, Plots.PlotMeasures
using BSON, Dates
using Statistics

using Distributed

const G = 1e-4u"T" # Gauss unit of magnetic flux density

const Smat_dir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\zeroB\Smat"
const gampwcs_dir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\zeroB\gampwcs"
const Ics_dir=raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Results\zeroB\Ics"

# parameters for ICs/matching
const lhs=3e0u"bohr"; const mid=5e1u"bohr"; const rhs=2e2u"bohr"; const rrhs=1e3u"bohr"

# convert wavenumber to corresponding temperature by 3/2 kᵦT
function corresponding_temp(k::typeof(0e0u"bohr^-1");μ=0.5*4.002602u"u")
    KE=1u"ħ^2"*k^2/(2*μ)
    T=2/(3u"k_au")*KE
    uconvert(u"μK",T)
end

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

#checks that save_pairwiseCS runs for 1e-12 Eh, 0.01T, lmax=0
function test_save_output(;ϵ=1e-12u"hartree",B=0u"T",lmax=0,
    lhs=lhs,mid=mid,rhs=rhs,rrhs=rrhs)
    output=S_matrix(ϵ,B,lmax,lhs,mid,rhs,rrhs)
    save_output(output)
    save_output(S2γ_output(output))
    save_output(S2I_output(output))
end

""" generate and save .Smat for different energies (in Eh), logarithmically spaced
    Inputs: log₁₀(austrip(Emin)), log₁₀(austrip(Emax)), n = number of different energies, lmax;
    B=0T magnetic field strength
    Outputs: / (saves files using save_Smat)"""
function gen_diffE_data(Emin_exp,Emax_exp,n::Integer,lmax::Integer;B=0u"T",
    lhs=lhs,mid=mid,rhs=rhs,rrhs=rrhs)
    Es=exp10.(LinRange(Emin_exp,Emax_exp,n))u"hartree" # energies
    println("lmax=$lmax, B=$B. Generating S_output for E/Eh= ")
    @sync @distributed for E in Es
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
    @sync @distributed for B in Bs
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
#=
#################################.gampwcs functions#############################
# saves γ_output as a .gampwcs file in ./gampwcs_dir/, E in Eh and B in Tk=
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
=#
#=
#checks that save_pairwiseCS runs for 1e-12 Eh, 0.01T, lmax=0
function test_save_gampwcs(;ϵ=1e-11u"hartree",B=0.01u"T",lmax=0)
    output=S2γ_output(S_matrix(ϵ,B,lmax))
    save_gampwcs(output)
end
=#
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
    Inputs: flag∈{"S","gam","I"}, Emin/max ~[E], Bmin/max ~ [Tesla], lmax
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

# test if load_pairwiseCS runs
test_load_Smat()=load_data("S",1e-12u"hartree",1e-8u"hartree",0.0u"T",1.0u"T",0)
test_load_gampwcs()=load_data("gam",1e-12u"hartree",1e-8u"hartree",0.0u"T",1.0u"T",0)
test_load_Ics()=load_data("I",1e-12u"hartree",1e-8u"hartree",0.0u"T",1.0u"T",0)


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

# finds unique kets in an array of S or γ outputs
function unqkets(datas)
    @assert all(x->typeof(x)==S_output,datas) ||
        all(x->typeof(x)==γ_output,datas) ||
        all(x->typeof(x)==I_output,datas) "datas not a list of S_outputs, γ_outputs, or I_outputs"
    @assert length(datas)>0 "unqkets is given an empty list of data"
    if typeof(datas[1])==S_output
        Sunq::Array{SmS_ket,1}=[]
        for d in datas # search all data
            for k in d.lookup
                if !(k in Sunq)  # if data has a as-yet-unseen γ_ket
                    push!(Sunq,k)
                end
            end
        end
        return Sunq
    else # γ_output data
        γunq::Array{γ_ket,1}=[]
        for d in datas # search all data
            for γ in d.γ_lookup
                if !(γ in γunq)  # if data has a as-yet-unseen γ_ket
                    push!(γunq,γ)
                end
            end
        end
        return γunq
    end
end

###################################Elastic σ plots##############################

# plot diagonal elements of pairwise σ (i.e. elastic cross sections) vs Energy
function diffE_gam_plot(Emin::Unitful.Energy,Emax::Unitful.Energy,B::Unitful.BField,lmax::Int)
    datas=load_data("gam",Emin,Emax,B,B,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.ϵ)) # sort by increasing energy
    unq = unqkets(reverse(datas)) # reverse to find uniques starting with high energy data
    pltdata=[] # array to store ([ϵ],[σ}) pairs for the different γ
    pltlabel=label_from_lookup(unq)
    for γ in unq
        datatuple::typeof(([0.0u"hartree"],[0.0u"bohr^2"]))=([],[])
        for d in datas # already sorted datas by energy
            if γ in d.γ_lookup
                dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
                push!(datatuple[1],d.ϵ) # store energy
                push!(datatuple[2],d.σ[dγ_index,dγ_index]) # store cross section
            end
        end
        push!(pltdata,datatuple)
    end
    plot(austrip.(pltdata[1][1]),austrip.(pltdata[1][2]),xlabel="Energy (Eh)", xscale=:log10,
    ylabel="σ (a₀²)",yscale=:log10, minorticks=true, label=pltlabel[1], legend=:outertopright,
    title="B=$B, lmax=$lmax",
    left_margin=5mm,bottom_margin=5mm,top_margin=5mm,
    linewidth=2,grid=false)
    if length(pltdata)>1 # plot rest of the γ_ket series
        for i in 2:length(pltdata)
            plot!(austrip.(pltdata[i][1]),austrip.(pltdata[i][2]),label=pltlabel[i],
            linewidth=2)
        end
    end
    hline!([4*pi*austrip((7.54u"nm")^2)],label="S=2 4πa²", linewidth=2)
end

"""
Plot diagonal elements of pairwise σ (i.e. elastic cross sections) vs wavenumber
Inputs: kmin~[L]⁻¹, kmax~[L]⁻¹ for plot x axis, B~[BField],lmax;
    Emin~[E]=0Eₕ, Emax~[E]=1Eₕ for finding appropriate data
Output: plot of elastic cross sections vs k"""
function diffk_gam_plot(kmin::typeof(0e0u"bohr^-1"),kmax::typeof(0e0u"bohr^-1"),
    B::Unitful.BField,lmax::Int)
    # load all data with correct B
    datas=load_data("gam",-(Inf)u"hartree",(Inf)u"hartree",B,B,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.ϵ)) # sort by increasing energy
    unq = unqkets(reverse(datas)) # reverse to find uniques starting with high energy data
    pltdata=[] # array to store ([k],[σ}) pairs for the different γ
    pltlabel=label_from_lookup(unq)
    for γ in unq
        datatuple::typeof(([0.0u"bohr^-1"],[0.0u"bohr^2"]))=([],[])
        for d in datas # already sorted datas by energy
            if γ in d.γ_lookup
                dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
                k=k∞(γ,d.ϵ,B) # the asymptotic wavenumber for this particular data
                imag(k)==0.0u"bohr^-1" || continue # don't store if the wavenumber is complex⟺channel closed
                kmin <= real(k) <= kmax || continue # don't store if the wavenumber is out of bounds
                push!(datatuple[1],k) # store wavenumber
                push!(datatuple[2],d.σ[dγ_index,dγ_index]) # store cross section
            end
        end
        push!(pltdata,datatuple)
    end
    println("Minimum values are:")
    println("S, mS = ",(x->"$(x.S), $(x.mS)").(unq))
    println(austrip.([pltdata[i][2][1] for i=1:length(pltdata)]))
    # plot first γ_ket
    plot(austrip.(pltdata[1][1]),austrip.(pltdata[1][2]),xlabel="Wavenumber (a₀⁻¹)", xscale=:log10,
    ylabel="σ (a₀²)",yscale=:log10, minorticks=true, label=pltlabel[1], legend=:outertopright,
    title="B=$B, lmax=$lmax",
    left_margin=5mm,bottom_margin=5mm, top_margin=5mm,
    linewidth=2, grid=false)
    if length(pltdata)>1 # plot rest of the γ_kets
        for i in 2:length(pltdata)
            plot!(austrip.(pltdata[i][1]),austrip.(pltdata[i][2]),label=pltlabel[i])
        end
    end
    hline!([4*pi*austrip((7.54u"nm")^2)],label="S=2 4πa²")
end

"""
    Plot diagonal elements of pairwise σ (i.e. elastic cross sections) vs Bfield
    at constant wavenumber
Inputs: Bmin~[BField], Bmax~[BField], k~[L]⁻¹ for plot x axis, lmax
Output: plot of elastic cross sections vs B"""
function diffB_gam_plot(Bmin::Unitful.BField, Bmax::Unitful.BField,
    k::typeof(0e0u"bohr^-1"),lmax::Integer)
    # load all data with correct B
    datas=load_data("gam",-(Inf)u"hartree",(Inf)u"hartree",Bmin,Bmax,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.B)) # sort by increasing Bfield
    unq = unqkets(reverse(datas)) # reverse to find uniques starting with high energy data
    pltdata=[] # array to store ([B],[σ}) pairs for the different γ
    pltlabel=label_from_lookup(unq)
    for γ in unq
        γdata::typeof(([0.0u"T"],[0.0u"bohr^2"]))=([],[]) # data for this γ
        for d in datas # already sorted datas by energy
            if γ in d.γ_lookup # in case γ is closed in this data
                dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
                # need match of energy to B field for this γ
                d.ϵ==E∞(γ,k,d.B) || continue
                push!(γdata[1],d.B) # store Bfield
                push!(γdata[2],d.σ[dγ_index,dγ_index]) # store el cross section
            end
        end
        push!(pltdata,γdata)
    end
    #=# for plotting individual mS states
    target_mS=2
    msindex=findall(x->x.S==2&&x.mS==target_mS,unq)[1]
    baseline_σ=median(austrip.(pltdata[msindex][2]))
    target_plt=plot(ustrip.(uconvert.(u"T",pltdata[msindex][1])),austrip.(pltdata[msindex][2]),
    ylabel="σ (a₀²)", xlabel="B (T)", title="elastic, |S=2, mₛ=$(target_mS)⟩",
    ylims=(baseline_σ-10,baseline_σ+10))
    target_plt=#
    # for plotting each channel in a subplot
    println("Plotting over $(length(pltdata[1][1])) different B fields")
    # plot S=0
    S0index = findall(x->x.S==0,unq)[1] # index of the S=0 ket
    S0base=median(austrip.(pltdata[S0index][2]))
    pltS0=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S0index][1])),austrip.(pltdata[S0index][2]),
    minorticks=true, label=pltlabel[S0index],legend=false,
    xlabel="B (G)", ylabel="σ (a₀²)",
    ylims=(S0base-10,S0base+10))
    #plot S=2,mS=-2
    S2mSm2index = findall(x->x.S==2&&x.mS==-2,unq)[1] # index of the S=2,mS=-2 ket
    S2mSm2base=median(austrip.(pltdata[S2mSm2index][2]))
    pltS2mSm2=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mSm2index][1])),austrip.(pltdata[S2mSm2index][2]),
    minorticks=true, label=pltlabel[S2mSm2index],legend=false,
    xlabel="B (G)", ylabel="σ (a₀²)",
    ylims=(S2mSm2base-1000,S2mSm2base+1000))
    # plot S=2,mS=-1
    S2mSm1index = findall(x->x.S==2&&x.mS==-1,unq)[1] # index of the S=2,mS=-2 ket
    S2mSm1base=median(austrip.(pltdata[S2mSm1index][2]))
    pltS2mSm1=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mSm1index][1])),austrip.(pltdata[S2mSm1index][2]),
    minorticks=true, label=pltlabel[S2mSm1index],legend=false,
    ylabel="σ (a₀²)", xlabel="B (G)",
    ylims=(S2mSm1base-1000,S2mSm1base+1000))
    # plot S=2,mS=0
    S2mS0index = findall(x->x.S==2&&x.mS==0,unq)[1] # index of the S=2,mS=-2 ket
    S2mS0base=median(austrip.(pltdata[S2mS0index][2]))
    pltS2mS0=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mS0index][1])),austrip.(pltdata[S2mS0index][2]),
    minorticks=true, label=pltlabel[S2mS0index],legend=false,
    ylabel="σ (a₀²)", xlabel="B (G)",
    ylims=(S2mS0base-1000,S2mS0base+1000))
    # plot S=2,mS=1
    S2mS1index = findall(x->x.S==2&&x.mS==1,unq)[1] # index of the S=2,mS=-2 ket
    S2mS1base=median(austrip.(pltdata[S2mS1index][2]))
    pltS2mS1=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mS1index][1])),austrip.(pltdata[S2mS1index][2]),
    minorticks=true, label=pltlabel[S2mS1index],legend=false,
    ylabel="σ (a₀²)", xlabel="B (G)",
    ylims=(S2mS1base-1000,S2mS1base+1000))
    # plot S=2,mS=2
    S2mS2index = findall(x->x.S==2&&x.mS==2,unq)[1] # index of the S=2,mS=-2 ket
    S2mS2base=median(austrip.(pltdata[S2mS2index][2]))
    pltS2mS2=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mS2index][1])),austrip.(pltdata[S2mS2index][2]),
    minorticks=true, label=pltlabel[S2mS2index],legend=false,
    ylabel="σ (a₀²)", xlabel="B (G)",
    ylims=(S2mS2base-1000,S2mS2base+1000))
    return pltS0, pltS2mSm2, pltS2mSm1, pltS2mS0, pltS2mS1, pltS2mS2
    #=# for plotting S=2 on top of each other
    S2indices=filter(x->x!=S0index, 1:length(unq))
    pltS2=plot(ustrip.(uconvert.(u"T",pltdata[S2indices[1]][1])),austrip.(pltdata[S2indices[1]][2]),
    xlabel="B (T)", ylabel="σ (a₀²)", minorticks=true, label=pltlabel[S2indices[1]], legend=:outertopright,
    ylims=(252000,261000))
    if length(S2indices)>1
        for i=2:length(S2indices)
            plot!(ustrip.(uconvert.(u"T",pltdata[S2indices[i]][1])),austrip.(pltdata[S2indices[i]][2]),
            label=pltlabel[S2indices[i]], legend=:outertopright)
        end
    end
    #hline!([4*pi*austrip((7.54u"nm")^2)],label="4π×(7.54nm)²") # S=2 theoretical σ
    # merge plots
    plot(pltS0, pltS2, layout=(2,1),title=["elastic, k=$k, lmax=$lmax" ""],
    left_margin=5mm,bottom_margin=5mm,top_margin=5mm,linewidth=2,grid=false)=#
end

###########################Ionisation σ plots###################################
"""Plot ionisation cross sections vs wavenumber"""
function diffk_I_plot(kmin::typeof(0e0u"bohr^-1"),kmax::typeof(0e0u"bohr^-1"),
    B::Unitful.BField, lmax::Integer)
    # load all data with correct B
    datas=load_data("I",-(Inf)u"hartree",(Inf)u"hartree",B,B,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.ϵ)) # sort by increasing energy
    unq = unqkets(reverse(datas)) # reverse to find uniques starting with high energy data
    pltdata=[] # array to store ([k],[σ}) pairs for the different γ
    pltlabel=label_from_lookup(unq)
    for γ in unq
        datatuple::typeof(([0.0u"bohr^-1"],[0.0u"bohr^2"]))=([],[])
        for d in datas
            if γ in d.γ_lookup
                dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
                k=k∞(γ,d.ϵ,B) # the asymptotic wavenumber for this particular data
                imag(k)==0.0u"bohr^-1" || continue # don't store if the wavenumber is complex⟺channel closed
                kmin <= real(k) <= kmax || continue # don't store if the wavenumber is out of bounds
                push!(datatuple[1],k) # store wavenumber
                push!(datatuple[2],d.σ[dγ_index]) # store cross section
            end
        end
        push!(pltdata,datatuple)
    end
    println("Minimum values are:")
    println("S, mS = ",(x->"$(x.S), $(x.mS)").(unq))
    println(austrip.([pltdata[i][2][1] for i=1:length(pltdata)]))
    # plot first γ_ket
    plt=plot(austrip.(pltdata[1][1]),austrip.(pltdata[1][2]),xlabel="Wavenumber (a₀⁻¹)", xscale=:log10,
    ylabel="σ(PI) (a₀²)", minorticks=true, label=pltlabel[1], legend=:outertopright,
    title="B=$B, lmax=$lmax",
    left_margin=5mm,bottom_margin=5mm,top_margin=5mm,
    linewidth=2,grid=false,
    yscale=:log10)
    if length(pltdata)>1 # plot rest of the γ_kets
        for i in 2:length(pltdata)
            plot!(austrip.(pltdata[i][1]),austrip.(pltdata[i][2]),label=pltlabel[i])
        end
    end
    plt
end

"""
    Plot ionisation σ vs Bfield at constant wavenumber
Inputs: Bmin~[BField], Bmax~[BField], k~[L]⁻¹ for plot x axis, lmax
Output: plot of ionisation cross sections vs B"""
function diffB_I_plot(Bmin::Unitful.BField, Bmax::Unitful.BField,
    k::typeof(0e0u"bohr^-1"),lmax::Integer)
    # load all data with correct B
    datas=load_data("I",-(Inf)u"hartree",(Inf)u"hartree",Bmin,Bmax,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.B)) # sort by increasing Bfield
    unq = unqkets(reverse(datas)) # reverse to find uniques starting with high energy data
    pltdata=[] # array to store ([B],[σ}) pairs for the different γ
    pltlabel=label_from_lookup(unq)
    for γ in unq
        γdata::typeof(([0.0u"T"],[0.0u"bohr^2"]))=([],[]) # data for this γ
        for d in datas # already sorted datas by energy
            if γ in d.γ_lookup # in case γ is closed in this data
                dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
                # need match of energy to B field for this γ
                d.ϵ==E∞(γ,k,d.B) || continue
                push!(γdata[1],d.B) # store Bfield
                push!(γdata[2],d.σ[dγ_index]) # store ionisation cross section
            end
        end
        push!(pltdata,γdata)
    end
    println("Plotting over $(length(pltdata[1][1])) different B fields")
    # for plotting all individually
    # for plotting each channel in a subplot
    println("Plotting over $(length(pltdata[1][1])) different B fields")
    # plot S=0
    S0index = findall(x->x.S==0,unq)[1] # index of the S=0 ket
    S0base=median(austrip.(pltdata[S0index][2]))
    pltS0=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S0index][1])),austrip.(pltdata[S0index][2]),
    minorticks=true, label=pltlabel[S0index], legend=false,
    xlabel="B (G)", ylabel="σ(PI) (a₀²)",
    ylims=(S0base-2000,S0base+2000))
    #plot S=2,mS=-2
    S2mSm2index = findall(x->x.S==2&&x.mS==-2,unq)[1] # index of the S=2,mS=-2 ket
    S2mSm2base=median(austrip.(pltdata[S2mSm2index][2]))
    pltS2mSm2=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mSm2index][1])),austrip.(pltdata[S2mSm2index][2]),
    minorticks=true, label=pltlabel[S2mSm2index], legend=false,
    xlabel="B (G)", ylabel="σ(PI) (a₀²)",
    ylims=(-10,1000))
    # plot S=2,mS=-1
    S2mSm1index = findall(x->x.S==2&&x.mS==-1,unq)[1] # index of the S=2,mS=-2 ket
    S2mSm1base=median(austrip.(pltdata[S2mSm1index][2]))
    pltS2mSm1=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mSm1index][1])),austrip.(pltdata[S2mSm1index][2]),
    minorticks=true, label=pltlabel[S2mSm1index], legend=false,
    ylabel="σ(PI) (a₀²)", xlabel="B (G)",
    ylims=(-10,1000))
    # plot S=2,mS=0
    S2mS0index = findall(x->x.S==2&&x.mS==0,unq)[1] # index of the S=2,mS=-2 ket
    S2mS0base=median(austrip.(pltdata[S2mS0index][2]))
    pltS2mS0=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mS0index][1])),austrip.(pltdata[S2mS0index][2]),
    minorticks=true, label=pltlabel[S2mS0index], legend=false,
    ylabel="σ(PI) (a₀²)", xlabel="B (G)",
    ylims=(-Inf,100))
    # plot S=2,mS=1
    S2mS1index = findall(x->x.S==2&&x.mS==1,unq)[1] # index of the S=2,mS=-2 ket
    S2mS1base=median(austrip.(pltdata[S2mS1index][2]))
    pltS2mS1=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mS1index][1])),austrip.(pltdata[S2mS1index][2]),
    minorticks=true, label=pltlabel[S2mS1index], legend=false,
    ylabel="σ(PI) (a₀²)", xlabel="B (G)",
    ylims=(-10,100))
    # plot S=2,mS=2
    S2mS2index = findall(x->x.S==2&&x.mS==2,unq)[1] # index of the S=2,mS=-2 ket
    S2mS2base=median(austrip.(pltdata[S2mS2index][2]))
    pltS2mS2=plot(10_000 .* ustrip.(uconvert.(u"T",pltdata[S2mS2index][1])),austrip.(pltdata[S2mS2index][2]),
    minorticks=true, label=pltlabel[S2mS2index], legend=false,
    ylabel="σ(PI) (a₀²)", xlabel="B (G)",
    ylims=(-10,100))
    return pltS0, pltS2mSm2, pltS2mSm1, pltS2mS0, pltS2mS1, pltS2mS2
    #=# for plotting S=0 and S=2 separately
    # plot S=0
    S0index = findall(x->x.S==0,unq)[1] # index of the S=0 ket
    pltS0=plot(ustrip.(uconvert.(u"T",pltdata[S0index][1])),austrip.(pltdata[S0index][2]),
    xlabel="B (T)", ylabel="σ (a₀²)", minorticks=true, label=pltlabel[S0index], legend=:outertopright)
    # plot S=2 states
    S2indices=filter(x->x!=S0index, 1:length(unq))
    pltS2=plot(ustrip.(uconvert.(u"T",pltdata[S2indices[1]][1])),austrip.(pltdata[S2indices[1]][2]),
    xlabel="B (T)", ylabel="σ (a₀²)", minorticks=true, label=pltlabel[S2indices[1]], legend=:outertopright)
    if length(S2indices)>1
        for i=2:length(S2indices)
            plot!(ustrip.(uconvert.(u"T",pltdata[S2indices[i]][1])),austrip.(pltdata[S2indices[i]][2]),
            label=pltlabel[S2indices[i]], legend=:outertopright)
        end
    end
    #hline!([4*pi*austrip((7.54u"nm")^2)],label="4π×(7.54nm)²") # S=2 theoretical σ
    # merge plots
    plot(pltS0, pltS2, layout=(2,1),title=["Ionisation, k=$k, lmax=$lmax" ""],
    left_margin=5mm,bottom_margin=5mm,top_margin=5mm,linewidth=2,grid=false)=#
end

########################Fano fitting############################################
using CurveFit
""" Returns least squares fit of Fano profile to elastic cross section data,
    given a state γ, wavenumber k, and lmax (along with Bmin, Bmax fitting region)"""
function nonlin_fano(γ::γ_ket,Bmin::Unitful.BField,Bmax::Unitful.BField,
    k::Union{typeof(0u"bohr^-1"),typeof(0e0u"bohr^-1")},lmax::Int)
    # load all data with correct B
    datas=load_data("gam",-(Inf)u"hartree",(Inf)u"hartree",Bmin,Bmax,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.B)) # sort by increasing Bfield
    Bs::Array{Unitful.BField,1}=[]
    σs::Array{typeof(0e0u"bohr^2"),1}=[]
    for d in datas # already sorted datas by energy
        if γ in d.γ_lookup # in case γ is closed in this data
            dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
            # need match of energy to B field for this γ
            d.ϵ==E∞(γ,k,d.B) || continue
            push!(Bs,d.B) # store Bfield
            push!(σs,d.σ[dγ_index,dγ_index]) # store elastic cross section
        end
    end
    #p[1≡q, p[2]≡Δ, p[3]≡B₀
    fano(x,p) = let
        ϵ=2*(x-p[3])/p[2]
        (ϵ+p[1])^2/(ϵ^2+1) # Fano eqn from wolfram demo
    end
    fano_err(xy, p) = xy[2] - fano(xy[1],p)
    # strip units before passing to fit
    Bs⁰=ustrip.(uconvert.(u"T",Bs)*10_000) # Bfields in G so Δ, B₀ in G
    σs⁰=austrip.(σs) # σ in a₀²
    # normalise peak height
    σs⁰ ./= maximum(σs⁰)
    # move min to zero
    xy=[Bs⁰ σs⁰]
    # parameter starting points
    backgnd=(σs⁰[1]+σs⁰[end])/2
    q0=30.1; Δ0=1e-4; B0=127.93 # Best guess (B field stuff in G)
    p0=[q0,Δ0,B0,backgnd] # initial parameters choice
    #plot(xy[:,1],xy[:,2]); return plot!(xy[:,1],[fano(i,p0) for i in xy[:,1]])
    eps=1e-7; maxiter=2000;
    coefs,converged,iter = nonlinear_fit(xy,fano_err,p0,eps,maxiter)
end

# linear fit
function lin_fano(γ::γ_ket,Bmin::Unitful.BField,Bmax::Unitful.BField,
    k::Union{typeof(0u"bohr^-1"),typeof(0e0u"bohr^-1")},lmax::Int)
    # load all data with correct B
    datas=load_data("gam",-(Inf)u"hartree",(Inf)u"hartree",Bmin,Bmax,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.B)) # sort by increasing Bfield
    Bs::Array{Unitful.BField,1}=[]
    σs::Array{typeof(0e0u"bohr^2"),1}=[]
    for d in datas # already sorted datas by energy
        if γ in d.γ_lookup # in case γ is closed in this data
            dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
            # need match of energy to B field for this γ
            d.ϵ==E∞(γ,k,d.B) || continue
            push!(Bs,d.B) # store Bfield
            push!(σs,d.σ[dγ_index,dγ_index]) # store elastic cross section
        end
    end
    #p[1≡q, p[2]≡Δ, p[3]≡B₀
    fano(x,p) = let
        ϵ=2*(x-p[3])/p[2]
        (ϵ+p[1])^2/(ϵ^2+1) # Fano eqn from wolfram demo
    end
    @. fano(𝐱,p) = (p[1]*p[2]/2 + 𝐱 - p[3])^2/((p[2]/2)^2+(𝐱-p[3])^2) + p[4]
    # strip units before passing to fit
    Bs⁰=ustrip.(uconvert.(u"T",Bs)*10_000) # Bfields in G so Δ, B₀ in G
    σs⁰=austrip.(σs) # σ in a₀²
    # normalise peak height
    σs⁰ ./= maximum(σs⁰)
    # estimate p[4]=background
    backgnd=(σs⁰[1]+σs⁰[end])/2
    # parameter starting points
    q0=1.0; Δ0=0.0001; B0=127.93 # Best guess (B field stuff in G)
    p0=[q0,Δ0,B0,backgnd] # initial parameters choice
    #return xy, fano, p0
    fit=LsqFit.curve_fit(fano,Bs⁰,σs⁰,p0)
    plot(Bs⁰,σs⁰)
    plot!(Bs⁰,fano(Bs⁰,coef(fit)))
end

using DelimitedFiles
# script for saving unit-less Bfields, σₑₗ's for python fitting
function saveBsσs(γ::γ_ket,Bmin::Unitful.BField,Bmax::Unitful.BField,
    k::Union{typeof(0u"bohr^-1"),typeof(0e0u"bohr^-1")},lmax::Int)
    # load all data with correct B
    datas=load_data("I",-(Inf)u"hartree",(Inf)u"hartree",Bmin,Bmax,lmax)
    @assert length(datas)>0 "Didn't find any suitable data"
    sort!(datas, by=(x->x.B)) # sort by increasing Bfield
    Bs::Array{Unitful.BField,1}=[]
    σs::Array{typeof(0e0u"bohr^2"),1}=[]
    for d in datas # already sorted datas by energy
        if γ in d.γ_lookup # in case γ is closed in this data
            dγ_index=findall(x->x==γ,d.γ_lookup)[1] # order of γ in σ array
            # need match of energy to B field for this γ
            d.ϵ==E∞(γ,k,d.B) || continue
            push!(Bs,d.B) # store Bfield
            push!(σs,d.σ[dγ_index]) # store elastic cross section
        end
    end
    # convert Bs to Gauss and σs to a₀²
    Bs=uconvert.(u"T",Bs)*10_000
    σs=auconvert.(σs)
    # strip units and form N×2 array
    Bs⁰=ustrip.(Bs)
    σs⁰=ustrip.(σs)
    Bσ=[Bs⁰ σs⁰]
    println("Saving $(length(Bs⁰)) data for S=$(γ.S), mS=$(γ.mS)")
    prevdir=pwd()
    cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Sem 2\October\mS-12-ion-decay")
    writedlm("S$(γ.S)mS$(γ.mS)_Bmin$(ustrip(uconvert(u"T",Bmin)*10_000))G_Bmax$(ustrip(uconvert(u"T",Bmax)*10_000))_el_vs_B.csv", Bσ, ',')
    cd(prevdir)
    plot(Bs⁰,σs⁰)
end

######################produce potentials .csv for DC to check###################
using Potentials
function gridpotentials(S::Int)
    @assert S in [0,1,2] "S not 0, 1, or 2"
    Rs=LinRange(3,200,100_000)u"bohr"
    if S==0
        Vs=Singlet.(Rs)
    elseif S==1
        Vs=Triplet.(Rs)
    else
        Vs=Quintet.(Rs)
    end
    prevdir=pwd()
    cd(raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Sem 2\October\BOpots-values")
    writedlm("S$S.csv", [austrip.(Rs) austrip.(Vs)], ',')
    cd(prevdir)
    plot(austrip.(Rs), austrip.(Vs))
end
