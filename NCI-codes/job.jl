#Bmin, Bmax, number of sims arguments
@assert length(ARGS)==3 "Not given 3 command-line arguments"
Bmin, Bmax, num = parse(Float64,ARGS[1]), parse(Float64,ARGS[2]), parse(Int,ARGS[3])

# the Julia script to be called in submit.sh (the shell script that is queued)
const Smat_dir = "/" * joinpath("home","111","th8512","resultsA","Smat")
const gampwcs_dir = "/" * joinpath("home","111","th8512","resultsA","gampwcs")
const Ics_dir =  "/" * joinpath("home","111","th8512","resultsA","Ics")

modDir = "/" * joinpath("home","111","th8512","Code","Modules")
union!(LOAD_PATH, [".",modDir])

include("./SimResults.jl")

# generate data
# EDIT THE LINE BELOW BETWEEN RUNS
gen_diffB_constk_data(Bmin*G,Bmax*G,num,1e-4u"bohr^-1",4)

# convert/data
S2γ_data()
S2I_data()
