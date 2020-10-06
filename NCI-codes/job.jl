# the Julia script to be called in submit.sh (the shell script that is queued)
const Smat_dir = "~/resultsA/Smat"
const gampwcs_dir = "~/resultsA/gampwcs"
const Ics_dir = "~/resultsA/Ics"

push!(LOAD_PATH,"~/Code/Modules")

include("./SimResults.jl")

# EDIT THE LINE BELOW BETWEEN RUNS
gen_diffB_constk_data(5458G,5462G,10000,1e-4u"bohr^-1",4)
