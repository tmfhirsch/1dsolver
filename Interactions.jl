#= Interactions as functions of |a⟩ kets, separation distance R and mass μ.=#

#module Interactions

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")
using HalfIntegers, Unitful, UnitfulAtomic

# rotational
using StateStructures: a_ket
include("./interactions/H_rot.jl")

 # electronic
using StateStructures: a12_ket, a_ket, asymmetric_eval
using Potentials
using Wigner9j
include("./interactions/H_el.jl")

 # spin-dipole
using StateStructures: a12_ket, α_ket, a_ket, asymmetric_eval
using WignerSymbols
include("./interactions/H_sd.jl")

#end #module
