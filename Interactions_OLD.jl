# Old code, deals with |a⟩ kets from Cocks et al (2019) and not the |SmSlml⟩
# kets from Beams et al (2006).

#= Interactions as functions of |a⟩ kets, separation distance R and mass μ.=#

module Interactions
export H_rot, H_el, H_sd_coeffs, H_sd_radial

using Traceur

using HalfIntegers, Unitful, UnitfulAtomic

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

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

#TODO hyperfine, zeeman

end #module
