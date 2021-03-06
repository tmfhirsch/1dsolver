#= Interactions as functions of |Φₐ⟩=|S₁S₂SmSlml⟩ states, radial dist. R and μ=#

module Interactions
export H_rot, H_el, H_sd_coeffs, H_sd_radial, H_zee

using Unitful, UnitfulAtomic

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using StateStructures

# rotational
include("./interactions/H_rot.jl")

# electronic
using Potentials
include("./interactions/H_el.jl")

# spin-dipole
using WignerSymbols, Wigner9j
include("./interactions/H_sd.jl")

# zeeman
include("./interactions/H_zee.jl")

end # module
