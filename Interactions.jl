#= Interactions as functions of |a⟩ kets, separation distance R and mass μ.=#

#module Interactions

push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")
using StateStructures

include(".\\interactions\\H_rot.jl") # rotational
include(".\\interactions\\H_el.jl") # electronic
include(".\\interactions\\H_sd.jl") # spin-dipole

#end #module
