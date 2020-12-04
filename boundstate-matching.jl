using Revise
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver\Modules")
using StateStructures, Interactions
using Unitful, UnitfulAtomic, LinearAlgebra

lookup=γ_lookup_generator()
mS2=lookup[end]

const G=1e-4u"T" # 1 gauss
Z=-H_zee(mS2,mS2,1G)/2/1G
#now H_zee = -Z*mₛ*B, B measured in G but actually able to take any units

const k=1e-4u"bohr^-1"
const μ=0.5*4.002602u"u"
KE=auconvert(1u"ħ^2"*k^2/(2*μ))

Bres(E::Unitful.Energy,mS::Int)=(KE-E)/(Z*mS)
