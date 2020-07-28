#=
Singlet, triplet and quintet Born-Oppenheimer potentials ¹Σ⁺g, ³Σ⁺u, and ⁵Σ⁺g.
Singlet and triplet from Müller and Przybytek.

Description last updates 27/07/2020
=#
module Potentials

using Unitful, UnitfulAtomic, Dierckx
using Plots

################################################################################
# Przybytek Quintet Potential
################################################################################
"""
    Damping function of Tang and Toennies [31 in Przybytek]
    Input: n::Int index of gₙ (note diff. to g₂ₙ in paper); x
    Output: gₙ(x) (Eₕ)
"""
function damping(n::Int, x)
    taylor = sum([x^i/(factorial(i)) for i=0:n])
    return (1-exp(-x)*taylor)u"hartree"
end

"""
    Przybytek analytic expression for the ⁵Σ⁺g He₂ potential, in atomic units
    Input: radial distance R (a₀)
    Output: potential V(R) (Eₕ)
"""
function Quintet(R)
    # Fitted constants from Table III. Here, Xᵈ ≡ X' as seen in the paper
    # constants all checked cf Table III 13/07/20
    A   = 1.2378383u"hartree" # [E]
    α   = 4.9176716e-1u"bohr^-1" # [L]⁻¹
    β   = 4.1575266e-2u"bohr^-2" # [L]⁻²
    Aᵈ  = 1.7014074e-2u"hartree" # [E]
    βᵈ  = 3.2614160e-2u"bohr^-2" # [L]⁻²
    η   = 7.6464764e-1u"bohr^-1" # [L]⁻¹
    C₆  = 3.276680e3u"bohr^6" # [L]⁶
    C₈  = 2.1056655e5u"bohr^8" # [L]⁸
    C₁₀ = 2.1786760e7u"bohr^10" # [L]¹⁰
    C₁₂ = (3.137e9 + 1.503e8)u"bohr^12" # [L]¹²
    C₁₁ = -1.671e8u"bohr^11" # [L]¹¹

    V = (A*exp(-α*R-β*R^2)
        +Aᵈ*exp(-βᵈ*R^2)
        -damping(6,η*R)*C₆/R^6
        -damping(8,η*R)*C₈/R^8
        -damping(10,η*R)*C₁₀/R^10
        -damping(12,η*R)*C₁₂/R^12
        -damping(11,η*R)*C₁₁/R^11)
    return V
end

"""
    Converts hartrees to cm⁻¹ (for use in test plot of Przybytek potential)
    Input: E in Eh
    Output: Corresponding ν in cm⁻¹ where E=hcν
"""
hartree2wavenumber(E)=uconvert(u"cm^-1", E/u"h*c")

#=#Test plot of przybytek potential
Rgrid = LinRange(5,26,10001)u"bohr"
Vgrid = Quintet.(Rgrid) # Eₕ
νgrid = hartree2wavenumber.(Vgrid) # cm⁻¹

plot(ustrip.(Rgrid), ustrip.(νgrid),
     xlabel="R (a₀)", ylabel="V(R) (cm⁻¹)",
     legend=false
     )=#

################################################################################
# Müller Singlet potential
################################################################################
"""
    Converts hartrees to eV (for test plots of Müller potentials)
    Input: E in Eₕ (really, any energy unit)
    Output: E in eV
"""
hartree2eV(E) = uconvert(u"eV", E)

"""
    Müller Singlet potential, interpolated with cubic spline (DC recommendation)
    and fitted with a decaying exponential to the Przybytek potential for R>14a₀.
    Input: Radial distance R (distance units, pref. a₀)
    Output: ¹Σ⁺g(R) (Eₕ)
"""
function Singlet(R)
    # Cubic spline interpolation of tabulated values
    interp_Rs = push!([i for i in 3:0.5:8],9,10,11,12,14) # interpolation R vals
    interp_Vs = [1877.8, # Values from Müller Table 3
                 1159.3, # Numbers checked by comparison with the table 27/07/20
                 461.1,
                 -78.9,
                 -422.4,
                 -602.4,
                 -665.1,
                 -652.5,
                 -595.7,
                 -517.1,
                 -431.6,
                 -275.0,
                 -161.7,
                 -91.1,
                 -50.4,
                 -15.8]
    interp_pot_unitless = Spline1D(interp_Rs, interp_Vs, k=3)
    interp_pot(R) = uconvert.(u"hartree",interp_pot_unitless(austrip(R))*1.0u"meV")
    # Exponential decay exchange function
    A₁, β₁ = 5.9784u"hartree", 0.7367u"bohr^-1" # values from Cocks 2019
    exch_pot(R) = A₁*exp(-β₁*R)
    if R < 3.0u"bohr" # within min range of Müller values
        error("R=$R, within radius of Müller tabulated values")
    elseif R < 14.0u"bohr" # within tabulated value range
        return interp_pot(R)
    else # asymptotic fit to quintet potential
        return Quintet(R) - exch_pot(R)
    end
end

# Test plot of Singlet potential, cf. Fig 6 and Table 4 of Müller
# Tested correct as of 27/07/20
#=Rs=LinRange(3,18,100)u"bohr"
vals=uconvert.(u"eV", Singlet.(Rs))
plot(ustrip.(Rs),ustrip.(vals))
vline!([6.15])
hline!([-0.6677])=#

################################################################################
# Müller Triplet potential
################################################################################

"""
    Müller Triplet potential, interpolated with cubic spline (DC recommendation)
    and fitted with a decaying exponential to the Przybytek potential for R>14a₀.
    Input: Radial distance R (distance units, pref. a₀)
    Output: ¹Σ⁺g(R) (Eₕ)
"""
function Triplet(R)
    # Cubic spline interpolation of tabulated values
    interp_Rs = push!([i for i in 3:0.5:8],9,10,11,12,14) # interpolation R vals
    interp_Vs = [2108.3, # from Table 3, checked 27/02/20
                 1100.4,
                 386.0,
                 -115.6,
                 -422.4,
                 -577.1,
                 -625.1,
                 -604.0,
                 -543.2,
                 -463.5,
                 -379.7,
                 -233.9,
                 -135.7,
                 -77.8,
                 -45.1,
                 -15.3]
    interp_pot_unitless = Spline1D(interp_Rs, interp_Vs, k=3)
    interp_pot(R) = uconvert.(u"hartree",interp_pot_unitless(austrip(R))*1.0u"meV")
    # Exponential decay exchange function
    A₃, β₃ = 1.7980u"hartree", 0.6578u"bohr^-1" # values from Cocks 2019
    exch_pot(R) = A₃*exp(-β₃*R)
    if R < 3.0u"bohr" # within min range of Müller values
        error("R=$R, within radius of Müller tabulated values")
    elseif R < 14.0u"bohr" # within tabulated value range
        return interp_pot(R)
    else # asymptotic fit to quintet potential
        return Quintet(R) - exch_pot(R)
    end
end

# Test plot of Triplet potential, cf. Fig 6 and Table 4 of Müller
# Tested correct as of 27/07/20
#=Rs=LinRange(3,18,100)u"bohr"
vals=uconvert.(u"eV", Triplet.(Rs))
plot(ustrip.(Rs),ustrip.(vals))
vline!([6.10])
hline!([-0.6256])=#

end # module
