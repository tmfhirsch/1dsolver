#=
Setting up Przybytek potential and getting to grips with Unitful

Description last updated 13/07/20
=#
using Unitful, UnitfulAtomic
using Plots
using Revise

"""
Damping function of Tang and Toennies [31 in Przybytek]
    Input: n::Int index of gₙ (note diff. to g₂ₙ in paper); x
    Output: gₙ(x)
"""
function damping(n::Int, x)
    taylor = sum([x^i/(factorial(i)) for i=0:n])
    return (1-exp(-x)*taylor)u"hartree"
end


"""
Converts hartrees to cm⁻¹
    Input: E in Eh
    Output: Corresponding ν in cm⁻¹ where E=hcν
"""
hartree2wavenumber(E)=uconvert(u"cm^-1", E/u"h*c")


"""
Przybytek analytic expression for the ⁵Σ⁺g He₂ potential, in atomic units
    Input: radial distance R (a₀)
    Output: potential V(R) (Eₕ)
"""
function przybytek(R)
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

Rgrid = LinRange(6.8,8.2,10001)u"bohr"
Vgrid = przybytek.(Rgrid) # Eₕ
νgrid = hartree2wavenumber.(Vgrid) # cm⁻¹

plot(ustrip.(Rgrid), ustrip.(νgrid),
     xlabel="R (a₀)", ylabel="V(R) (cm⁻¹)",
     legend=false,
     yticks=(-940):-20:(-1040),
     xticks=6.8:0.2:8.2
     )
