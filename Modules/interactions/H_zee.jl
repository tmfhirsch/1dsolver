const gₛ = 2.00231930436256 # Electron spin g-facctor

""" Zeeman Interaction ̂H_zee=gₛμ_B/ħ*B*S_z, from Venturi & Whittingham 1999)
    Inputs: bra, ket, B~[Tesla] = magnetic field strength
    Outputs: ̂⟨bra|H_zee|ket⟩ ~ [E], default units of Eₕ """
function H_zee(bra::SmS_ket,ket::SmS_ket,B)
    @assert dimension(B)==dimension(1u"T") "B not in units of magnetic field"
    # diagonal in SmS kets
    bra == ket || return 0.0u"hartree"
    S_z = ket.mS
    return 0.5u"e_au/me_au"*gₛ*B*S_z
end
