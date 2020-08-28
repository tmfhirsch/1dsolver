# Old code, deals with |a⟩ kets from Cocks et al (2019) and not the |SmSlml⟩
# kets from Beams et al (2006).
"""Rotational interaction
    Input: ⟨a'|, |a⟩, R~[L], μ~[M]
    Output: ⟨a'|̂H_rot|a⟩(R) ~ [E]"""
function H_rot(bra::a_ket,ket::a_ket,R,μ)
    # delta function
    bra == ket || return 0.0u"hartree"
    l=ket.l
    return uconvert(u"hartree", l*(l+1)*1.0u"ħ^2"/(2*μ*R^2))
end
