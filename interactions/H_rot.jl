"""Rotational interaction for |Φₐ⟩≡|SmS⟩|lml⟩ kets.
    Input: ⟨Φₐ'|, |Φₐ⟩, R~[L], μ~[M]
    Output: ⟨Φₐ'|̂H_rot|Φₐ⟩(R) ~ [E]"""
function H_rot(bra::SmS_ket,ket::SmS_ket,R,μ)
    bra == ket || return 0.0u"hartree" # delta function
    l=ket.l
    return uconvert(u"hartree", l*(l+1)*1.0u"ħ^2"/(2*μ*R^2))
end
