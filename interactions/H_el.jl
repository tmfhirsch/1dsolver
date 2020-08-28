"""Electronic interaction for |Φₐ⟩≡|SmS⟩|lml⟩ kets.
    Input: ⟨Φₐ'|, |Φₐ⟩, R~[L]
    Output: ⟨Φₐ'|̂H_el|Φₐ⟩(R) ~ [E]"""
function H_el(bra::SmS_ket,ket::SmS_ket,R)
    bra == ket || return 0.0u"hartree" # δ function; the |Φₐ⟩ basis diagonalises
    S=ket.S
    @assert S==0 || S==1 || S==2 "S !∈ {0,1,2}"
    if S==0
        return Singlet(R)
    elseif S==1
        return Triplet(R)
    elseif S==2
        return Quintet(R)
    end
end
