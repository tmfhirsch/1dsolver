# Old code, deals with |a⟩ kets from Cocks et al (2019) and not the |SmSlml⟩
# kets from Beams et al (2006).
"""Asymmetrised states ̂Hₑₗ
    Input: ⟨a₁₂'|, |a₁₂⟩, R~[L]
    Output: ⟨a₁₂'|̂H_rot|a₁₂⟩(R) ~ [E]"""
function H_el_asym(bra::a12_ket,ket::a12_ket,R)
    # unpack quantum numbers
    S₁_, i₁_, f₁_ = bra.Γ.α₁.S, bra.Γ.α₁.i, bra.Γ.α₁.f
    S₁, i₁, f₁ = ket.Γ.α₁.S, ket.Γ.α₁.i, ket.Γ.α₁.f
    S₂_, i₂_, f₂_ = bra.Γ.α₂.S, bra.Γ.α₂.i, bra.Γ.α₂.f
    S₂, i₂, f₂ = ket.Γ.α₂.S, ket.Γ.α₂.i, ket.Γ.α₂.f
    f_, l_, J_, mJ_ = bra.f, bra.l, bra.J, bra.mJ
    f, l, J, mJ = ket.f, ket.l, ket.J, ket.mJ
    # Δ(S₁,S₂,i₁,i₂,f,l,J,mJ)
    (S₁_==S₁ && S₂_==S₂ && i₁_==i₁ && i₂_==i₂ && f_==f && l_==l && J_==J && mJ_==mJ) || return 0.0u"hartree"
    # build solution expression
    # 9-j coupling prefactors
    x = sqrt((2*f₁_+1)*(2*f₂_+1)*(2*f₁+1)*(2*f₂+1))
    Si_sum = 0.0u"hartree"
    for S in abs(S₁-S₂):1:(S₁+S₂), i in abs(i₁-i₂):1:(i₁+i₂)
        term=(2*S+1)*(2*i+1)
        term*=wigner9j(S₁,S₂,S,i₁,i₂,i,f₁_,f₂_,f)
        term*=wigner9j(S₁,S₂,S,i₁,i₂,i,f₁,f₂,f)
        # Born-Oppenheimer potential
        @assert S==0 || S==1 || S==2 "S !∈ {0,1,2}"
        if S==0
            term*=Singlet(R)
        elseif S==1
            term*=Triplet(R)
        elseif S==2
            term*=Quintet(R)
        end
        Si_sum+=term
    end
    x*= Si_sum
    return x
end

"""Symmetrised states ̂Hₑₗ
    Input: ⟨a'|, |a⟩, R~[L]
    Output: ⟨a'|̂H_el|a⟩(R) ~ [E]"""
function H_el(bra::a_ket,ket::a_ket,R)
    return asymmetric_eval(H_el_asym, bra, ket, R)
end
