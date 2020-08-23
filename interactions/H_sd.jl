"""Asymmetrised |a₁₂⟩ states ̂H_sd
    Input: ⟨a₁₂'|, |a₁₂⟩, R~[L]
    Output: ⟨a₁₂'|̂H_sd|a₁₂⟩(R) ~ [E]"""
function H_sd_asym(bra::a12_ket, ket::a12_ket, R)
    # unpack quantum numbers
    S₁_, i₁_, f₁_ = bra.Γ.α₁.S, bra.Γ.α₁.i, bra.Γ.α₁.f
    S₁, i₁, f₁ = ket.Γ.α₁.S, ket.Γ.α₁.i, ket.Γ.α₁.f
    S₂_, i₂_, f₂_ = bra.Γ.α₂.S, bra.Γ.α₂.i, bra.Γ.α₂.f
    S₂, i₂, f₂ = ket.Γ.α₂.S, ket.Γ.α₂.i, ket.Γ.α₂.f
    f_, l_, J_, mJ_ = bra.f, bra.l, bra.J, bra.mJ
    f, l, J, mJ = ket.f, ket.l, ket.J, ket.mJ
    # Δ(i₁,i₂)
    (i₁_==i₁ && i₂_==i₂) || return 0.0u"hartree"
    result=sqrt((2*f₁_+1)*(2*f₂+1)*(2*f₁+1)*(2*f₂+1))
    mfml_sum = 0
    for mf_ in -f_:1:f_, ml_ in -l_,1:l, mf in -f:1:f, ml in -l:1:l
        mfml_term = clebschgordan_lax(f_,mf_,l_,ml_,J_,mJ_)
        mfml_term*= clebschgordan_lax(f,mf,l,ml,J,mJ)
        Si_sum = 0
        #TODO up to here (2.72 in my current LyX document output)
end

# ξ factor (Cocks et al (2019) equation (7))
const μ_ratio = 1.00115965 # ratio of e⁻ mag mom to Bohr magneton (Cocks 2019)
const α_fs = 0.0072973525693 # fine structure constant (Wikipedia)
const ξ = α_fs^2 * μ_ratio^2 * 1.0u"hartree*bohr^3"

"""Asymmetrised |α⟩ states ̂H_sd
    Input: ⟨α'|, |α⟩, R~[L]
    Output: ⟨a₁₂'|̂H_sd|a₁₂⟩(R) ~ [E]"""
function H_sd_α(bra::α_ket, ket::α_ket, R)
    # unpack quantum numbers
    S₁_, S₂_, S_, mS_, l_, ml_ = bra.S₁, bra.S₂, bra.S, bra.mS, bra.l, bra.ml
    S₁, S₂, S, mS, l, ml = ket.S₁, ket.S₂, ket.S, ket.mS, ket.l, ket.ml
    #Δₛ₁ₛ₂ factor
    (S₁==S₁_ && S₂==S₂_) || return 0
    # Sum over rotated frame angular momenta
    Ωsum=0
    for Ωₛ in max(-S,-S_):1:min(S,S_), Ω₁ in -S₁:1:S₁, Ω₂ in -S₂:1:S₂
        term = clebschgordan_lax(l_,0,2,0,l,0)
        term *= clebschgordan_lax(l_,ml_,2,mS_-mS,l,ml)
        term *= clebschgordan_lax(S,mS,2,mS_-mS,S_,mS_)
        term *= clebschgordan_lax(S,Ωₛ,2,0,S_,Ωₛ)
        term *= Ω₁*Ω₂
        term *= clebschgordan_lax(S₁,Ω₁,S₂,Ω₂,S_,Ωₛ)
        term *= clebschgordan_lax(S₁,Ω₁,S₂,Ω₂,S,Ωₛ)
        Ωsum += term
    end
    return (-15*ξ/R^3)*(1/(2*S_+1))*sqrt((2l_+1)/(2l+1))*Ωsum
end

"""Symmetrised states ̂H_sd
    Input: ⟨a'|, |a⟩, R~[L]
    Output: ⟨a'|̂H_sd|a⟩(R) ~ [E]"""
function H_sd(bra::a_ket,ket::a_ket,R)
    return asymmetric_eval(H_sd_asym, bra, ket, R)
end
