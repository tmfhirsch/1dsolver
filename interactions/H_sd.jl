"""Redefined CG, returns 0 for unphysical (jᵢ,mᵢ) combination (instead of error)"""
function clebschgordan_lax(j₁,m₁,j₂,m₂,j₃,m₃=m₁+m₂)
    for (jᵢ,mᵢ) in ((j₁, m₁), (j₂, m₂), (j₃, m₃))
        if !WignerSymbols.ϵ(jᵢ,mᵢ) # unphysical jᵢ,mᵢ entered
            return 0::Int
        end
    end
    return clebschgordan(j₁,m₁,j₂,m₂,j₃,m₃) # if all good, send on to CG calc
end

# ξ factor (Cocks et al (2019) equation (7))
const μ_ratio = 1.00115965 # ratio of e⁻ mag mom to Bohr magneton (Cocks 2019)
const α_fs = 0.0072973525693 # fine structure constant (Wikipedia)
const ξ = α_fs^2 * μ_ratio^2 * 1.0u"hartree*bohr^3"

################################################################################
# My code - super slow w/ huge memory allocations.

"""Asymmetrised |a₁₂⟩ states ̂H_sd
    Input: ⟨a₁₂'|, |a₁₂⟩, R~[L]
    Output: ⟨a₁₂'|̂H_sd|a₁₂⟩(R) ~ [E]
    Based off eqn at end of 2.4.3 in my Thesis (23/08/20)"""
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
    result=sqrt((2*f₁_+1)*(2*f₂_+1)*(2*f₁+1)*(2*f₂+1))
    mfml_sum = 0
    for mf_ in -f_:1:f_, ml_ in -l_:1:l_, mf in -f:1:f, ml in -l:1:l
        mfml_term = clebschgordan_lax(f_,mf_,l_,ml_,J_,mJ_)
        mfml_term*= clebschgordan_lax(f,mf,l,ml,J,mJ)
        Si_sum = 0
        for S_ in abs(S₁_-S₂_):1:(S₁_+S₂_), S in abs(S₁-S₂):1:(S₁+S₂), i in abs(i₁-i₂):1:(i₁+i₂)
            Si_term = sqrt((2*S_+1)*(2*S+1))
            Si_term*= 2*i+1
            Si_term*= wigner9j(S₁_,S₂_,S_,i₁,i₂,i,f₁_,f₂_,f_)
            Si_term*= wigner9j(S₁,S₂,S,i₁,i₂,i,f₁,f₂,f)
            mSmi_sum = 0
            for mS_ in -S_:1:S_, mS in -S:1:S, mi in -i:1:i
                mSmi_term = clebschgordan_lax(S_,mS_,i,mi,f_,mf_)
                mSmi_term = clebschgordan_lax(S,mS,i,mi,f,mf)
                bra_α = α_ket(S₁_,S₂_,S_,mS_,l_,ml_)
                ket_α = α_ket(S₁,S₂,S,mS,l,ml)
                mSmi_term*= H_sd_α(bra_α,ket_α,R)
                mSmi_sum += mSmi_term
            end
            Si_term *= mSmi_sum
            Si_sum += Si_term
        end
        mfml_term *= Si_sum
        mfml_sum += mfml_term
    end
    result *= mfml_sum
    return (-15*ξ/R^3)*result
end

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
    return (1/(2*S_+1))*sqrt((2l_+1)/(2l+1))*Ωsum
end

################################################################################
# Scheme from Cocks et al. (2019). Coded up by me. Simpler than my monstrosity.

""" ⟨γ',S'|𝐓²|γ,S⟩/ħ² from (36) in Cocks et al (2019)
    Inputs: γ'={S1',S2'}, S', γ={S1,S2}, S
    Outputs: ⟨γ',S'|𝐓²|γ,S⟩/ħ²
    Tested for γ=(1,1) against ¶ below (37) in Cocks (2019) 10/08/20"""
function TTensor(γ_,S_,γ,S)
    S1_, S2_ = γ_[1], γ_[2]
    S1, S2 = γ[1], γ[2]
    γ_ == γ || return 0 #δᵧ_ᵧ
    x = sqrt(5*S1*(S1+1)*S2*(S2+1))
    x*= sqrt((2*S1+1)*(2*S2+1)*(2*S+1))
    x*= wigner9j(S1,S2,S,1,1,2,S1,S2,S_)
    return x
end

"""Coupling coefficients for *symmetrised* states ̂H_sd
    Input: ⟨a'|, |a⟩, R~[L]
    Output: ⟨a'|̂H_sd|a⟩(R)/Vₚ(R)"""
function H_sd_coeffs(bra::a_ket,ket::a_ket)
    return asymmetric_eval(H_sd_a12_coeffs, bra, ket, nothing)
end

""" Coupling coefficients of the symmetrised |a₁₂⟩ states ̂H_sd.
    Input: ⟨a₁₂'|, |a₁₂⟩, R~[L]
    Output: ⟨a₁₂'|̂H_sd|a₁₂⟩/Vₚ(R) ~ 1
    Based off Cocks et al (2019) eqn (39)"""
function H_sd_a12_coeffs(bra::a12_ket, ket::a12_ket)
    # unpack quantum numbers
    S₁_, i₁_, f₁_ = bra.Γ.α₁.S, bra.Γ.α₁.i, bra.Γ.α₁.f
    S₁, i₁, f₁ = ket.Γ.α₁.S, ket.Γ.α₁.i, ket.Γ.α₁.f
    S₂_, i₂_, f₂_ = bra.Γ.α₂.S, bra.Γ.α₂.i, bra.Γ.α₂.f
    S₂, i₂, f₂ = ket.Γ.α₂.S, ket.Γ.α₂.i, ket.Γ.α₂.f
    f_, l_, J_, mJ_ = bra.f, bra.l, bra.J, bra.mJ
    f, l, J, mJ = ket.f, ket.l, ket.J, ket.mJ
    # Δλ
    (S₁_==S₁ && S₂_==S₂ && i₁_==i₁ && i₂_==i₂ && J_==J && mJ==mJ) || return 0
    # Dₐ₁₂₋ₐ₁₂
    result=(-1)^(l_+J)
    result*=sqrt((2*f₁_+1)*(2*f₂_+1)*(2*f_+1)
                 *(2*f₁+1)*(2*f₂+1)*(2*f+1)
                 *(2*l_+1))
    result*=sqrt((2*l+1)/(2*l_+1))*clebschgordan_lax(l,0,2,l_,0)
    result*=wigner6j(f,2,f_,l_,J,l)
    Si_sum=0
    for S_ in abs(S₁_-S₂_):1:(S₁_+S₂_), S in abs(S₁-S₂):1:(S₁+S₂), i in abs(i₁-i₂):1:(i₁+i₂)
        Si_term=(-1)^(-S_-i)
        Si_term*=(2*S_+1)*(2*i+1)
        Si_term*=sqrt(2*S+1)
        Si_term*=wigner6j(f,2,f_,S_,i,S)
        Si_term*=TTensor((S₁_,S₂_),S_, (S₁,S₂),S)
        Si_term*=wigner9j(S₁,S₂,S_,i₁,i₂,i,f₁_,f₂_,f_)
        Si_term*=wigner9j(S₁,S₂,S,i₁,i₂,i,f₁,f₂,f)
        Si_sum += Si_term
    end
    result*=Si_sum
    return result
end

""" Radial factor in H_sd, Vₚ(R)
    Input: R~[L]
    Output: Vₚ(R)~[E]
    Based off Cocks et al (2019) eqn (39)"""
function H_sd_radial(R)
    # radial factor
    b=-sqrt(6)*ξ # ħ⁻² not written since it cancels with the TTensor factor
    Vₚ=b/R^3
    return Vₚ
end
