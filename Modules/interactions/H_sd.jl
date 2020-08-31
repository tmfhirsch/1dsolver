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

""" Radial factor in H_sd, Vₚ(R)
    Input: R~[L]
    Output: Vₚ(R)~[E]
    Based off Cocks et al (2019) eqn (39)=Beams et al (2006) below (9)"""
function H_sd_radial(R)
    # radial factor
    b=-sqrt(6)*ξ # ħ⁻² not written since it cancels with the TTensor factor
    Vₚ=b/R^3
    return Vₚ
end

"""Coefficients for spin-dipole interaction for |Φₐ⟩≡|SmS⟩|lml⟩ kets.
    Input: ⟨Φₐ'|, |Φₐ⟩, R~[L],
    Output: ⟨Φₐ'|̂H|Φₐ⟩/Vₚ(R) ~ 1"""
function H_sd_coeffs(bra::SmS_ket,ket::SmS_ket)
    # unpack quantum numbers
    S₁_, S₂_, S_, mS_, l_, ml_ = bra.Γ.S₁, bra.Γ.S₂, bra.S, bra.mS, bra.l, bra.ml
    S₁, S₂, S, mS, l, ml = ket.Γ.S₁, ket.Γ.S₂, ket.S, ket.mS, ket.l, ket.ml
    # δ function
    (mS+ml)==(mS_+ml_) || return 0
    term=(-1)^(mS_-mS)
    term*= clebschgordan_lax(S,mS,2,mS_-mS,S_,mS_)
    term*=clebschgordan_lax(l,ml,2,ml_-ml,l_,ml_)
    term*=TTensor(S₁,S₂,S,S₁_,S₂_,S_)
    term*=CTensor(l,l_)
    term
end

"""⟨Γ'S'||𝐓²||ΓS⟩ from Beams et al. (2006) eqn (19)"""
function TTensor(S₁,S₂,S,S₁_,S₂_,S_)
    (S₁==S₁_ && S₂==S₂_) || return 0
    𝐓² = sqrt(S₁*(S₁+1)*S₂*(S₂+1))
    𝐓²*= sqrt(5*(2*S₁+1)*(2*S₂+1)*(2*S+1))
    𝐓²*= wigner9j(S₁,S₂,S,1,1,2,S₁,S₂,S_)
    𝐓²
end

"""⟨Γ'S'||𝐂²||ΓS⟩ from Beams et al. (2006) eqn (20)"""
CTensor(l,l_)=sqrt((2*l+1)/(2*l_+1))*clebschgordan_lax(l,0,2,0,l_,0)

##############################|a⟩ kets ######################################
"""̂H_sd Coupling coefficients for *symmetrised* |a⟩≡|ΓflJmJ⟩ states
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
    result*=CTensor(l,l)
    result*=wigner6j(f,2,f_,l_,J,l)
    Si_sum=0
    for S_ in abs(S₁_-S₂_):1:(S₁_+S₂_), S in abs(S₁-S₂):1:(S₁+S₂), i in abs(i₁-i₂):1:(i₁+i₂)
        Si_term=(-1)^(-S_-i)
        Si_term*=(2*S_+1)*(2*i+1)
        Si_term*=sqrt(2*S+1)
        Si_term*=wigner6j(f,2,f_,S_,i,S)
        Si_term*=TTensor(S₁_,S₂_,S_,S₁,S₂,S)
        Si_term*=wigner9j(S₁,S₂,S_,i₁,i₂,i,f₁_,f₂_,f_)
        Si_term*=wigner9j(S₁,S₂,S,i₁,i₂,i,f₁,f₂,f)
        Si_sum += Si_term
    end
    result*=Si_sum
    return result
end
