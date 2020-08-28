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
