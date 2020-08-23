#= Module containing structs used to hold different states (w/ differerent
quantum numbers), and functions to convert between them. For example, the |a⟩
states and the |a₁₂⟩ states.=#
module StateStructures

using HalfIntegers

""" |a⟩=|Γ,f,l,J,mJ,XS⟩ structure"""
struct α_nos
    S::HalfInteger
    i::HalfInteger
    f::HalfInteger
end
struct Γ_nos; α₁::α_nos; α₂::α_nos; end
struct a_ket # symmetrised states
    Γ::Γ_nos
    f::HalfInteger
    l::Integer
    J::HalfInteger
    mJ::HalfInteger
    XN::HalfInteger
end


"""|a₁₂⟩ unsymmetrised version of |a⟩; used to evaluate  ̂Hₑₗ and ̂H_sd"""
struct a12_ket
    Γ::Γ_nos
    f::HalfInteger
    l::Integer
    J::HalfInteger
    mJ::HalfInteger
end
"""Generates unsymmetrised state |a₁₂⟩
    Input: |a⟩ state containing Γ={α₁,α₂}
    Output: related |a₁₂⟩"""
function a12_maker(a::a_ket)
    return a12_ket(a.Γ,a.f,a.l,a.J,a.mJ)
end
"""Generates related unsymmetrised state |a₂₁⟩. Does NOT do the proper equality
with phase factors.
    Input: |a⟩ state containing Γ={α₁,α₂}
    Output: related |a₂₁⟩"""
function a21_maker(a::a_ket)
    Γ=Γ_nos(a.Γ.α₂, a.Γ.α₁) # flips α₁ and α₂
    return a12_ket(Γ,a.f,a.l,a.J,a.mJ)
end

"""Given |a⟩-type states, and a function representing the interaction that
operates on |a₁₂⟩-type states, calculates the equivalent evaluation using
unsymmetrised states. *Assumes mass is irrelevant*
    Input: bra::|a⟩, ket::|a⟩, H_asymmetric ::|a₁₂⟩ × |a₁₂⟩ × R → [E]
    Output: ⟨bra|̂H|ket⟩, evaluated by expanding into unsymmetrised basis"""
function asymmetric_eval(H12, bra::a_ket, ket::a_ket, R)
    # 1/sqrt(2*(1+δ(α₁α₂))) prefactor
    prefac_ = bra.Γ.α₁==bra.Γ.α₂ ? 1/2 : 1/sqrt(2)
    prefac  = ket.Γ.α₁==ket.Γ.α₂ ? 1/2 : 1/sqrt(2)
    # (-1)^(Xₙ+l+f₁+f₂-f) phase factors
    phase_=(-1)^(bra.XN + bra.l + bra.Γ.α₁.f + bra.Γ.α₂.f - bra.f)
    phase =(-1)^(ket.XN + ket.l + ket.Γ.α₁.f + ket.Γ.α₂.f - ket.f)
    # construct asymmetric states
    a12_, a21_ = a12_maker(bra), a21_maker(bra)
    a12,  a21  = a12_maker(ket), a21_maker(ket)
    # return expansion in terms of unsymmetrised states # see notebook 4 14/8/20
    return prefac_*prefac*(H12(a12_,a12,R)
                           +phase_*H12(a21_,a12,R)
                           +phase*H12(a12_,a21,R)
                           +phase_*phase*H12(a21_,a21,R)
                           )
end

"""Generates all |a⟩ states up to and including lmax
    Input: lmax=4
    Output: vector of a_kets
    Tested 12/08/20, I believe it is working correctly"""
function a_lookup_generator(lmax)
    lookup=Vector{a_ket}()
    for S₁ in HalfInt.(1:1), i₁ in HalfInt.(0:0), S₂ in HalfInt.(1:1), i₂ in HalfInt.(0:0)
        for f₁ in abs(S₁-i₁):1:(S₁+i₁), f₂ in abs(S₂-i₂):1:(S₂+i₂)
            α₁, α₂ = α_nos(S₁,i₁,f₁), α_nos(S₂,i₂,f₂)
            Γ = Γ_nos(α₁,α₂)
            for f in abs(f₁-f₂):1:(f₁+f₂), l in HalfInt.(0:1:lmax)
                for J in abs(f-l):1:(f+l)
                    for mJ in (-J):1:J
                        @assert i₁==i₂ "i₁!=i₂, symmetrised states not suitable"
                        XN = i₁== 0 ? 0 : 1 # bosonic or fermionic
                        if α₁==α₂ # check that state exists
                            XN = i₁== 0 ? 0 : 1 # bosonic or fermionic
                            phase=(-1)^(XN+l+f₁+f₂-f)
                            phase == 1 || continue # skip if state doesn't exist
                        end
                        push!(lookup, a_ket(Γ,f,l,J,mJ,XN))
                    end
                end
            end
        end
    end
    return lookup
end

################################################################################
# |α⟩ kets
"""|α⟩=|S₁ S₂ S mₛ⟩|l mₗ⟩ states. Used in evaluation of ̂H_sd"""
struct α_ket
    S₁::HalfInteger
    S₂::HalfInteger
    S::HalfInteger
    mS::HalfInteger
    l::Integer
    ml::Integer
end


end #module
