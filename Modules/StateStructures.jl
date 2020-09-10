module StateStructures

#########= |S₁S₂Smₛ⟩|lmₗ⟩  kets, from Beams et al. (2006)=########################
export SmS_ket, SmS_lookup_generator, γ_ket, γ_ket_convert

# S₁S₂ quantum numbers (=1,1 for He*)
struct SmS_Γ_nos # not using halfintegers since working in with metastable helium
    S₁ :: Int
    S₂ :: Int
end

# γ=|S₁S₂Smₛ⟩ kets, used to express cross sections after summming over |lmₗ⟩
struct γ_ket
    Γ :: SmS_Γ_nos
    S :: Int
    mS :: Int
end

# |Smₛ⟩=|S₁S₂Smₛlmₗ⟩ kets, used as channels for calculating wavefunctions
struct SmS_ket
    Γ :: SmS_Γ_nos
    S :: Int
    mS :: Int
    l :: Int
    ml :: Int
end

"""Generates all |...SmS⟩|lml⟩ states up to and including l=lmax"""
function SmS_lookup_generator(lmax)
    lookup = Vector{SmS_ket}()
    for S₁ in 1:1, S₂ in 1:1
        Γ=SmS_Γ_nos(S₁,S₂)
        for S in abs(S₁-S₂):1:(S₁+S₂), l in 0:1:lmax
            mod(l+S,2)==0 || continue # eigenvalue of XN=(-1)^(l+S)==1 (Bosonic)
            for mS in -S:1:S, ml in -l:1:l
                push!(lookup, SmS_ket(Γ,S,mS,l,ml))
            end
        end
    end
    lookup
end

# convert SmS_ket to γ_ket by stripping l,mₗ quantum numbers
function γ_ket_convert(ket::SmS_ket)
    S1S2=ket.Γ; S=ket.S; mS=ket.mS
    γ_ket(S1S2,S,mS)
end

##########################|a⟩≡|ΓflJmJ⟩ kets#####################################
export a_ket, a12_ket, a_ket, asymmetric_eval,  α_ket, a_lookup_generator
using HalfIntegers

""" |a⟩=|Γ,f,l,J,mJ,XS⟩ structure"""
struct α_nos
    S::HalfInteger
    i::HalfInteger
    f::HalfInteger
end
struct a_Γ_nos; α₁::α_nos; α₂::α_nos; end
struct a_ket # symmetrised states
    Γ::a_Γ_nos
    f::HalfInteger
    l::Integer
    J::HalfInteger
    mJ::HalfInteger
    XN::HalfInteger
end


"""|a₁₂⟩ unsymmetrised version of |a⟩; used to evaluate  ̂Hₑₗ and ̂H_sd"""
struct a12_ket
    Γ::a_Γ_nos
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
    Γ=a_Γ_nos(a.Γ.α₂, a.Γ.α₁) # flips α₁ and α₂
    return a12_ket(Γ,a.f,a.l,a.J,a.mJ)
end

"""Given |a⟩-type states, and a function representing the interaction that
    operates on |a₁₂⟩-type states, calculates the equivalent evaluation using
    unsymmetrised states. *Assumes mass is irrelevant*
    Input: bra::|a⟩, ket::|a⟩, ̂O12 ::|a₁₂⟩ × |a₁₂⟩ × R → [E]
    Output: ⟨bra|̂̂̂̂O|ket⟩ evaluated by expanding into unsymmetrised basis"""
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
    # R=nothing case (used to evaluate H_sd coupling coefficients)
    if R==nothing
        return prefac_*prefac*(H12(a12_,a12)
                               +phase_*H12(a21_,a12)
                               +phase*H12(a12_,a21)
                               +phase_*phase*H12(a21_,a21)
                               )
    end
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
            Γ = a_Γ_nos(α₁,α₂)
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

end # module
