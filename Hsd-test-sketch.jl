#= Comparing my H_sd coupling scheme with the coupling scheme of Cocks
et al. (2019)=#
using WignerSymbols

"""Structure to hold quantum numbers of the form
    |S₁ S₂ S mS l ml ⟩
    Initialise by entering the numbers in that order"""
mutable struct α
    S₁ :: Int
    S₂ :: Int
    S :: Int
    mS :: Int
    l :: Int
    ml :: Int
end

"""My coupling scheme, current as of 4/8/20.
    Inputs: α',α
    Outputs: ⟨α'|̂Hsd|α'⟩×(-R^3/ξ)
"""
function my_scheme(α_::α, α::α)
    #unpack quantum numbers
    S₁_,S₂_,S_,mS_,l_,ml_=α_.S₁,α_.S₂,α_.S,α_.mS,α_.l,α_.ml
    S₁, S₂, S, mS, l, ml = α.S₁, α.S₂, α.S, α.mS, α.l, α.ml
    #Δₛ₁ₛ₂ factor
    if S₁!=S₁_ || S₁!=S₂_
        return false
    end
    # First, coupling term
    Ωsum=0
    for Ωₛ in max(-S,-S_):1:min(S,S_), Ω₁ in -S₁:1:S₁, Ω₂ in -S₂:1:S₂ # outer sum
        # evaluate ∑C innermost sum
        Csum=0
        for C in abs(S-S_):1:(S+S_) # inner sum
            Cterm = (2C+1)
            Cterm *= clebschgordan(l_,0,C,0,l,0)
            Cterm *= clebschgordan(l_,ml_,C,mS-mS_,l,ml)
            Cterm *= clebschgordan(S_,mS_,C,mS-mS_,S,mS)
            Cterm *= clebschgordan(S_,Ωₛ,C,0,S,Ωₛ)
            Csum += Cterm
            Cterm = 0
        end
        term = Ω₁*Ω₂
        term *= clebschgordan(S₁,Ω₁,S₂,Ω₂,S_,Ωₛ)
        term *= clebschgordan(S₁,Ω₁,S₂,Ω₂,S,Ωₛ)
        term *= Csum
        Ωsum += term
        term = 0
    end
    coupling_term=Ωsum*3/(2*S+1)*sqrt((2l_+1)/(2l+1))
    # Second, diagonal term
    if S != S_ || mS != mS_ || l != l_ || ml != ml_
        diag_term = 0
    else
        diag_term = -0.5*(S*(S+1)-S₁*(S₁+1)-S₂*(S₂+1))
    end
    return coupling_term + diag_term
end


# Test that my scheme runs
function test_my_scheme()
    testα_=α(1,1,1,1,0,0)
    testα=α(1,1,2,0,1,1)
    my_scheme(testα_,testα)
end

""" Coupling scheme from Cocks et al. (2019)
    Inputs: (α', α) quantum numbers in that order
    Outputs: √6/ħ^2 *Dₐₐ_ as defined in (35) from the paper
"""
function cocks2019_scheme(α_::α, α::α)
    #unpack quantum numbers
    S₁_,S₂_,S_,mS_,l_,ml_=α_.S₁,α_.S₂,α_.S,α_.mS,α_.l,α_.ml
    S₁, S₂, S, mS, l, ml = α.S₁, α.S₂, α.S, α.mS, α.l, α.ml
    #δ(mₛ_+mₗ_,mₛ+mₗ) factor
    if mS_+ml_ != mS+ml
        return 0
    end
    term=(-1)^(mS_-mS) #(-1)^(mₛ₋-mₛ) factor
    term *= clebschgordan(S,mS,2,mS_-mS,S_,mS_)
    term *= clebschgordan(l,ml,2,ml_-ml,l_,ml_)
    #TODO tensor factors
end

"""Wigner 9j calculator, based off expression (1998) doi: 10.1063/1.168745"""
function wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    #TODO
end
