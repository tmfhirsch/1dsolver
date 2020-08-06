#= Comparing my H_sd coupling scheme with the coupling scheme of Cocks
et al. (2019)=#
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using WignerSymbols, Wigner9j

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

"""Redefined CG, returns 0 for unphysical (jᵢ,mᵢ) combination (instead of error)"""
function clebschgordan_lax(j₁,m₁,j₂,m₂,j₃,m₃=m₁+m₂)
    for (jᵢ,mᵢ) in ((j₁, m₁), (j₂, m₂), (j₃, m₃))
        if !WignerSymbols.ϵ(jᵢ,mᵢ) # unphysical jᵢ,mᵢ entered
            return 0::Int
        end
    end
    return clebschgordan(j₁,m₁,j₂,m₂,j₃,m₃) # if all good, send on to CG calc
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
            Cterm *= clebschgordan_lax(l_,0,C,0,l,0)
            Cterm *= clebschgordan_lax(l_,ml_,C,mS-mS_,l,ml)
            Cterm *= clebschgordan_lax(S_,mS_,C,mS-mS_,S,mS)
            Cterm *= clebschgordan_lax(S_,Ωₛ,C,0,S,Ωₛ)
            Csum += Cterm
            Cterm = 0
        end
        term = Ω₁*Ω₂
        term *= clebschgordan_lax(S₁,Ω₁,S₂,Ω₂,S_,Ωₛ)
        term *= clebschgordan_lax(S₁,Ω₁,S₂,Ω₂,S,Ωₛ)
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

#=
# Test that my scheme runs
function test_my_scheme()
    testα_=α(1,1,1,1,0,0)
    testα=α(1,1,2,0,1,1)
    my_scheme(testα_,testα)
end
=#


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
    eval=(-1)^(mS_-mS) #(-1)^(mₛ₋-mₛ) factor
    eval *= clebschgordan(S,mS,2,mS_-mS,S_,mS_)
    eval *= clebschgordan(l,ml,2,ml_-ml,l_,ml_)
    eval *= TTensor((S₁_,S₂_),S_, (S₁,S₂),S)
    eval *= sqrt((2*l+1)/(2*l_+1))*clebschgordan(l,0,2,0,l_,0)
    return sqrt(6)*eval
end

""" ⟨γ',S'|𝐓²|γ,S⟩/ħ² from (36) in Cocks et al (2019)
    Inputs: γ'={S1',S2'}, S', γ={S1,S2}, S
    Outputs: ⟨γ',S'|𝐓²|γ,S⟩/ħ²"""
function TTensor(γ_,S_,γ,S)
    S1_, S2_ = γ_[1], γ_[2]
    S1, S2 = γ[1], γ[2]
    if γ_ != γ #δᵧ_ᵧ term
        return 0
    end
    eval = sqrt(5*S1*(S1+1)*S2*(S2+1))
    eval*= sqrt((2*S1+1)*(2*S2+1)*(2*S+1))
    eval*= wigner9j(S1,S2,S,1,1,2,S1,S2,S_)
    return eval
end

# test my scheme against Cocks 2019
function tester(tol=1e-10)
    α1=α(1,1,2,-1,2,0)
    α2=α(1,1,1,0,0,0)
    α3=α(1,1,0,0,1,-1)
    α4=α(1,1,2,1,2,1)
    αs=[α1,α2,α3,α4]
    for α_ in αs, α in αs
        println("α'=$α_, α=$α")
        if abs(my_scheme(α_,α)-cocks2019_scheme(α_,α))>=tol
            @info α_, α, "Disagreement between my scheme and Cocks (2019)"
            @info "My scheme produces" my_scheme(α_,α)
            @info "Cocks (2019) scheme produces" cocks2019_scheme(α_,α)
            return
        end
    end
end

#playing around outside fn
α1=α(1,1,2,-1,2,0)
α2=α(1,1,1,0,0,0)
α3=α(1,1,0,0,1,-1)
α4=α(1,1,2,1,2,1)
