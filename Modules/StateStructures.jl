#= |S₁S₂Smₛ⟩|lmₗ⟩  kets, from Beams et al. (2006)=#
module StateStructures
export SmS_ket, SmS_lookup_generator

struct Γ_nos # not using halfintegers since working in with metastable helium
    S₁ :: Integer
    S₂ :: Integer
end

struct SmS_ket
    Γ :: Γ_nos
    S :: Integer
    mS :: Integer
    l :: Integer
    ml :: Integer
end

"""Generates all |...SmS⟩|lml⟩ states up to and including l=lmax"""
function SmS_lookup_generator(lmax)
    lookup = Vector{SmS_ket}()
    for S₁ in 1:1, S₂ in 1:1
        Γ=Γ_nos(S₁,S₂)
        for S in abs(S₁-S₂):1:(S₁+S₂), l in 0:1:lmax
            mod(l+S,2)==0 || continue # eigenvalue of XN=(-1)^(l+S)==1 (Bosonic)
            for mS in -S:1:S, ml in -l:1:l
                push!(lookup, SmS_ket(Γ,S,mS,l,ml))
            end
        end
    end
    lookup
end

end # module
