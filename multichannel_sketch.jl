#= Multichannel equations solver. Uses symmetrised |a⟩ basis states and all
of the interactions except for the Zeeman interaction.
All states up to and including the lmax parameter are used
Description last updated 12/08/20=#
push!(LOAD_PATH,raw"C:\Users\hirsc\OneDrive - Australian National University\PHYS4110\Code\1dsolver")

using Revise
using HalfIntegers, LinearAlgebra, StaticArrays, OrdinaryDiffEq, WignerSymbols
using Wigner9j, Potentials # my modules
using Unitful, UnitfulAtomic
using Plots

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
    l::HalfInteger
    J::HalfInteger
    mJ::HalfInteger
    XN::HalfInteger
end
struct a12_ket # unsymmetrised states, used to evaluate ̂Hₑₗ and ̂H_sd
    Γ::Γ_nos
    f::HalfInteger
    l::HalfInteger
    J::HalfInteger
    mJ::HalfInteger
end

"""Generates unsymmetrised state |a₁₂⟩
    Input: |a⟩ state containing Γ={α₁,α₂}
    Output: related |a₁₂⟩"""
function a12_maker(a::a_ket)
    return a12_ket(a.Γ,a.f,a.l,a.J,a.mJ)
end
"""Generates unsymmetrised state |a₂₁⟩
    Input: |a⟩ state containing Γ={α₁,α₂}
    Output: related |a₂₁⟩"""
function a21_maker(a::a_ket)
    Γ=Γ_nos(a.Γ.α₂, a.Γ.α₁) # flips α₁ and α₂
    return a12_ket(Γ,a.f,a.l,a.J,a.mJ)
end

"""Generates all |a⟩ states up to and including lmax
    Input: lmax=4
    Output: vector of a_kets
    Tested 12/08/20, I believe it is working correctly"""
function lookup_generator(lmax)
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
                        push!(lookup, a_ket(Γ,f,l,J,mJ,XN))
                    end
                end
            end
        end
    end
    return lookup
end

#TODO TISE solver. Input: lookup, lhs, rhs. Output: sol(R)
"""multichannel TISE solver
    Input: lookup vector, IC~SA{{[L],...,1,...}}, pot~|a⟩×|a⟩×[L]×[M]→[E],
    energy~[E], lhs~[L], rhs~[L]
    Output: sol(R) [where R ∈ (lhs,rhs)] ~ IC"""
function solver(lookup, IC, pot, ϵ, lhs, rhs; μ=0.5*4.002602u"u")
    # Check units of ϵ, μ, lhs, rhs
    @assert dimension(ϵ)==dimension(1u"J") "ϵ not an energy"
    @assert dimension(μ)==dimension(1u"g") "μ not a mass"
    @assert dimension(lhs)==dimension(1u"m") "lhs not a length"
    @assert dimension(rhs)==dimension(1u"m") "rhs not a length"
    # Check length and units of IC
    n = length(lookup) # number of channels
    @assert size(IC)[1]==2*n "Initial condition has wrong number of channels"
    if length(size(IC))==1 # IC a vector
        for i=1:n # check wavefunction entries
            @assert dimension(IC[i])==dimension(1u"m") "IC[$i] not a length"
        end
        for i=(n+1):2*n # check derivative entries
            @assert dimension(IC[i])==dimension(1) "IC[$i] not dimensionless"
        end
    elseif length(size(IC))==2 # IC a matrix of initial condition vectors
        for i=1:n, j=1:size(IC)[2] # check wavefunction entries
            @assert dimension(IC[i,j])==dimension(1u"m") "IC[$i,$j] not a length"
        end
        for i=(n+1):2*n, j=1:size(IC)[2] # check derivative entries
            @assert dimension(IC[i,j])==dimension(1u"m") "IC[$i,$j] not dimensionless"
        end
    else # IC not a 1D vector or 2D array
        error("IC not 1D or 2D")
    end
    # strip units from constants
    ϵ⁰, μ⁰ = austrip(ϵ), austrip(μ)
    ħ⁰ = austrip(1.0u"ħ")
    lhs⁰, rhs⁰ = austrip(lhs), austrip(rhs)
    # TISE differential equation
    function TISE(u,p,x)
        # Construct V(R) matrix
        V = zeros(n,n)u"hartree" # initialise
        for i=1:n, j=1:n
            V[i,j] = pot(lookup[i], lookup[j], x*1u"bohr", μ)
        end
        V⁰=austrip.(V) # strip units from V
        M = (-2μ⁰/ħ⁰^2)*(ϵ⁰*I-V⁰) # double derivative matix
        D = SMatrix{2*n,2*n}([0*I I
                              M 0*I])
        D*u # ⃗u' = D . ⃗u
    end
    # strip units from IC
    IC⁰ = austrip.(IC)
    # solve
    prob=ODEProblem(TISE,IC⁰,(lhs⁰,rhs⁰))
    sol_unitless=solve(prob,Tsit5())
    # add units back
    units = SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol = x -> sol_unitless(austrip(x)).*units
    return sol
end

""" test function for solver - runs for zero potential
    Tested successfully 12/08/2020"""
function test_solver()
    V(bra,ket,R,μ)=0u"hartree" # free
    lhs, rhs = 1.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=lookup_generator(0) # s-wave only
    n=length(lookup)
    IC=SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol=solver(lookup, IC, V, ϵ, lhs, rhs)
    Rs=LinRange(lhs,rhs,1000); vals = getindex.(sol.(Rs),1)
    plot(austrip.(Rs), austrip.(vals),title="If you see this, solver runs for the trivial case")
end

"""Rotational interaction
    Input: ⟨a'|, |a⟩, R~[L], μ~[M]
    Output: ⟨a'|̂H_rot|a⟩(R) ~ [E]"""
function H_rot(bra,ket,R,μ)
    # delta function
    bra == ket || return 0.0u"hartree"
    l=ket.l
    return uconvert(u"hartree", l*(l+1)*1.0u"ħ^2"/(2*μ*R^2))
end

# test of rotational interaction
function test_H_rot()
    V(bra,ket,R,μ)=H_rot(bra,ket,R,μ)
    lhs, rhs = 1.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=lookup_generator(1)
    n=length(lookup)
    IC=SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol=solver(lookup, IC, V, ϵ, lhs, rhs)
    Rs=LinRange(lhs,rhs,1000)
    vals1 = getindex.(sol.(Rs),1); vals10 = getindex.(sol.(Rs),10)
    vals2 = getindex.(sol.(Rs),2); vals11 = getindex.(sol.(Rs),11)
    vals=hcat(vals1,vals10,vals2,vals11)
    plot(austrip.(Rs), austrip.(vals),title="If you see this, solver runs for H_rot")
end


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
    return prefac_*prefac*(H_el_asym(a12_,a12,R)
                        +phase_*H_el_asym(a21_,a12,R)
                        +phase*H_el_asym(a12_,a21,R)
                        +phase_*phase*H_el_asym(a21_,a21,R)
                        )
end

# test of electronic interaction
function test_H_el()
    V(bra,ket,R,μ)=H_el(bra,ket,R)
    lhs, rhs = 3.0u"bohr", 100u"bohr"
    ϵ=1e-5u"hartree"
    lookup=lookup_generator(1) # s-wave only
    n=length(lookup)
    IC=SVector{2*n}(vcat(fill(1.0u"bohr",n),fill(1.0,n)))
    sol=solver(lookup, IC, V, ϵ, lhs, rhs)
    Rs=LinRange(lhs,rhs,1000)
    vals1 = getindex.(sol.(Rs),1); vals10 = getindex.(sol.(Rs),10)
    vals2 = getindex.(sol.(Rs),2); vals11 = getindex.(sol.(Rs),11)
    vals=hcat(vals1,vals10,vals2,vals11)
    plot(austrip.(Rs), austrip.(vals),title="If you see this, solver runs for H_el")
end

#TODO Spin-Dipole interaction
"""Symmetrised states ̂H_sd
    Input: ⟨a'|, |a⟩, R~[L]
    Output: ⟨a'|̂H_sd|a⟩(R) ~ [E]"""
function H_sd(bra::a_ket,ket::a_ket,R)
    #TODO D_a12'a12 coupling factor
    #TODO Vₚ(R) radial factor
end


#TODO hyperfine interaction

#TODO in2out_solver. Inputs: defaults=(see notebook). Output: matrix of IC results
