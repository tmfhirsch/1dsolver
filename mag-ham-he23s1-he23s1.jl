#=This code aims to write the magnetic hamiltonian H=-Bμg(S1z+S2z)/ħ
in the bases |m1 m2> and |s m2> for the s1=s2=1 case
=#
using Plots, SparseArrays, WignerSymbols

#Hamiltonian scalars, all defined to be 1
const B=1
const μ=1
const g=1
const ħ=1

################################################################################
#Individual spin basis |m1 m2>
m1m2Basis=[(m1,m2) for m2=-1:1, m1=-1:1]

#Hamiltonian matrix element: Function <m1' m2'| H |m1 m2>==<b1| H |b2>
function m1m2H(b1,b2)
    return Int(-B*μ*g*(sum(b2))*(b1==b2 ? 1 : 0)/ħ)
end

m1m2HamCheck=[(b1,b2) for b1 in vec(m1m2Basis), b2 in vec(m1m2Basis)]
m1m2Ham=[m1m2H(b1,b2) for b1 in vec(m1m2Basis), b2 in vec(m1m2Basis)]

#Plot m1m2 Hamiltonian
#spy(sparse(m1m2Ham),legend=nothing,axis=nothing)

################################################################################
#|s m> basis states (dimension one array)
smBasis=[]
for s=0:2
    for m=-s:s
        push!(smBasis,tuple(s,m))
    end
end

#|s m> basis Hamiltonian
function smH(b1,b2)
    summand(si,mi,sj,mj,m1,m2)=-B*g*μ/ħ*conj(clebschgordan(1,m1,1,m2,si,mi))*
    clebschgordan(1,m1,1,m2,sj,mj)*(m1+m2)
    sum=0
    for m1m2 in m1m2Basis
        m1,m2 = m1m2[1],m1m2[2]
        si,mi = b1[1],b1[2]
        sj,mj = b2[1],b2[2]
        sum += summand(si,mi,sj,mj,m1,m2)
    end
    return sum
end

smHam=[smH(b1,b2) for b1 in vec(smBasis), b2 in vec(smBasis)]
