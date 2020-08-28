#= Module for the Wigner 9j symbol, implemented following the formula decsribed
in Wei (1998) (doi: 10.1063/1.168745)
Description last updated 6/08/20=#

module Wigner9j
export wigner9j

using HalfIntegers, WignerSymbols

""" Wigner 9-j symbol calculator, based off Wei (1998)\n
    Inputs: a,...,j into \n⌈a b c⌉
                            \n{d e f}
                            \n⌊g h j⌋\n
    Outputs: Evaluation of the Wigner 9-j symbol
    Tested successfully against Wei Table 1 10/08/20"""
function wigner9j(a,b,c,d,e,f,g,h,j)
    # Checks that each jᵢ is a (half-)integer, as all angular momenta should be
    for i in (a,b,c,d,e,f,g,h,j)
        i = try convert(HalfInt,i)
        catch err
            throw(DomainError(i,"Not a Integer or Half Integer"))
        end
    end
    # convert each argument into (half-)integer
    (a,b,c,d,e,f,g,h,j)=convert.(HalfInt,(a,b,c,d,e,f,g,h,j))
    # Δ factors
    eval = Δ(a,b,c)*Δ(d,e,f)*Δ(g,h,j)*Δ(a,d,g)*Δ(b,e,h)*Δ(c,f,j)
    # Construct sum
    I₁=max(abs(h-d),abs(b-f),abs(a-j))
    I₂=min(h+d,b+f,a+j)
    sum=0
    for k in I₁:1:I₂
        # construct sum term
        term=(-1)^(2*k)
        term*=(2*k+1)
        term*=Wei6(a,b,c,f,j,k)*Wei6(f,d,e,h,b,k)*Wei6(h,j,g,a,d,k)
        # add term to sum
        sum+=term
        # reset term
        term=0
    end
    # multiply by sum
    eval*= sum
    return eval
end

"""     Recursive binomial function, to reduce overflow when calculating large
        9j symbols. Uses ⁿCₖ=ⁿCₖ₋₁*((n+1)-k)/k
        Inputs: n, k integers
        Outputs: Binomial (ⁿCₖ)"""
function WeiBinomial(n::Int, k::Int)
    0 <= k <= n || return 0
    x=big(1)
    while k>0
        x*=((n+1)-k)//k
        k-=1
    end
    #k=0
    return x
end

""" Test recursive binomial # Tested successfully 10/08/20"""
function WeiBinomTest()
    rands1=convert.(Int, ceil.(rand(100).*20))
    rands2=convert.(Int, ceil.(rand(100).*20))
    for n in rands1, k in rands2
        WeiBinomial(n,k)==binomial(n,k) || return("Fail for $n, $k")
    end
    println("Success!")
end


"""     DEPRICATED: function Δ in WignerSymbols does this
        Δ(a b c) as defined below eqn (1) in Wei (1998)"""
function WeiΔ(a,b,c)
    if !δ(a,b,c) # first check triangle condition is upheld
        throw(DomainError((a,b,c),"Do not satisfy the triangle relation"))
    end
    return sqrt(factorial(a+b-c)
                *factorial(a-b+c)
                *factorial(-a+b+c)
                /factorial(a+b+c+1))
end

""" [m1 m2 m3
    \n   m4 m5 m6] from (2) in Wei (1998)"""
function Wei6(m1,m2,m3,m4,m5,m6)
    p=max(m1+m5+m6,
          m2+m4+m6,
          m3+m4+m5,
          m1+m2+m3)
    q=min(m1+m2+m4+m5,
          m1+m3+m4+m6,
          m2+m3+m5+m6)
    sum=0
    for n in p:1:q
        # build construct term to add to sum
        term=(-1)^n
        # define binomial to deal with the HalfInt data type
        bn(x,y)=WeiBinomial(convert(Int,x),convert(Int,y))
        term*=bn(n+1,n-m1-m5-m6)
        term*=bn(m1+m5-m6,n-m2-m4-m6)
        term*=bn(m1-m5+m6,n-m3-m4-m5)
        term*=bn(-m1+m5+m6,n-m1-m2-m3)
        # add term to sum
        sum += term
        # reset term/
        term = 0
    end
    return sum
end

end #module
