#=
Chebyshev interpolation in binary intervals.

We want to evaluate a function in an interval (a=0, b) but assume that we need a resolution
that around x is (roughly) proporional to x. What we do is divide the interval into halves,
further refining the left half. Then in each of the resulting intervals we approximate the function
with a sum of Chebyshev polinomials.
=#


"""
A struct to contain info to divide the interval (0, xmax) into approximately log-spaced
binary intervals.
"""
struct BinaryIntervals{T, N}
    "Intervals will be indexed from 0 to k"
    k::Int

    "The root interval is (0, xmax)"
    xmax::T
end


function interval(b::BinaryIntervals, i)
    if i == 0
        (l, r) = interval(b, 1)
        return (zero(l), l)
    end
    
    (;k, xmax) = b
    l = 2.0^(i - k - 1)
    r = 2 * l

    return (l * xmax, r * xmax)

end

"""
Compute Chebyshev nodes in the interval (`a`, `b`) for order `n`.
"""
function chebnodes(n::Int, a=-1, b=1)
    return [(a + b)/2 + (b - a)/2 * cos((2k + 1) * π / (2n)) for k in 0:(n - 1)]
end

"""
Compute Chebyshev nodes in interva `i` in the binary division `b`, with order `n`.
"""
function chebnodes(b::BinaryIntervals, i, n)
    (a, b) = interval(b, i)
    return chebnodes(n, a, b)
end


"""
Evaluate Chebyshev polynomials up to `N` and return their values as a tuple (T_0(ξ), T_1(ξ),...)
"""
@generated function chebval(ξ, ::Val{N}) where N
    if N == 1
        expr = quote
            return (one(ξ),)
        end
    elseif N == 2
        expr = quote
            return (one(ξ), ξ)
        end
    elseif N > 2
        expr = quote
            tpl = (one(ξ), ξ)
        end
        
        for i in 3:N
            push!(expr.args, quote
                      tpl = (tpl..., 2 * ξ * tpl[end] - tpl[end - 1])
                  end)
        end
        push!(expr.args, quote
                  return tpl
                  end)
    end
    return expr        
end

# This way of defining chebval kept resulting in type instabilities in some cases.
# I don't know way the compiler was not able of inferring.
# chebval(ξ, v::Val{1}, t::Nothing=nothing) = return (one(ξ),)
# chebval(ξ, v::Val{2}, t::Nothing=nothing) = return (one(ξ), ξ)

# function chebval(ξ, v::Val{N}, t::Nothing=nothing) where {N}
#     return chebval(ξ, v, (one(ξ), ξ))
# end


# function chebval(ξ::T, v::Val{N}, t::NTuple{M, T}) where {T, N, M}
#     @assert M >= 2

#     if N == M
#         return t
#     else    
#         u = 2 * ξ * t[end] - t[end - 1]
#         return chebval(ξ, v, (t..., u))
#     end
# end

"""
Evaluate at `x` a chebysev expansion contained in `a` for the binary division specified by `b`.
If multiple evaluations are requred for the same `x` part of the computations may be stored in
`pre`, which is computed from `precheb`.
"""
function chebeval(x::T, b::BinaryIntervals{T, N}, a, pre=nothing) where {T, N}
    if isnothing(pre)
        pre = precheb(x, b)
    end

    (i, t) = pre

    a1 = a[i + 1]
    
    return sum(ntuple(j -> a1[j] * t[j], Val(N)))
end

function precheb(x::T, b::BinaryIntervals{T, N})
    (;k, xmax) = b
    x1 = x / xmax
    (s, l) = frexp(x1)
    i = x1 == 0 ? 0 : l + k

    if i > 0
        ξ = 4s - 3
    else
        ξ = 2^(k + 1) * x1 - 1
        i = 0
    end
    
    t = chebval(ξ, Val{N}())

    return (i, t)    
end

"""
Find a Chebyshev expansion for function `f` in each of the intervals of `b` and return it as a
vector of SVectors, one for each interval.
"""
function fit(f, b::BinaryIntervals{T, N}) where {N, T}
    (;k) = b
    a = SVector{N, T}[]

    for i in 0:k
        (l, r) = interval(b, i)
        x = chebnodes(N, l, r)
        fx = f.(x)
        ξ = @.((2 * x - (l + r)) / (r - l))
        A = hcat(SVector{N}.(chebval.(ξ, Val(N)))...)'
        coeffs = A \ fx
        push!(a, SVector{N}(coeffs))
    end

    return a
end
