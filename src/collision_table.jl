abstract type AbstractCollisionTable{T, C <: Tuple}; end

"""
Struct with a set of collisions evaluated at the same energy grid.

The rate bound can either be a number or another struct. It is used in the form

ratebound(::CollisionTable, energy).

If rate bound is a number the  assumption is that the total collision rate is independent of
energy (because the space is filled with NullCollisions).

The rate bound at a given energy must be an upper limit of the energy-dependent rate during a timestep.
"""
@kwdef struct CollisionTable{T, C <: Tuple, V <: AbstractRange{T},
                             A <: AbstractMatrix{T}, TR <: Union{Vector{T}, T}} <: AbstractCollisionTable{T, C}
    "Each type of collision"
    proc::C

    "Energy grid"
    energy::V

    "Array with collision rates for each process."
    rate::A

    "Upper limit of the total collision rate during a time-step"
    ratebound::TR
end

Base.length(c::CollisionTable) = length(c.proc)

ratebound(c::CollisionTable, eng, pre=nothing) = ratebound(c.ratebound, c, eng, pre)
ratebound(x::Number, c, eng, pre=nothing) = x

function ratebound(v::Vector, c, eng, pre::Nothing)
    pre = presample(c, nothing, eng)
    ratebound(v, c, eng, pre)
end

function ratebound(v::Vector, c, eng, pre)
    k, w = pre
    return w * v[k] + (1 - w) * v[k + 1]
end

maxenergy(c::CollisionTable) = last(c.energy)

# Checking for collisions involves many tests but some computations are
# common to all of them. Here we store them in a generic way. energy is passed
# as an optimization.
presample(c::CollisionTable, state, energy) = indweight(c.energy, energy)

# Obtain probability rate of process j
function rate(c::CollisionTable, j, pre)
    k, w = pre

    return w * c.rate[j, k] + (1 - w) * c.rate[j, k + 1]
end


"""
A collision table based on Chebyshev expansions in a set of binary intervals.
"""
@kwdef struct ChebyshevCollisionTable{T, N, C <: Tuple, A <: AbstractArray, B} <: AbstractCollisionTable{T, C}
    "Each type of collision"
    proc::C

    "Binary divisions"
    b::BinaryIntervals{T}

    "Chebyshev expansions for each interval for each process"
    rate::A

    "Chebyshev expansions for the upper limit of the total collision rate during a time-step"
    ratebound::B
end


Base.length(c::ChebyshevCollisionTable) = length(c.proc)

maxenergy(c::ChebyshevCollisionTable) = c.b.xmax

presample(c::ChebyshevCollisionTable{T, N}, state, energy) where {T, N} = precheb(energy, c.b, Val{N}())

# Obtain probability rate of process j
function rate(c::ChebyshevCollisionTable{T, N}, j, pre) where {T, N}
    (i, t) = pre

    # i = interval
    # j = process
    # k = order
    @inbounds return sum(ntuple(k -> @inbounds(c.rate[k, j, i + 1] * t[k]), Val(N)))    
end

function ratebound(c::ChebyshevCollisionTable, eng, pre::Nothing=nothing)
    pre = presample(c, nothing, eng)
    ratebound(c, eng, pre)
end

function ratebound(c::ChebyshevCollisionTable{T, N}, eng, pre) where {T, N}
    (i, t) = pre

    # i = interval
    # j = process
    # k = order
    return sum(ntuple(k -> c.ratebound[k, i + 1] * t[k], Val(N)))    
end


"""
Build a (Chebyshev) collision table from a list of `processes` for a given particle type `particle`.
`processes` must be a list of tuples `(coeff, process)` where `coeff` is the multiplicative factor
that converts cross section per atom to a inverse mean-free-path (typically this is a density of
scatterers).
"""
function collision_table_from_processes(processes, particle, Fdt;
                                        safety=1.1,
                                        nintervals=32,
                                        order=3,
                                        # This upper limit ensures that 1 keV is at an interval
                                        # boundary and discontinuities
                                        # there are ok.
                                        max_energy = 1e3 * co.eV * 2^18)
    nprocs = length(processes)

    b = BinaryIntervals{Float64}(nintervals, max_energy)
    a = Array{Float64, 3}(undef, (order, nprocs, nintervals + 1))

    for j in 1:nprocs
        (coef, proc) = processes[j]
        a1 = chebfit(eng -> coef * speed(particle, eng) * totalcs(proc, eng), b, order)
        a[:, j, :] .= reinterpret(reshape, Float64, a1)
    end
    
    rb = compratebound(particle, a, b, Fdt, safety)
        
    eng = LogLinRange(log(1e-2 * co.eV), log(0.9999 * max_energy), 100_000)
    
    v1 = speed.(Ref(particle), eng)
    p = sortperm(1:nprocs, by=i ->    
                 processes[i][1] * maximum(totalcs.(Ref(processes[i][2]), eng) .* v1))
    
    proc = tuple([processes[p[i]][2] for i in 1:nprocs]...)
    a = a[:, p, :]
    
    return ChebyshevCollisionTable{Float64, order, typeof(proc), typeof(a),
                                   typeof(rb)}(;proc=proc, b=b, rate=a, ratebound=rb)    
end


function compratebound(particle, a, b, Fdt, safety)
    mc2 = co.electron_mass * co.c^2
    γ(::Type{Electron}, eng) = 1 + eng / mc2
    γ(::Type{Positron}, eng) = 1 + eng / mc2
    γ(::Type{Photon}, eng) = 0

    rb = chebfit(b, size(a, 1)) do eng
        df = sum(j -> Particulator.chebdiff(eng, b, @view a[:, j, :]), axes(a, 2))
        f = sum(j -> Particulator.chebeval(eng, b, @view a[:, j, :]), axes(a, 2))
        return f + γ(particle, eng)^3 * speed(particle, eng) * Fdt * abs(df) * safety
    end

    return rb
end
