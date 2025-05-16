##
## Particle pushing routines
##
abstract type AbstractForcing; end
abstract type AbstractPusher; end

# If a force is not defined for a given particle it does not act on it
force(::AbstractForcing, s) = zero(momentum(s))

# The simplest force does not act on anything
struct NullForcing <: AbstractForcing; end


struct CombinedForcing{T<:Tuple} <: AbstractForcing
    tpl::T    
end

CombinedForcing(a, b...) = CombinedForcing((a, b...))

force(c::CombinedForcing, s) = _force(c.tpl, s)

_force(tpl::Tuple, s) = force(first(tpl), s) + _force(Base.tail(tpl), s)
_force(tpl::Tuple{}, s) = zero(momentum(s))


# A restricted forcing acts only in one type of particle. This is mostly for debugging or if you
# want to put cont. losses only for positrons.
struct RestrictedForcing{T, F} <: AbstractForcing
    forcing::F
end
RestrictedForcing{T}(f) where T = RestrictedForcing{T, typeof(f)}(f)

force(r::RestrictedForcing{T}, s::T) where T = force(r.forcing, s)
force(r::RestrictedForcing, s) = zero(momentum(s))


struct RK2Pusher{F}
    forcing::F
end

function advance_particle(psh::RK2Pusher, y::ParticleState, Δt)
    (;forcing) = psh

    v1 = velocity(y)
    f1 = force(forcing, y)
    
    t2 = y.t + 2 * Δt / 3
    
    y2 = setproperties(y,
                       x = y.x + 2 * Δt * v1 / 3,
                       p = y.p + 2 * Δt * f1 / 3,
                       t = y.t + 2 * Δt / 3)
    v2 = velocity(y2)

    f2 = force(forcing, y2) 
    
    yf = setproperties(y,
                       x = y.x + Δt * (v1 / 4 + 3 * v2 / 4),
                       p = y.p + Δt * (f1 / 4 + 3 * f2 / 4),
                       t = y.t + Δt)
    
    return yf    
end


# A restricted pusher pushes only one type of particle.
struct RestrictedPusher{T, P}
    pusher::P
end
RestrictedPusher{T}(p) where T = RestrictedPusher{T, typeof(p)}(p)

advance_particle(psh::RestrictedPusher{T}, y::T, Δt) where T = advance_particle(psh.pusher, y, Δt)
advance_particle(psh::RestrictedPusher, y, Δt) = setproperties(y, t = y.t + Δt)

struct NullPusher; end
advance_particle(psh::NullPusher, y, Δt) = setproperties(y, t = y.t + Δt)
