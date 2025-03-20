const Positron = ParticleType{:positron}

struct PositronState{T} <: ParticleState{T}
    x::SVector{3, T}
    p::SVector{3, T}
    w::T

    t::T

    # s is the time to the next collision normalized by 1 / r
    s::T

    # r is a bound on the max. rate that the particle will experience during a timestep
    r::T

    active::Bool

    PositronState{T}(x::SVector{3, T}, v::SVector{3, T}, w::T=1.0, t::T=0.0, s::T=nextcoll(),
                     r::T=0.0, active=true) where T = new{T}(x, v, w, t, s, r, active)
    PositronState(x::SVector{3, T}, v::SVector{3, T}, w::T=1.0, t::T=0.0, s::T=nextcoll(),
                  r::T=0.0, active=true) where T = new{T}(x, v, w, t, s, r, active)
end


particle_type(::Type{PositronState{T}}) where T = Positron
particle_type(::PositronState) = Positron
new_particle(::Type{Positron}, x, v) = PositronState(x, v, 1.0, nextcoll(), true)

mass(p::PositronState) = co.electron_mass
mass(::Type{Positron}) = co.electron_mass
mass(::Positron) = co.electron_mass
" Charge in units of the elementary charge. "
charge(::Type{Positron}) = +1
speed(::Type{Positron}, eng) = co.c * sqrt(1 - (co.electron_mc2 / (co.electron_mc2 + eng))^2)
charge(::PositronState) = +1
# Valid for scalar or vector v
momentum_from_v(::Type{Positron}, v) = co.electron_mass * v / sqrt(1 - abs(v)^2 / co.c^2)
momentum_norm_from_kin(::Type{Positron}, kin) = sqrt((kin + co.electron_mc2)^2 - co.electron_mc2^2) / co.c



# advance_free_boris is defined in electron.jl and it works also for positrons.
@inline advance_free(s::PositronState, efield, bfield, Δt) = advance_free_boris(s, efield, bfield, Δt)

function lincomb(a::PositronState{T}, b::PositronState{T}, w::Number) where T
    PositronState{T}(a.x * w + b.x * (1 - w),
                     a.v * w + b.v * (1 - w),
                     a.w * w + b.w * (1 - w),
                     a.t * w + b.t * (1 - w))
end
