const Positron = ParticleType{:positron}

struct PositronState{T} <: ParticleState{T}
    x::SVector{3, T}
    v::SVector{3, T}
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
new_particle(::Type{Positron}, x, v) = PositronState(x, v, 1.0, nextcoll(), true)

mass(p::PositronState) = co.electron_mass
mass(::Type{Positron}) = co.electron_mass
mass(::Positron) = co.electron_mass

" Charge in units of the elementary charge. "
charge(::Type{Positron}) = 1
charge(::PositronState) = 1

gamma(p::PositronState) = 1 / (sqrt(1 - dot(p.v, p.v) / co.c^2))
momentum(p::PositronState) = gamma(p) * mass(p) * p.v
speed(::Type{Positron}, eng) = co.c * sqrt(1 - (co.electron_mc2 / (co.electron_mc2 + eng))^2)

# Only kinetic energy
energy(p::PositronState) = (gamma(p) - 1) * mass(p) * co.c^2

# advance_free_boris is defined in electron.jl and it works also for positrons.
@inline advance_free(p::PositronState, efield, bfield, Δt) = advance_free_boris(p, efield, bfield, Δt)
