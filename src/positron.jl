const Positron = ParticleType{:positron}

struct PositronState{T} <: ParticleState{T}
    x::SVector{3, T}
    v::SVector{3, T}
    w::T
    s::T    
    t::T
    active::Bool
end

PositronState(x, v, w=1.0, s=nextcoll(), t=0.0, active=true) = PositronState(x, v, w, s, t, active)


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

# Only kinetic energy
energy(p::PositronState) = (gamma(p) - 1) * mass(p) * co.c^2

# advance_free_boris is defined in electron.jl and it works also for positrons.
@inline advance_free(p::PositronState, efield, bfield, Δt) = advance_free_boris(p, efield, bfield, Δt)
