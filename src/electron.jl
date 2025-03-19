#
# Relativistic dynamics
#
const Electron = ParticleType{:electron}

# Classical electron radius, about 2.8e-15 m
const r_e = co.elementary_charge^2 / (co.electron_mass * co.c^2) / (4π * co.epsilon_0)

# Bohr radius
const a_0 = co.hbar / (co.electron_mass * co.c * co.fine_structure)

    
struct ElectronState{T} <: ParticleState{T}
    x::SVector{3, T}
    p::SVector{3, T}
    w::T

    t::T

    # s is the time to the next collision normalized by 1 / r
    s::T

    # r is a bound on the max. rate that the particle will experience during a timestep
    r::T

    active::Bool

    # When creating ElectronStates it is better to leave out the s, r, active fields so they get
    # apropriate values
    ElectronState{T}(x::SVector{3, T}, p::SVector{3, T}, w::T=1.0, t::T=0.0, s::T=nextcoll(),
                     r::T=0.0, active=true) where T = new{T}(x, p, w, t, s, r, active)
    ElectronState(x::SVector{3, T}, p::SVector{3, T}, w::T=1.0, t::T=0.0, s::T=nextcoll(),
                  r::T=0.0, active=true) where T = new{T}(x, p, w, t, s, r, active)
end



particle_type(::Type{ElectronState{T}}) where T = Electron
particle_type(::ElectronState) = Electron
new_particle(::Type{Electron}, x, p) = ElectronState(x, p)

mass(p::ElectronState) = co.electron_mass
mass(::Type{Electron}) = co.electron_mass
mass(::Electron) = co.electron_mass
" Charge in units of the elementary charge. "
charge(::Type{Electron}) = -1
speed(::Type{Electron}, eng) = co.c * sqrt(1 - (co.electron_mc2 / (co.electron_mc2 + eng))^2)
charge(::ElectronState) = -1
# Valid for scalar or vector v
momentum_from_v(::Type{Electron}, v) = co.electron_mass * v / sqrt(1 - norm(v)^2 / co.c^2)
momentum_norm_from_kin(::Type{Electron}, kin) = sqrt((kin + co.electron_mc2)^2 - co.electron_mc2^2) / co.c
# kinetic energy
kinenergy(s::ElectronState) = (sqrt(co.electron_mc2^2 + co.c^2 * dot(s.p, s.p)) - co.electron_mc2)
gamma(s::ElectronState) = 1 + co.c^2 * dot(s.p, s.p) / co.electron_mc2^2
momentum(s::ElectronState) = s.p
velocity(s::ElectronState) = s.p * (1 / (co.electron_mass * gamma(s)))


@inline function advance_free_boris(s, efield, bfield, Δt)
    γ = gamma(s)
    v = velocity(s)
    u = γ * v
    q = co.elementary_charge * charge(s)
    m = mass(s)
    
    u1 = u + ((q * Δt) / 2m) * efield(s.x)
    γB = sqrt(1 + dot(u1, u1) / co.c^2)

    Ω = (q / m) * bfield(s.x)
    u2 = ((1 - (norm(Ω) * Δt / 2γB)^2) * u1
          + cross(u1, Ω * Δt)
          + ((Δt / γB)^2 / 2) * dot(u1, Ω) * Ω) / (1 + (norm(Ω) * Δt / 2γB)^2)

    uf = u2 + ((q * Δt) / 2m) * efield(s.x)
    γf = sqrt(1 + dot(uf, uf) / co.c^2)
    v1 = uf / γf
    p1 = momentum_from_v(particle_type(s), v1)
    
    typeof(s)(s.x + Δt * v1, p1, s.w, s.t + Δt, s.s, s.r, s.active)
end

# Use only if B = 0
@inline function advance_free_leapfrog(s, efield, bfield, Δt)
    # Leapfrog integration. Note that x and v are not synchronous.
    Δv = -(Δt * co.elementary_charge / mass(s)) .* efield(s.x)
    v1 = velocity(s) .+ Δv
    p1 = momentum_from_velocity(particle_type(s), v1)
    typeof(p)(s.x .+ Δt * v1, p1, s.w, s.t + Δt, s.s, s.r, s.active)
end


const w0 = -2^(1/3) / (2 - 2^(1/3))
const w1 = 1 / (2 - 2^(1/3))
const c1 = w1 / 2
const c2 = (w0 + w1) / 2
const c3 = c2
const c4 = c1
const d1 = w1
const d2 = w0
const d3 = d1

"""
Yoshida 4th order integrator.
Use only if B = 0.
"""
@inline function advance_free_yoshida(s, efield, bfield, Δt)
    v = velocity(s)
    x1 = s.x + c1 * Δt * v
    v1 = v + d1 * Δt * (charge(s) / mass(s)) * efield(x1)
    x2 = x1 + c2 * Δt * v1
    v2 = v1 + d2 * Δt * (charge(s) / mass(s)) * efield(x2)
    x3 = x2 + c3 * Δt * v2
    v3 = v2 + d3 * Δt * (charge(s) / mass(s)) * efield(x3)
    x4 = x3 + c4 * Δt * v3
    v4 = v3

    p4 = momentum_from_velocity(particle_type(s), v4)
    typeof(s)(x4, p4, s.w, s.t + Δt, s.s, s.r, s.active)
end

@inline advance_free(p::ElectronState, efield, bfield, Δt) = advance_free_boris(p, efield, bfield, Δt)

"""
Produce a linear combination between two `ElectronState`s, `a` and `b`, where `a` is assigned a weight
`w` and `b` a weight `(1 - w)`.
"""
function lincomb(a::ElectronState{T}, b::ElectronState{T}, w::Number) where T
    ElectronState{T}(a.x * w + b.x * (1 - w),
                     a.p * w + b.p * (1 - w),
                     a.w * w + b.w * (1 - w),
                     a.t * w + b.t * (1 - w))
end
