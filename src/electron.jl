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
    v::SVector{3, T}
    w::T
    s::T    
    active::Bool
end

ElectronState(x, v, w=1.0, s=nextcoll(), active=true) = ElectronState(x, v, w, s, active)

particle_type(::Type{ElectronState{T}}) where T = Electron
new_particle(::Type{Electron}, x, v) = ElectronState(x, v, 1.0, nextcoll(), true)

mass(p::ElectronState) = co.electron_mass
mass(::Type{Electron}) = co.electron_mass
mass(::Electron) = co.electron_mass
" Charge in units of the elementary charge. "
charge(::Type{Electron}) = -1
charge(::ElectronState) = -1

gamma(p::ElectronState) = 1 / (sqrt(1 - dot(p.v, p.v) / co.c^2))
momentum(p::ElectronState) = gamma(p) * mass(p) * p.v

# Only kinetic energy
energy(p::ElectronState) = (gamma(p) - 1) * mass(p) * co.c^2

@inline function advance_free_boris(p::ElectronState, efield, bfield, Δt)
    γ = gamma(p)
    u = γ * p.v
    q = co.elementary_charge * charge(p)
    m = mass(p)
    
    u1 = u + ((q * Δt) / 2m) * efield(p.x)
    γB = sqrt(1 + dot(u1, u1) / co.c^2)

    Ω = (q / m) * bfield(p.x)
    u2 = ((1 - (norm(Ω) * Δt / 2γB)^2) * u1
          + cross(u1, Ω * Δt)
          + ((Δt / γB)^2 / 2) * dot(u1, Ω) * Ω) / (1 + (norm(Ω) * Δt / 2γB)^2)

    uf = u2 + ((q * Δt) / 2m) * efield(p.x)
    γf = sqrt(1 + dot(uf, uf) / co.c^2)
    v1 = uf / γf

    ElectronState(p.x + Δt * v1, v1, p.w, p.s, p.active)
end

# Use only if B = 0
@inline function advance_free_leapfrog(p::ElectronState, efield, bfield, Δt)
    # Leapfrog integration. Note that x and v are not synchronous.
    Δv = -(Δt * co.elementary_charge / mass(p)) .* efield(p.x)
    v1 = p.v .+ Δv

    ElectronState(p.x .+ Δt * v1, v1, p.w, p.s, p.active)
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
@inline function advance_free_yoshida(p::ElectronState, efield, bfield, Δt)
    x1 = p.x + c1 * Δt * p.v
    v1 = p.v + d1 * Δt * (charge(p) / mass(p)) * efield(x1)
    x2 = x1 + c2 * Δt * v1
    v2 = v1 + d2 * Δt * (charge(p) / mass(p)) * efield(x2)
    x3 = x2 + c3 * Δt * v2
    v3 = v2 + d3 * Δt * (charge(p) / mass(p)) * efield(x3)
    x4 = x3 + c4 * Δt * v3
    v4 = v3
    
    ElectronState(x4, v4, p.w, p.s, p.active)
end

@inline advance_free(p::ElectronState, efield, bfield, Δt) = advance_free_boris(p, efield, bfield, Δt)

