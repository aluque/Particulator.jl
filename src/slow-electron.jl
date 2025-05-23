#
# Non-relativistic dynamics
#
const SlowElectron = ParticleType{:slow_electron}

" Charge in units of the elementary charge. "
charge(::Type{SlowElectron}) = -1

struct SlowElectronState{T} <: ParticleState{T}
    x::SVector{3, T}
    v::SVector{3, T}
    w::T
    s::T    
    active::Bool
end

SlowElectronState(x, v, w=1.0, s=nextcoll(), active=true) = SlowElectronState(x, v, w, s, active)

particle_type(::Type{SlowElectronState{T}}) where T = SlowElectron
new_particle(::Type{SlowElectron}, x, v) = SlowElectronState(x, v, 1.0, nextcoll(), true)

mass(p::SlowElectronState) = co.electron_mass
mass(::Type{SlowElectron}) = co.electron_mass
mass(::SlowElectron) = co.electron_mass
charge(p::SlowElectronState) = -co.elementary_charge

kinenergy(p::SlowElectronState) = 0.5 * mass(p) * (p.v[1]^2 + p.v[2]^2 + p.v[3]^2)

@inline function advance_free_leapfrog(p::SlowElectronState, efield, Δt)
    # Leapfrog integration. Note that x and v are not synchronous.
    Δv = -(Δt * co.elementary_charge / mass(p)) .* efield(p.x)
    v1 = p.v .+ Δv

    SlowElectronState(p.x .+ Δt * v1, v1, p.w, p.s, p.active)
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
"""
@inline function advance_free_yoshida(p::SlowElectronState, efield, Δt)
    x1 = p.x + c1 * Δt * p.v
    v1 = p.v + d1 * Δt * (charge(p) / mass(p)) * efield(x1)
    x2 = x1 + c2 * Δt * v1
    v2 = v1 + d2 * Δt * (charge(p) / mass(p)) * efield(x2)
    x3 = x2 + c3 * Δt * v2
    v3 = v2 + d3 * Δt * (charge(p) / mass(p)) * efield(x3)
    x4 = x3 + c4 * Δt * v3
    v4 = v3
    
    SlowElectronState(x4, v4, p.w, p.s, p.active)
end

@inline advance_free(p::SlowElectronState, efield, Δt) = advance_free_yoshida(p, efield, Δt)


#
# Collisional processes and collision outcomes
#
struct Excitation{T} <: CollisionProcess
    threshold::T
end

struct Ionization{T} <: CollisionProcess
    threshold::T
end

struct Attachment{T} <: CollisionProcess
    threshold::T
end

struct Elastic{T} <: CollisionProcess
    mass_ratio::T
end

# This photo-emission. The energy loss and ionization are already included in
# other processes, so this only creates a new photon.
struct PhotoEmission{T} <: CollisionProcess;
    log_νmin::T
    log_νmax::T

    # `photon_multiplier` is used to produce a smoother distribution of photons.
    # The cross-section is multiplied by this factor and the weight is divided
    # by it.
    photon_multiplier::T

    # This is to reduce noise even beyond real physical noise, which can be
    # useful for comparing with fluid simulations.
    photon_weight::T
    
    function PhotoEmission(νmin, νmax, weight_scale, photon_weight)
        lmin, lmax = promote(log(νmin), log(νmax))
        new{typeof(lmin)}(lmin, lmax, weight_scale, photon_weight)
    end
end


function collide(c::Excitation, p::SlowElectronState{T}, energy) where T
    E1 = max(0, energy - c.threshold)
    vabs = sqrt(2 * E1 / mass(p))
    v1 = randsphere() .* vabs

    StateChangeOutcome(SlowElectronState{T}(p.x, v1, p.w, p.s, p.active))
end

function collide(c::Ionization, p::SlowElectronState{T}, energy) where T
    # Energy equipartitiom
    E1 = max(0, energy - c.threshold) / 2
    vabs = sqrt(2 * E1 / mass(p))

    v  = randsphere() .* vabs
    v1 = randsphere() .* vabs

    p1 = SlowElectronState{T}(p.x,  v, p.w, p.s, p.active)
    p2 = SlowElectronState{T}(p.x, v1, p.w, nextcoll(), p.active)
    NewParticleOutcome(p1, p2)
end

function collide(c::Attachment, p::SlowElectronState, energy)
    RemoveParticleOutcome(p)
end

function collide(c::Elastic, p::SlowElectronState{T}, energy) where T
    # Velocity in the center-of-mass frame
    v_cm = (c.mass_ratio / (1 + c.mass_ratio)) .* p.v
    v_final = randsphere() .* norm(p.v - v_cm) .+ v_cm

    StateChangeOutcome(SlowElectronState{T}(p.x, v_final, p.w, p.s, p.active))
end

function collide(c::PhotoEmission, p::SlowElectronState{T}, energy) where T
    nphot = rand(Poisson(p.w / c.photon_multiplier / c.photon_weight))
    MultiplePhotonOutcome(nphot, p, c)
end
