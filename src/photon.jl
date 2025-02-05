#
# Relativistic dynamics
#
const Photon = ParticleType{:photon}

struct PhotonState{T} <: ParticleState{T}
    "Photon position"
    x::SVector{3, T}

    "Photon momentum"
    p::SVector{3, T}

    "Weight"
    w::T

    "Normalized time to collision"
    s::T    

    "Time of last update"
    t::T
    
    "Active flag"
    active::Bool
end

PhotonState(x, p, w=1.0, s=nextcoll(), t=0.0, active=true) = PhotonState(x, p, w, s, t, active)

particle_type(::Type{PhotonState{T}}) where T = Photon
new_particle(::Type{Photon}, x, p) = PhotonState(x, p, 1.0, nextcoll(), 0.0, true)

mass(p::PhotonState) = 0
mass(::Type{Photon}) = 0
mass(::Photon) = 0
" Charge in units of the elementary charge. "
charge(::Type{Photon}) = 0
charge(::PhotonState) = 0

gamma(p::PhotonState) = Inf
momentum(p::PhotonState) = p.p
energy(p::PhotonState) = norm(p.p) * co.c

@inline function advance_free(p::PhotonState, efield, bfield, Δt)
    v = (co.c / norm(p.p)) * p.p
    PhotonState(p.x + Δt * v, p.p, p.w, p.s, p.t + Δt, p.active)
end
