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

    "Time of last update"
    t::T

    "Normalized time to collision"
    s::T    

    "Bound on the max. rate experienced during a timestep"
    r::T
    
    "Active flag"
    active::Bool

    PhotonState{T}(x, p, w=1.0, t=0.0, s=nextcoll(), r=0.0, active=true) where T = new{T}(x, p, w, t, s, r, active)

end


particle_type(::Type{PhotonState{T}}) where T = Photon
new_particle(::Type{Photon}, x, p) = PhotonState(x, p, 1.0, 0.0, nextcoll(), 0.0, true)

mass(p::PhotonState) = 0
mass(::Type{Photon}) = 0
mass(::Photon) = 0
" Charge in units of the elementary charge. "
charge(::Type{Photon}) = 0
charge(::PhotonState) = 0

speed(::Type{Photon}, eng) = co.c

gamma(p::PhotonState) = Inf
momentum(p::PhotonState) = p.p
energy(p::PhotonState) = norm(p.p) * co.c

@inline function advance_free(p::PhotonState{T}, efield, bfield, Δt) where T
    v = (co.c / norm(p.p)) * p.p
    PhotonState{T}(p.x + Δt * v, p.p, p.w, p.t + Δt, p.s, p.r, p.active)
end
