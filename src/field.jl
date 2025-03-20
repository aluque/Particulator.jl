#=
Useful definitions for fields.
=#
struct HomogeneousField{T}
    v::SVector{3, T}
end

(f::HomogeneousField)(x, t) = f.v

struct DoubleLayerField{T}
    z1::T
    z2::T
    v::SVector{3, T}
end

(f::DoubleLayerField)(x, t) = f.z1 < x[3] < f.z2 ? f.v : zero(f.v)

struct ElectromagneticField{E, B} <: AbstractForcing
    e::E
    b::B
end

"""
Compute the instantaneous force on a particle.
"""
function force(em::ElectromagneticField, s::Union{ElectronState,PositronState})
    e = em.e(s.x, s.t)
    b = em.b(s.x, s.t)
    v = velocity(s)

    return charge(s) * co.elementary_charge * (e + cross(v, b))
end

force(em::ElectromagneticField, s::PhotonState) = zero(s.p)

