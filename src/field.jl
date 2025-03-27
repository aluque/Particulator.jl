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


"""
An electric field with a close form that approximate the electric field of a double layer
but is confined on x and y.
"""
struct ConfinedDoubleLayerField{T}
    sx::T
    sy::T
    sz::T
    ez0::T
end

function (f::ConfinedDoubleLayerField)(r, t)
    (;ez0, sx, sy, sz) = f
    (x, y, z) = r
    
    ex = -(ez0 * exp(-(x^2 / (2*sx^2) + y^2 / (2*sy^2) + z^2 / (2*sz^2))) * x * z) / sx^2
    ey = -(ez0 * exp(-(x^2 / (2*sx^2) + y^2 / (2*sy^2) + z^2 / (2*sz^2))) * y * z) / sy^2
    ez = (ez0 * exp(-(x^2 / (2*sx^2) + y^2 / (2*sy^2) + z^2 / (2*sz^2))) 
          - (ez0 * exp(-(x^2 / (2*sx^2) + y^2/(2*sy^2) + z^2 / (2*sz^2))) * z^2) / sz^2)

    return SA[ex, ey, ez]
end

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

