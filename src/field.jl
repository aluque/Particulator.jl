#=
Useful definitions for fields.
=#
struct HomogeneousField{T}
    v::SVector{3, T}
end

(f::HomogeneousField)(x) = f.v

struct DoubleLayerField{T}
    z1::T
    z2::T
    v::SVector{3, T}
end

(f::DoubleLayerField)(x) = f.z1 < x[3] < f.z2 ? f.v : zero(f.v)
