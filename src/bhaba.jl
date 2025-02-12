#=
Bhaba scattering (ionization) for positrons.

Note: Currently I am not using continuous slowing-down (CSDA) so one has to set a low tcut to include
all energy losses.
=#
struct Bhaba{T}
    Z::Int
    tcut::T
end

function collide(b::Bhaba, pos::PositronState{T}, eng) where T
    # This is the energy of the liberated electron.
    mc2 = co.electron_mc2

    E2 = sample_secondary_energy(b, eng)

    E0 = eng
    E1 = E0 - E2

    p0 = sqrt(E0^2 + 2 * mc2 * E0) / co.c
    p1 = sqrt(E1^2 + 2 * mc2 * E1) / co.c
    p2 = sqrt(E2^2 + 2 * mc2 * E2) / co.c

    γ1 = 1 + (E1 / mc2)
    γ2 = 1 + (E2 / mc2)

    vnorm = norm(pos.v)
    v1norm = p1 / γ1 / co.electron_mass
    v2norm = p2 / γ2 / co.electron_mass

    # Lehtinen 1999
    cosθ1 = sqrt(E1 * (E0 + 2mc2) / (E0 * (E1 + 2mc2)))    
    cosθ2 = sqrt(E2 * (E0 + 2mc2) / (E0 * (E2 + 2mc2)))
    
    ϕ = 2π * rand()

    v1 = turn(pos.v, cosθ1,  ϕ, v1norm)
    v2 = turn(pos.v, cosθ2, -ϕ, v2norm)

    NewParticleOutcome(PositronState{T}(pos.x, v1, pos.w, pos.t),
                       ElectronState{T}(pos.x, v2, pos.w, pos.t))
end

function totalcs(b::Bhaba, eng)
    (;Z, tcut) = b
    
    mc2 = co.electron_mc2
    γ = 1 + eng / mc2

    β = sqrt((γ^2 - 1) / γ^2)
    x = tcut / eng
    y = 1 / (γ + 1)
    (B1, B2, B3, B4) = bhaba_bs(y)
    
    A = (1 / x - 1) / β^2 + B1 * log(x) + B2 * (1 - x) - B3 * (1 - x^2) / 2 + B4 * (1 - x^3) / 3

    return max(0, 2π * co.r_e^2 * Z * A / (γ - 1))
end

"""
Sample energy of the liberated electron in Bhaba scattering.
Based on the GEANT4 Phys. Ref. Manual v. 11.3 p. 130.
"""
function sample_secondary_energy(b::Bhaba, eng)
    (;tcut) = b

    ϵ0 = tcut / eng
    mc2 = co.electron_mc2
    γ = 1 + eng / mc2

    β = sqrt((γ^2 - 1) / γ^2)
    y = 1 / (γ + 1)

    B0 = γ^2 / (γ^2 - 1)
    (B1, B2, B3, B4) = bhaba_bs(y)

    local ϵ
    accept = false
    while !accept
        r = rand()
        ϵ = ϵ0 / (1 - r + ϵ0 * r)
        g1 = B0 + B1 * ϵ + B2 * ϵ^2 + B3 * ϵ^3 + B4 * ϵ^4
        g2 = B0 + B1 * ϵ0 + B2 * ϵ0^2 + B3 * ϵ0^3 + B4 * ϵ0^4
        g = g1 / g2

        accept = rand() < g
    end

    return ϵ * eng
end

    
function bhaba_bs(y)
    B1 = 2 - y^2
    B2 = (1 - 2y) * (3 + y^2)
    B3 = (1 - 2y)^2 + (1 - 2y)^3
    B4 = (1 - 2y)^3

    return (B1, B2, B3, B4)
end

