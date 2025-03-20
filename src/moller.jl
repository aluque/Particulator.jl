#=
Moller scattering.

Moller is an alternative to RBEB scattering, which becomes a better approximation for energies far
above the ionization threshold. It can also be complemented by continuum losses.
=#

struct Moller{T}
    Z::Int
    tcut::T
end

function collide(m::Moller, el::ElectronState{T}, eng) where T
    # This is the energy of the liberated electron.
    mc2 = co.electron_mc2

    E2 = sample_secondary_energy(m, eng)

    E0 = eng
    E1 = E0 - E2

    p0 = sqrt(E0^2 + 2 * mc2 * E0) / co.c
    p1 = sqrt(E1^2 + 2 * mc2 * E1) / co.c
    p2 = sqrt(E2^2 + 2 * mc2 * E2) / co.c

    # Lehtinen 1999
    cosθ1 = sqrt(E1 * (E0 + 2mc2) / (E0 * (E1 + 2mc2)))    
    cosθ2 = sqrt(E2 * (E0 + 2mc2) / (E0 * (E2 + 2mc2)))
    
    ϕ = 2π * rand()

    p1vec = turn(el.p, cosθ1,  ϕ, p1)
    p2vec = turn(el.p, cosθ2, -ϕ, p2)

    NewParticleOutcome(ElectronState{T}(pos.x, p1vec, pos.w, pos.t),
                       ElectronState{T}(pos.x, p2vec, pos.w, pos.t))
end


function totalcs(m::Moller, eng)
    (;Z, tcut) = m
    
    mc2 = co.electron_mc2
    γ = 1 + eng / mc2

    β2 = (γ^2 - 1) / γ^2
    x = tcut / eng
    
    A = (((γ - 1)^2 / γ^2) * (1/2 - x)
         + 1 / x
         - 1 / (1 - x)
         - ((2γ - 1) / γ^2) * log((1 - x) / x))

    return max(0, 2π * co.r_e^2 * Z * A / (γ - 1) / β2)
end


"""
Sample energy of the liberated electron in Moller scattering.
Based on the GEANT4 Phys. Ref. Manual v. 11.3 p. 130.
"""
function sample_secondary_energy(m::Moller, eng)
    (;tcut) = m

    ϵ0 = tcut / eng
    mc2 = co.electron_mc2
    γ = 1 + eng / mc2

    β = sqrt((γ^2 - 1) / γ^2)
    y = 1 / (γ + 1)

    local ϵ
    accept = false
    while !accept
        r = rand()
        ϵ = ϵ0 / (1 - r + 2 * ϵ0 * r)

        g = 4 / (9γ^2 - 10γ + 5) * ((γ - 1)^2 * ϵ^2
                                    - (2γ^2 + 2γ - 1) * (ϵ / (1 - ϵ))
                                    + γ^2 / (1 - ϵ)^2)
        accept = rand() < g
    end

    return ϵ * eng
end
