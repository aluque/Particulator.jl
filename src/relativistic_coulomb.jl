#=
Collection of formulas and sampling algos for relativistic screened Coulomb scattering.
=#

struct RelativisticCoulomb{T} <: CollisionProcess
    "Target atomic number"
    Z::T
end

function collide(c::RelativisticCoulomb, electron::ElectronState{T}, energy) where T
    ϕ = 2π * rand()
    β = norm(electron.v) / co.c

    a = 1.3413 * c.Z^(-1//3) * a_0
    p = momentum(electron)

    # Screening factor in the denominator, do not confuse with fine structure α
    α = co.hbar^2 / (4 * dot(p, p) * a^2)

    cosθ = sample_rel_sr(α, β)
    vnew = turn(electron.v, cosθ, ϕ)
    StateChangeOutcome(ElectronState{T}(electron.x, vnew, electron.w, electron.s, electron.t,
                                        electron.active))
end


"""
Sample azimuthal angles for a cross section proportional to

(1 - β^2 sin^2(θ/2)) / (sin^2(θ/2) + α)^2.

The strategy is to (1) change variable to x = sin^2(θ/2) and (2) decompose the resulting pdf into
an easy to sample envelope and a remaining part that is accounted for by rejection sampling. For high
energies the sample is almost never rejected.
"""
function sample_rel_sr(α, β)
    local x
    accept = false
    while !accept
        u = rand()
        x = α * u / (α + 1 - u)

        # Test if accept
        z = rand()
        accept = (z < (1 - β^2 * x))
    end

    # We have
    # θ = 2 * asin(sqrt(x))
    # but as an optimization we actually return cos(θ), which, as x = sin^2(θ/2) is
    return 1 - 2x
end


"""
Total elastic scattering cross-section at (kinteic) energy `K`
for a single species.
"""
function totalcs(C::RelativisticCoulomb, K)
    (;Z) = C

    # Avoid NaNs when K is exactly 0
    K += 1e-4 * co.eV
    
    a = 1.3413 * Z^(-1//3) * co.a_0
    γ = 1 + K / (co.electron_mass * co.c^2)
    p = sqrt(K * (K + 2 * (co.electron_mass * co.c^2))) / co.c
    β = p / (γ * co.electron_mass * co.c)

    # Screening factor in the denominator, do not confuse with fine structure α
    α = co.hbar^2 / (4 * p^2 * a^2)
    A = (Z * co.r_e / (2 * β^2 * γ))^2

    return (π * r_e^2 * Z^2 / (β^4 * γ^2) *
            ((1 + α * β^2) / (α * (1 + α)) + β^2 * log(α / (1 + α))))
end

