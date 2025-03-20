##
## Continuum energy losses for e- and e+
## See GEANT4 Phys. Ref. Manual 11.3, chap. 11.
##

struct ContinuumLoss{T} <: AbstractForcing
    # Electron density in the material
    nel::T

    # Mean excitation energy
    I::T

    # Min energy cut
    Tcut::T
end

function force(cl::ContinuumLoss, s::Union{ElectronState, PositronState})
    f = energy_loss(cl, s)
    return s.p * (-f / norm(s.p))
end

force(cl::ContinuumLoss, s::PhotonState) = zero(s.p)


"""
Compute energy loss (aka friction force) due to ionization below Tcut for electrons and positrons.
"""
function energy_loss(cl::ContinuumLoss, s::ParticleState, eng)
    (;nel, I, Tcut) = cl
    mc2 = co.electron_mc2
    r_e = co.r_e
    
    τ = eng / mc2
    τc = Tcut / mc2
    τmax = taumax(p, τ)
    γ = 1 + τ
    y = 1 / (γ + 1)
    β2 = 1 - 1 / γ^2    
    τup = min(τc, τmax)

    F = _F(p, τ, τup, β2)

    # Compute δ
    x = log(γ^2 * β2) / log(10) / 2
    hνp = co.hbar * co.c * sqrt(4π * nel * r_e)
    C = 1 + 2 * log(I / hνp)
    xa = C / log(10) / 2
    m = 3
    (x0, x1) = x0x1(C)
    a = 2 * log(10) * (xa - x) / (x1 - x0)^m

    if x < x0
        δ = 0.0
    elseif x < x1
        δ = 2 * log(10) * x - C + a * (x1 - x)^m
    else
        δ = 2 * log(10) * x - C
    end

    return (2π * r_e^2 * mc2 * nel / β2) * (log((2 * (γ + 1)) / (I / mc2)^2) + F - δ)
end

energy_loss(cl, s) = energy_loss(cl, s, energy(s))


taumax(::PositronState, τ) = τ
taumax(::ElectronState, τ) = τ / 2

function _F(::PositronState, τ, τup, β2)
    y = 1 / (2 + τ)
    return (log(τ * τup)
            - (τup^2 / τ) * (τ * 2τup - 3τup^2 * y / 2 - (τup - τup^3 / 3) * y^2
                             - (τup^2 / 2 - τ * τup^3 / 3 + τup^4 / 4) * y^3))
end

function _F(::ElectronState, τ, τup, β2)
    γ = 1 + τ
    return  (-1 - β2 + log((τ - τup) * τup) + τ / (τ - τup)
             + (τup^2 / 2 + (2τ + 1) * log(1 - τup / τ)) / γ^2)
end

function x0x1(C)
    if C < 10
        return (1.6, 4.0)
    elseif C < 10.5
        return (1.7, 4.0)
    elseif C < 11.0
        return (1.8, 4.0)
    elseif C < 11.5
        return (1.9, 4.0)
    elseif C < 12.25
        return (2.0, 4.0)
    elseif C < 13.804
        return (2.0, 5.0)
    else
        return (0.326 * C - 2.5, 5.0)
    end
end
        
