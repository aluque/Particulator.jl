#=
Collection of formulas and sampling algos for RBEB scattering.
=#

struct RBEB{T} <: CollisionProcess
    "Bounding energy of the orbital"
    B::T

    "Mean kinetic energy of the orbital"
    U::T

    "Electron occupation number"
    N::Int
end

function Base.show(io::IO, r::RBEB)
    uev = r.U / co.eV
    bev = r.B / co.eV
    
    print(io, "RBEB(B = $bev eV, U = $uev eV, N = $(r.N))")
end

# B and U of molecular orbitals,
# Hwang at al. J. Chem. Phys. 104, 2956 (1996)
N2_ORBITALS = [
    # B (eV)   U (eV)   N
    # Santos J. Phys. B: At. Mol. Opt. Phys. 36 (2003) 4211–4224 for K-Shell orbitals
    # 4 electrons (2 in each atomic 1s orbital)
    RBEB(0.4095e3 * co.eV,   0.6033e3 * co.eV, 4),
    
    # Hwang at al. J. Chem. Phys. 104, 2956 (1996) for molecular orbitals
    RBEB(41.72 * co.eV,    71.13 * co.eV,    2),
    RBEB(21.00 * co.eV,    63.18 * co.eV,   2),
    RBEB(17.07 * co.eV,    44.30 * co.eV,   4),
    RBEB(15.58 * co.eV,    54.91 * co.eV,   2)]

O2_ORBITALS = [
    # B (eV)   U (eV)   N
    # Santos J. Phys. B: At. Mol. Opt. Phys. 36 (2003) 4211–4224 for K-Shell orbitals
    # 4 electrons (2 in each atomic 1s orbital)
    RBEB(0.5438e3 * co.eV,   0.7962e3 * co.eV, 4),

    # Hwang at al. J. Chem. Phys. 104, 2956 (1996) for molecular orbitals
    RBEB(46.19 * co.eV,    79.73 * co.eV,   2),
    RBEB(29.82 * co.eV,    90.92 * co.eV,   2),
    RBEB(19.64 * co.eV,    59.89 * co.eV,   4),
    RBEB(19.79 * co.eV,    71.84 * co.eV,   2),
    RBEB(12.07 * co.eV,    84.88 * co.eV,   2)]


const ORBITALS = Dict("N2" => N2_ORBITALS,
                      "O2" => O2_ORBITALS)

function collide(c::RBEB, electron::ElectronState{T}, eng) where T
    (;B) = c
    
    mc2 = co.electron_mc2

    E0 = eng
    E2 = rbeb_sample(E0, B)
    E1 = E0 - E2 - B

    @assert E2 < E1
    
    p0 = sqrt(E0^2 + 2 * mc2 * E0) / co.c
    p1 = sqrt(E1^2 + 2 * mc2 * E1) / co.c
    p2 = sqrt(E2^2 + 2 * mc2 * E2) / co.c

    γ1 = 1 + (E1 / mc2)
    γ2 = 1 + (E2 / mc2)

    vnorm = norm(electron.v)
    v1norm = p1 / γ1 / co.electron_mass
    v2norm = p2 / γ2 / co.electron_mass
    
    # cosθ1 = (p0^2 + p1^2 - p2^2) / (2 * p0 * p1)
    # cosθ2 = (p0^2 + p2^2 - p1^2) / (2 * p0 * p2)

    # Lehtinen 1999
    cosθ1 = sqrt(E1 * (E0 + 2mc2) / (E0 * (E1 + 2mc2)))    
    cosθ2 = sqrt(E2 * (E0 + 2mc2) / (E0 * (E2 + 2mc2)))
    
    ϕ = 2π * rand()

    v1 = turn(electron.v, cosθ1, ϕ, v1norm)
    v2 = turn(electron.v, cosθ2, -ϕ, v2norm)

    NewParticleOutcome(ElectronState{T}(electron.x, v1, electron.w, electron.t),
                       ElectronState{T}(electron.x, v2, electron.w, electron.t))
end


"""
Differential cross section for RBEB scattering with secondary kinetic energy `W`, primary
kinetic energy `T`, bounding energy of the orbital `B` and mean kinetic energy of the orbital `U`.
"""
function rbeb_dσdW(W, T, B, U)
    mc2 = co.electron_mc2
    t1 = T / mc2
    b1 = B / mc2
    u1 = U / mc2

    βt2 = 1 - 1 / (1 + t1)^2
    βb2 = 1 - 1 / (1 + b1)^2
    βu2 = 1 - 1 / (1 + u1)^2

    t = T / B
    w = W / B

    _rbeb_dσdW(w, b1, t, t1, βb2, βt2, βu2)
end

function _rbeb_dσdW(w, b1, t, t1, βb2, βt2, βu2)
    α = co.fine_structure

    # Note that βs are always squared
    return ((2 * co.a_0^2 * π * α^4 / (b1 * (βb2 + βt2 + βu2))) * 
        (b1^2 / (1 + t1)^2 + 1 / (t - w)^2 + 1 / (1 + w)^2 -
        ((1 + 2 * t1) * (1 / (t - w) + 1 / (1 + w))) / ((1 + t) * (1 + t1)^2) + 
        (1 / (t - w)^3 + 1 / (1 + w)^3) * (-βt2 - log(2 * b1) + log(βt2 / (1 - βt2)))))
end


"""
Total ionization cross section for RBEB scattering of an orbital `orb`.
"""
function totalcs(C::RBEB, T)
    (;B, U, N) = C
    return N * rbeb_σI(T, B, U) * (T > B)
end


"""
Total ionization cross section for RBEB scattering with primary kinetic energy `T`,
bounding energy of the orbital `B` and mean kinetic energy of the orbital `U`.
"""
function rbeb_σI(T, B, U)
    mc2 = co.electron_mc2
    t1 = T / mc2
    b1 = B / mc2
    u1 = U / mc2

    βt2 = 1 - 1 / (1 + t1)^2
    βb2 = 1 - 1 / (1 + b1)^2
    βu2 = 1 - 1 / (1 + u1)^2

    t = T / B

    _rbeb_σI(b1, t, t1, βb2, βt2, βu2)
end

function _rbeb_σI(b1, t, t1, βb2, βt2, βu2)
    α = co.fine_structure

    return (2 * co.a_0^2 * π * α^4 / (b1 * (βb2 + βt2 + βu2)) *
        (1 - 1 / t +
        (b1^2 * (t - 1))/ (2 * (1 + t1)^2) -
        ((1 + 2 * t1) * log(t)) / ((1 + t) * (1 + t1)^2) -
        ((t^2 - 1)*(βt2 - log(βt2) + log(2 * b1 * (1 - βt2))))/(2 * t^2)))
end


"""
Sample the secondary energy in an RBEB collision with primary energy `T` and bounding energy `B`.
"""
function rbeb_sample(T, B)
    mc2 = co.electron_mass * co.c^2
    t1 = T / mc2
    b1 = B / mc2

    βt2 = 1 - 1 / (1 + t1)^2

    t = T / B

    w = B * _rbeb_sample(b1, t, t1, βt2)
    return w
end


function _rbeb_sample(b1, t, t1, βt2)
    # We use rejection sampling. We bound the probability from above by
    # (2B + 2C + (t + 1)^2 * B * M / 4) / (1 + w)^2
    
    A = -(1 + 2 * t1) / (t + 1) / (1 + t1)^2
    B = 1
    C = (log(βt2 / (1 - βt2)) - βt2 - log(2 * b1))
    M = b1^2 / (1 + t1)^2

    accept = false

    g(t, w, q) = 1 / (w + 1)^q + 1 / (t - w)^q
    local w
    while !accept
        # First, sample a distribution ~1/(1+w)^2
        u = rand()
        w = u / ((t + 1) / (t - 1) - u)
        
        # Check if we accept
        pb = (2B + 2C + (t + 1)^2 * B * M / 4) / (1 + w)^2        
        p0 = A * g(t, w, 1) + B * (g(t, w, 2) + M) + C * g(t, w, 3)

        accept = rand() * pb < p0
    end

    return w
end
