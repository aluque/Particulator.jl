#=
Collection of formulas and sampling algos for bremsstrahlung using the Seltzer-Berger data
and Geant4 scattering models.
=#

"""
Struct to sample Bremmstrahlung data from Seltzer and Berger data tables.
"""
struct SeltzerBerger{T, V <: AbstractVector{T}}
    "log of the electron (primary) energies contained in the tables"
    log_energy::Vector{T}

    "Total cross section for each of the energy values"
    totalcs::Vector{T}
    
    "pre-computed cumulative probabilities"
    pcum::V

    "Data with the photon (secondary) log-energy corresponding to a cum. probability for a
    primary log-energy"
    data::Matrix{T}

    function SeltzerBerger(d::RawG4Physics2DVector{T}, pcum::V, Z;
                           gamma_min=1e2 * co.eV, gamma_max=5e7 * co.eV,
                           energy_scale=1e6 * co.eV) where {T, V}
        # primary energies
        log_energy = log.(exp.(d.y) * energy_scale)
        data = zeros(T, (length(pcum), length(log_energy)))
        mc2 = co.electron_mass * co.c^2
        totalcs = zeros(T, length(log_energy))
        
        for i in eachindex(log_energy)
            # d.x = fraction of prim. energy taken by sec.
            T1 = exp(log_energy[i])
            k = T1 .* d.x
            β = sqrt(1 - 1 / (1 + T1 / mc2)^2)

            logk = log.(k)
            data[:, i] = findcumvalues(logk, d.value[i, :], pcum, log(gamma_min), log(gamma_max))
            totalcs[i] = ((β^2 / Z^2) * 1e-31 *
                scaledcs(logk, d.value[i, :], log(gamma_min), log(gamma_max)))
        end
        
        new{T, V}(log_energy, totalcs, pcum, data)
    end

    function SeltzerBerger(T, fname::AbstractString, pcum::V, Z; kw...) where V
        d = Seltzer.RawG4Physics2DVector{T}(fname)
        SeltzerBerger(d, pcum, Z; kw...)
    end
end


function collide(c::SeltzerBerger, electron::ElectronState{T}, eng) where T
    mc2 = co.electron_mass * co.c^2
    m2c4 = mc2^2
    
    k = sample_secondary_energy(sb, eng)
    # Photom momentum magnitude
    photon_p_norm = k / co.c
        
    cosθ = sample_photon_cos_theta(sb, eng)
    ϕ = 2π * rand()
    
    # Photon momentum
    p_ph = turn(electron.v, cosθ, ϕ, photon_p_norm)

    # Electron momentum
    p0 = momentum(electron)
    p_e = p0 - p_ph
    γ2 = (1 + dot(p_e, p_e) / m2c4)
    β = sqrt(1 - 1 / γ2)
    v_e = co.c * β * p_e / norm(p_e)
    
    NewParticleOutcome(ElectronState{T}(electron.x, v_e, electron.w, electron.s, electron.active),
                       PhotonState{T}(electron.x, p_ph, electron.w, electron.s, electron.active))
end



"""
Sample energy of secundary from primary kinetic energy `T` using the data in `sb`.
"""
function sample_secondary_energy(sb::SeltzerBerger, T)
    x = rand()
    y = log(T)
    u = sb.data
    
    i2 = searchsortedfirst(sb.pcum, x)
    i1 = i2 - 1
        
    j2 = searchsortedfirst(sb.log_energy, y)
    j1 = j2 - 1
    
    # Bi-linear interp
    x1 = sb.pcum[i1]
    x2 = sb.pcum[i2]
    y1 = sb.log_energy[j1]
    y2 = sb.log_energy[j2]

    
    A = (x2 - x1) * (y2 - y1)
    S = (u[i1, j1] * (x2 - x) * (y2 - y) +
        u[i1, j2] * (x - x1) * (y2 - y) +
        u[i2, j1] * (x2 - x) * (y - y1) +
        u[i2, j2] * (x - x1) * (y - y1))

    return exp(S / A)
end

"""
Sample the azimuthal angle θ given the primary kinetic energy `T`.
This method is obtained from the GEANT4 physics reference manual, (release 10.4 p. 78).

Note: currently the SeltzerBerger data or the secondary energy are not used.
"""
function sample_photon_azimuth(::SeltzerBerger, T)
    a = 0.625
    d = 27
    accept = false
    factor = (1 + T / (co.electron_mass * co.c^2))
    
    local u
    while !accept
        r1 = rand()
        b = r1 < 9 / (9 + d) ? a : 3a
        
        r2 = rand()
        r3 = rand()
        u = -log(r2 * r3) / b

        accept = u <= factor * π
    end

    return u / factor
end


"""
Sample the cosine of azimuthal angle θ given the primary kinetic energy `T`.
This method is obtained from the GEANT4 source code:

source/processes/electromagnetic/standard/src/G4ModifiedTsai.cc
(commit c07cea1fe028470cd9050371f165c7c815eefb23 in github)

The resulting theta distribution is the same as for sample_azimuth for small theta but differs
for larger thetas (they are the same if cos θ = 1 - θ^2/2).

Note: currently the SeltzerBerger data or the secondary energy are not used.
"""
function sample_photon_cos_theta(::SeltzerBerger, T)
    umax = 2 * (1 + T / (co.electron_mass * co.c^2))
    a1 = 1.6
    a2 = a1 / 3
    border = 0.25
    
    accept = false
    
    local u
    while !accept
        uu = -log(rand() * rand())
        u = border > rand() ? uu * a1 : uu * a2
        accept = u <= umax
    end

    return 1 - 2 * u^2 / umax^2
end


"""
Find the values of `x` for which the cumulative probability of a given (unnormalized) probability
density `p` reach the values in `pcum`. Store the results in `s`
"""
function findcumvalues!(s, x, p, pcum, xmin, xmax; rtol=1e-6)
    itp = interpolate(x, p, BSplineOrder(2))
    cumint = integral(itp)
    @assert all(itp.(knots(itp)) .> 0)

    # Data points to normalize the comumative distribution as 0 for kmin, 1 for kmax
    cum0 = cumint(max(xmin, minimum(x)))
    cum1 = cumint(min(xmax, maximum(x)))

    deriv = diff(cumint)
    knt = knots(cumint)
    fknt = @. ((cumint(knt) - cum0) / (cum1 - cum0))

    for i in eachindex(pcum)
        j = searchsortedlast(fknt, pcum[i])

        xsol = 0.5 * (knt[max(firstindex(knt), j)] + knt[min(lastindex(knt), j + 1)])
        dx = Inf
        while abs(dx / xsol) > rtol
            f = (cumint(xsol) - cum0) / (cum1 - cum0)
            df = deriv(xsol) / (cum1 - cum0)
            dx = (f - pcum[i]) / df 
            xsol = xsol - dx
        end
        
        s[i] = xsol
    end

    return s
end

findcumvalues(k, p, pcum, kmin, kmax; rtol=1e-6) = 
    findcumvalues!(Vector{eltype(pcum)}(undef, length(pcum)), k, p, pcum, kmin, kmax; rtol)

"""
Compute the total scaled cross section from tabulated `logk` values and `s` = k * dσ/dk with integration
interval (`kmax`, `kmin`).
"""
function scaledcs(logk, s, logkmin, logkmax)
    itp = interpolate(logk, s, BSplineOrder(2))
    cumint = integral(itp)
    @assert all(itp.(knots(itp)) .> 0)

    # Data points to normalize the comumative distribution as 0 for kmin, 1 for kmax
    cum0 = cumint(max(logkmin, minimum(logk)))
    cum1 = cumint(min(logkmax, maximum(logk)))

    return cum1 - cum0
end


"""
Compute total cross section for bremsstrahlung process modeled by `sb` with primary kinetic energy `K`
"""
function totalcs(sb::SeltzerBerger, K)
    logK = log(K)

    # We allow non-uniform sampling of energies because this method is supposed to be used
    # only for initialization
    i = searchsortedlast(sb.log_energy, logK)
    w = (sb.log_energy[i + 1] - logK) / (sb.log_energy[i + 1] - sb.log_energy[i])
    
    return w * sb.totalcs[i] + (1 - w) * sb.totalcs[i + 1]    
end
