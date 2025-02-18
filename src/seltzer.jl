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

    function SeltzerBerger(d::RawG4Physics2DVector{T}, Z;
                           ncum=1000,
                           pcum=LinRange(0.0, 1.0, ncum),
                           gamma_min=1e2 * co.eV,
                           gamma_max=5e7 * co.eV,
                           energy_scale=1e6 * co.eV) where {T}
        # primary energies
        log_energy = log.(exp.(d.y) * energy_scale)
        data = zeros(T, (length(pcum), length(log_energy)))
        mc2 = co.electron_mc2
        totalcs = zeros(T, length(log_energy))
        
        for i in eachindex(log_energy)
            # d.x = fraction of prim. energy taken by sec.
            T1 = exp(log_energy[i])
            k = d.x
            β = sqrt(1 - 1 / (1 + T1 / mc2)^2)

            logk = log.(k)
            data[:, i] = findcumvalues(logk, d.value[i, :], pcum,
                                       log(gamma_min / T1), log(gamma_max / T1))
            totalcs[i] = ((Z^2 / β^2) * 1e-31 *
                scaledcs(logk, d.value[i, :], log(gamma_min), log(gamma_max)))
        end
        
        new{T, typeof(pcum)}(log_energy, totalcs, pcum, data)
    end

    function SeltzerBerger(T::Type, fname::AbstractString, Z; kw...)
        d = RawG4Physics2DVector{T}(fname)
        SeltzerBerger(d, Z; kw...)
    end

    function SeltzerBerger(T::Type, Z::Int; kw...)
        fname = joinpath(DATA_DIR, "brem_SB", "br$(Z)")
        
        SeltzerBerger(T, fname, Z; kw...)
    end

    SeltzerBerger(fname::AbstractString, args...; kw...) = SeltzerBerger(Float64, fname, args...; kw...)
    SeltzerBerger(Z::Int, args...; kw...) = SeltzerBerger(Float64, Z, args...; kw...)
end


function collide(sb::SeltzerBerger, electron::ElectronState{T}, eng) where T
    mc2 = co.electron_mc2
    m2c2 = co.electron_mass^2 * co.c^2
    m2c4 = mc2^2
    
    k = eng * sample_secondary_energy(sb, eng)
    @assert k < eng "Secondary energy surpasses primary: $(k / co.eV) eV ~ $(eng / co.eV) eV"
    
    # Photon momentum magnitude
    photon_p_norm = k / co.c
        
    cosθ = sample_cos_theta(sb, eng)
    ϕ = 2π * rand()
    
    # Photon momentum
    p_ph = turn(electron.v, cosθ, ϕ, photon_p_norm)

    # Electron momentum
    p0 = momentum(electron)
    p_e = p0 - p_ph
    γ2 = (1 + dot(p_e, p_e) / m2c2)
    β = sqrt(1 - 1 / γ2)
    v_e = co.c * β * p_e / norm(p_e)
    
    NewParticleOutcome(ElectronState{T}(electron.x, v_e, electron.w, electron.t),
                       PhotonState{T}(electron.x, p_ph, electron.w, electron.t))
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
Note: currently the SeltzerBerger data or the secondary energy are not used.
"""
sample_cos_theta(::SeltzerBerger, T) = sample_modified_tsai_cos_theta(T)


"""
Find the values of `x` for which the cumulative probability of a given (unnormalized) probability
density `p` reach the values in `pcum`. Store the results in `s`
"""
function findcumvalues!(s, x, p, pcum, xmin, xmax; rtol=1e-6)
    itp = BSplineKit.interpolate(x, p, BSplineOrder(2))
    cumint = integral(itp)
    @assert all(itp.(BSplineKit.knots(itp)) .> 0)

    # Data points to normalize the comumative distribution as 0 for kmin, 1 for kmax
    cum0 = cumint(max(xmin, minimum(x)))
    cum1 = cumint(min(xmax, maximum(x)))

    deriv = diff(cumint)
    knt = BSplineKit.knots(cumint)
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
    itp = BSplineKit.interpolate(logk, s, BSplineOrder(2))
    cumint = integral(itp)
    @assert all(itp.(BSplineKit.knots(itp)) .> 0)

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
    
    if (i == 0)
        return zero(logK)
    end
    
    w = (sb.log_energy[i + 1] - logK) / (sb.log_energy[i + 1] - sb.log_energy[i])
    
    return w * sb.totalcs[i] + (1 - w) * sb.totalcs[i + 1]
end
