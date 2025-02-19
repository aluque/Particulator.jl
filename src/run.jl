function run!(mpopl, efield, bfield, tfinal, dt, callback; output_dt=tfinal / 20)
    t = 0.0
    nxt = t
    
    while t < tfinal
        advance!(mpopl, efield, bfield, t + dt, callback)
        foreach(droplow!, mpopl)
        t += dt
        onstep(callback, mpopl, t)
        if !isnothing(output_dt) && t >= nxt
            nxt += output_dt
            msg = _msg(mpopl, t, 100 * t / tfinal)
            @info msg
        end
    end
end

function _msg(mpopl, t, perc)
    nstr = join(map(x -> @sprintf("%15s => %-10s", string(first(x)), repr(last(x))), 
                    map(((a, b),) -> a => nparticles(b), Tuple(pairs(mpopl)))), "\n")
    lstr = join(map(x -> @sprintf("%15s => %-10s", string(first(x)), repr(last(x))), 
                    map(((a, b),) -> a => spread(b)[1], Tuple(pairs(mpopl)))), "\n")
    
    return (@sprintf("%.2f", perc) * "% complete\ntime = " * @sprintf("%.3f", (t / 1e-9)) * " ns\n" * 
        "# of particles:\n" * nstr * 
        "\ncentroid locations:\n" * lstr)
end
