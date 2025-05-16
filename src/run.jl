function run!(mpopl, pusher, tfinal, dt, callback; output_dt=tfinal / 20, verbosity=1)
    t = 0.0
    nxt = t
    isave = 0
    
    while t < tfinal
        advance!(mpopl, pusher, t + dt, callback)
        foreach(droplow!, mpopl)
        t += dt
        cont = onstep(callback, mpopl, t)

        if !isnothing(output_dt) && t >= nxt            
            nxt += output_dt
            isave += 1
            if verbosity > 0
                msg = _msg(mpopl, t, 100 * t / tfinal)
                @info msg
            end
            onoutput(callback, mpopl, t, isave)
        end

        if !cont
            if verbosity > 0
                @info "Early stop due to onstep() callback"
            end
            break
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
