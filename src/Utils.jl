__precompile__()

module Utils

    export @idx0
    using MacroTools: postwalk, @capture

    macro idx0(ex)
        return esc(postwalk(ex) do x
            @capture(x, A_[idxs__]) || return x
            return :( $A[ $([ :( $idx .+ 1 ) for idx in idxs ]...) ] )
        end)
    end


    export prettyshow

    function prettyshow(M)
        show(STDOUT, MIME"text/plain"(), M)
        println()
    end


    export pause

    function pause()
        println(STDOUT, "Press ENTER to continue!")
        readline(STDIN)
    end


    export gnuplotexport

    function gnuplotexport(filename, xlabel, ylabel, data; lt=<=, xmap=identity, ymap=identity)
        open(filename, "w") do fp
            println(fp, "# $(xlabel) $(ylabel)")
            K = collect(keys(data))
            sort!(K; lt=lt)
            for k in K
                v = data[k]
                println(fp, "$(xmap(k)) $(ymap(v))")
            end
        end
    end


    export findmappedmax

    function findmappedmax(f, lst)
        lst_ = collect(lst)
        res = map(f, lst_)
        (max_val, max_idx) = findmax(res)
        return (max_val, lst_[max_idx])
    end


    export myround_towardszero

    # Float64 has 52bits for significand; myround_towardszero() restricts resolution to 2^(52-k) values per [ 2^i , 2^(i+1) ) interval
    myround_towardszero(x::Float64, k::UInt) = reinterpret(Float64, reinterpret(UInt64, x) & (0xffffffffffffffff<<k))

end
