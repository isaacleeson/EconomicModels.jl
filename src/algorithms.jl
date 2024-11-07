function repeatedly_while(maxiters, fs...; tol=1e-10)
    _err = 0.0
    for f in fs
        err = Inf
        i = 1
        while i < maxiters && err isa Real && err > tol
            err = f()
            i += 1
        end
        _err = err
    end
    return _err
end
