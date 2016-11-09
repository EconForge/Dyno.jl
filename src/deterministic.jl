using NLsolve
using Dolang
using DataFrames

function newton(f, init; tol=1e-8, maxit=500, backsteps=5, verbose=false)

    line = 0.5.^(0:(backsteps-1))
    x0 = init[:]
    x = copy(x0)
    N = length(x0)
    it = 0
    err = 1.0
    while (it<maxit) & (err>tol)
        it += 1
        ress, jacc = f(reshape(x0,size(init)...))
        res = reshape(ress,N)
        jac = reshape(jacc,N,N)
        err = maximum(abs(ress))
        ul = 0.1
        if err>tol
            dx = -jac\res
            nerr = 2.0
            for l in line
                ul = l
                try
                    x = x0+l*dx
                    nr, junk = f(reshape(x,size(init)...))
                    nerr = maximum(abs(nr))
                catch
                    println("skip")
                end
                if nerr<err
                    break
                end
            end
            if nerr>err
                println("Backtracking didn't work")
            end
        end
        x0 = x
        if verbose
            println(it, ": ", err, ":", ul)
        end
    end
    return reshape(x0,size(init)...)
end


function deterministic(model, exogenous_series; verbose=false, N=20)


    symbols = model.symbols
    variables = symbols[:variables]
    exogenous = symbols[:exogenous]
    parameters = symbols[:parameters]

    it = Dolang.IncidenceTable(model.equations)

    fut_syms = it.by_date[1]
    cur_syms = it.by_date[0]
    past_syms = it.by_date[-1]

    fut_inds = [i for (i,v) in enumerate(variables) if v in it.by_date[1]]
    cur_inds = [i for (i,v) in enumerate(variables) if v in it.by_date[0]]
    past_inds = [i for (i,v) in enumerate(variables) if v in it.by_date[-1]]
    shocks_inds = [i for (i,v) in enumerate(exogenous) if v in it.by_date[0]]
    parms_inds = [i for (i,v) in enumerate(parameters) if v in it.by_date[0]]

    P = length(symbols[:variables])
    indices = fut_inds
    indices = cat(1,indices, P+cur_inds)

    indices = cat(1,indices, 2*P+past_inds)
    findices = cat(1,indices, 3*P+shocks_inds)

    ne = length(model.symbols[:exogenous])

    T = size(exogenous_series,1)
    exogenous = zeros(N, ne)
    exogenous[1:T,:] = exogenous_series
    for t=T+1:N
        exogenous[t,:] = exogenous[T,:]
    end

    f_s = model.functions.f_static
    f_d = model.functions.f_dynamic


    y0 = [model.calibration[e] for e in model.symbols[:variables]]
    e0 = [model.calibration[e] for e in model.symbols[:exogenous]]
    p = [model.calibration[e] for e in model.symbols[:parameters]]



    e_start = exogenous[1,:]
    e_final = exogenous[end,:]

    function fobj_start(x, out)
        out[:] = f_s(Der{0},cat(1, x, e_start) , p)
    end

    function fobj_final(x, out)
        out[:] = f_s(Der{0},cat(1, x, e_final) , p)
    end


    fdiff_start = NLsolve.DifferentiableMultivariateFunction(fobj_start)
    sol_start = nlsolve(fdiff_start,y0)
    y_start = sol_start.zero
    fdiff_final = NLsolve.DifferentiableMultivariateFunction(fobj_final)
    sol_final = nlsolve(fdiff_final,y0)
    y_final = sol_final.zero

    y_init_0 = [(1-l)*y_start + l*y_final for l in linspace(0,1,N)]
    y_init = cat(1, [e' for e in y_init_0]...) :: Array{Float64, 2}

    n_y = length(model.symbols[:variables])
    n_e = length(model.symbols[:exogenous])


    function residual(y)

        a = y[1:end-2,:]   # y{t-1}
        b = y[2:end-1,:]   # y_t
        c = y[3:end,:]     # y_{t+1}
        d = exogenous[1:end-2,:]

        x = cat(2,c,b,a,d)
        xx = x[:,findices]

        N,P = size(y)

        res = zeros(N,P)
        jac0 = zeros(N,P,3*P)
        for n = 2:N-1
            res[n,:] = f_d(Der{0},xx[n-1,:],p)
            jac0[n,:,indices] = f_d(Der{1},xx[n-1,:],p)[:,1:length(indices)]
        end
        jac = zeros(N,P,N,P)
        for n in 2:N-1
            jac[n,:,n,:] += jac0[n,:,(P+1):2*P]  # y_t
            jac[n,:,n+1,:] += jac0[n,:,1:P] # y_{t+1}
            jac[n,:,n-1,:] +=  jac0[n,:,2*P+1:3*P]  # y_{t-1}
        end
        jac[1,:,1,:] = eye(P)
        jac[end,:,end,:] = eye(P)
        res[1,:] = y[1,:]-y_start
        res[end,:] = y[end,:]-y_final
        return res, jac
    end
    a,b = residual(y_init)
    println("Initial: ", maximum(abs(a)))
    sol = newton(residual, y_init, verbose=verbose)

    tvec = collect(0:(size(sol,1)-1))
    sol = cat(2,tvec,sol,exogenous)

    colnames = cat(1, [:t], model.symbols[:variables], model.symbols[:exogenous])
    columns = Dict(colnames[i] => sol[:,i] for i=1:size(sol,2))
    return DataFrame(columns)
end
