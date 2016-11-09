
type LinearSolution
    names::Dict{Symbol, Array{Symbol, 1}}
    g_0::Array{Float64,1}
    g_y::Array{Float64,2}
    g_e::Array{Float64,2}
    sigma::Array{Float64,2}
end

function solve(model)

    symbols = model.symbols

    n_v = length(model.symbols[:variables])
    n_e = length(model.symbols[:exogenous])

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

    y_ss = model.calibration_grouped[:variables]
    e_ss = model.calibration_grouped[:exogenous]
    p_ss = model.calibration_grouped[:parameters]

    v = cat(1, y_ss[fut_inds], y_ss[cur_inds], y_ss[past_inds], e_ss[shocks_inds])
    p = p_ss[parms_inds]



    f = model.functions.f_dynamic

    # res = f(Der{0}, v,p)
    # check residuals are zero

    jac = f(Der{1}, v,p)


    n_v = length(y_ss)
    n_e = length(e_ss)

    F = zeros(n_v, n_v)
    G = zeros(n_v, n_v)
    H = zeros(n_v, n_v)

    K = zeros(n_v, n_e)

    n1 = length(fut_inds)
    n2 = length(cur_inds)
    n3 = length(past_inds)

    F[:,fut_inds] = jac[:,1:n1]
    G[:,cur_inds] = jac[:,n1+1:n1+n2]
    H[:,past_inds] = jac[:,n1+n2+1:n1+n2+n3]
    K[:,shocks_inds] = jac[:,n1+n2+n3+1:end]


    g_y = solve_second_order(F,G,H)
    g_e = -(F*g_y+G)\K

    names = Dict(
        :variables => model.symbols[:variables],
        :exogenous => model.symbols[:exogenous],
    )

    g_0 = model.calibration_grouped[:variables]
    sigma = model.exogenous.Σ.mat

    return LinearSolution(names, g_0, g_y, g_e, sigma)

end


function solve_second_order(F,G,H)

    n_v = size(F,1)
    N = size(F,1)

    D = [ zeros(N,N) F ;
                eye(N) zeros(N,N)]

    E = [ -H           -G ;
          zeros(N, N) eye(N)] ;

    sf = schurfact(D, E)
    crit = 1.0
    select = abs(sf.alpha).>crit*sf.beta
    ordschur!(sf, select)

    T = sf.S
    S = sf.T

    eigvals = diag(sf.S)./diag(sf.T)

    Z = sf.Z'
    Q = sf.Q

    Z_21 = Z[n_v+1:end, 1:n_v]
    Z_22 = Z[n_v+1:end, n_v+1:end]

    g_1 = - Z_22 \ Z_21
    g_1 = real(g_1)

    return g_1
end


function simulate(sol::LinearSolution; N=20)
    n_v = size(sol.g_0,1)
    n_e = size(sol.sigma, 1)
    y_vec = zeros(n_v, (N+1))
    e_vec = zeros(n_e, (N+1))
    rnrm = MultivariateNormal(zeros(n_e), sol.sigma)
    A = sol.g_y
    B = sol.g_e
    y = sol.g_0
    y_0 = sol.g_0
    y_vec[:,1] = y
    e_vec[:,2:end] = rand(rnrm,N)
    for t=2:N+1
        y = y_0 + A*(y-y_0) + B*e_vec[:,t]
        y_vec[:,t] = y
    end
    res = [y_vec; e_vec]

    tvec = collect(1:size(res,2)) - 1
    res = [tvec';res]

    colnames = cat(1, [:t], sol.names[:variables], sol.names[:exogenous])
    columns = Dict(colnames[i] => res[i,:] for i=1:size(res,1))

    return DataFrame(columns)

end

function simulate(model::Model; N=20)
    sol = solve(model)
    return simulate(sol, N=N)
end



#### higher orders

#
# using Dolang
# import Dyno
# solve_second_order = Dyno.solve_second_order
#
# symbols = model.symbols
#
# n_v = length(model.symbols[:variables])
# n_e = length(model.symbols[:exogenous])
#
# variables = symbols[:variables]
# exogenous = symbols[:exogenous]
# parameters = symbols[:parameters]
#
# it = Dolang.IncidenceTable(model.equations)
#
# fut_syms = it.by_date[1]
# cur_syms = it.by_date[0]
# past_syms = it.by_date[-1]
#
# fut_inds = [i for (i,v) in enumerate(variables) if v in it.by_date[1]]
# cur_inds = [i for (i,v) in enumerate(variables) if v in it.by_date[0]]
# past_inds = [i for (i,v) in enumerate(variables) if v in it.by_date[-1]]
# shocks_inds = [i for (i,v) in enumerate(exogenous) if v in it.by_date[0]]
# parms_inds = [i for (i,v) in enumerate(parameters) if v in it.by_date[0]]
#
# y_ss = model.calibration_grouped[:variables]
# e_ss = model.calibration_grouped[:exogenous]
# p_ss = model.calibration_grouped[:parameters]
#
# v = cat(1, y_ss[fut_inds], y_ss[cur_inds], y_ss[past_inds], e_ss[shocks_inds])
# p = p_ss[parms_inds]
#
#
#
# f = model.functions.f_dynamic
#
# # res = f(Der{0}, v,p)
# # check residuals are zero
#
# jac = f(Der{1}, v,p)
#
#
# n_v = length(y_ss)
# n_e = length(e_ss)
#
# F = zeros(n_v, n_v)
# G = zeros(n_v, n_v)
# H = zeros(n_v, n_v)
#
# K = zeros(n_v, n_e)
#
# n1 = length(fut_inds)
# n2 = length(cur_inds)
# n3 = length(past_inds)
#
# F[:,fut_inds] = jac[:,1:n1]
# G[:,cur_inds] = jac[:,n1+1:n1+n2]
# H[:,past_inds] = jac[:,n1+n2+1:n1+n2+n3]
# K[:,shocks_inds] = jac[:,n1+n2+n3+1:end]
#
#
# g_y = solve_second_order(F,G,H)
# g_e = -(F*g_y+G)\K
#
# names = Dict(
#     :variables => model.symbols[:variables],
#     :exogenous => model.symbols[:exogenous],
# )
#
# g_0 = model.calibration_grouped[:variables]
# sigma = model.exogenous.Σ.mat
#
# hes_sparse = model.functions.f_static(Der{2}, v,p)
# hes = full(hes_sparse)
# nn = n1+n2+n3+n_e
# hes = reshape(hes, size(hes,1), nn, nn)
