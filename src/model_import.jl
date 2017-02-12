immutable MyType end

type Functions
    f_static
    f_dynamic
end

type Model
    symbols:: OrderedDict{Symbol,Array{Symbol,1}}
    equations:: Array{Expr,1}
    calibration:: Dict{Symbol, Float64}
    calibration_grouped:: Dict{Symbol, Array{Float64,1}}
    exogenous # type to be defined
    functions:: Functions
end

function Model(symbols, equations, calibration, exogenous)

    calibration_grouped = OrderedDict{Symbol,Array{Float64,1}}()
    for vg in [:endogenous,:exogenous,:parameters]
        if vg == :parameters
            default = NaN
        else
            default = 0.0
        end
        calibration_grouped[vg] =  [get(calibration,v,default) for v in symbols[vg]]
    end

    # create functions

    sym_vars = symbols[:endogenous]
    sym_parms = symbols[:parameters]
    sym_exo = symbols[:exogenous]

    v_args = cat(1, sym_vars, sym_exo)


    # analyze equations
    it = Dolang.IncidenceTable(equations)

    # we list the variables of the dynamic equations
    # only variables which appear symbolically are added to that vector
    v_args = cat(1,
                [(v,1) for v in sym_vars if v in it.by_date[1]],
                [(v,0) for v in sym_vars if v in it.by_date[0]],
                [(v,-1) for v in sym_vars if v in it.by_date[-1]],
                [(v,0) for v in sym_exo if v in it.by_date[0]]
            )
    ss_args = cat(1,
                [(v,0) for v in sym_vars],
                [(v,0) for v in sym_exo]
            )
    p_args = sym_parms


    # BUG: compatibility fix: exp is converted to E by SymEngine
    E = e

    # steady-state equations
    ss_equations = [Dolang.steady_state(eq) for eq in equations]

    # create functions
    code_1 = make_method(ss_equations, ss_args, p_args, funname=:f_s, orders=[0,1])
    fun_temp_s = eval(code_1)
    code_2 = make_method(equations, v_args, p_args, funname=:f_d, orders=[0,1])
    fun_temp_d = eval(code_2)

    functions = Functions(fun_temp_s,fun_temp_d)

    m = Model(symbols, equations, calibration, calibration_grouped, exogenous, functions)

end


function parse_equation(eq_string)
    if '=' in eq_string
        a, b = split(eq_string,'=')
        lhs, rhs = parse(a), parse(b)
        return :($rhs - $lhs)
    else
        return parse(eq_string)
    end
end


function read_modfile(modfile)
    root_dir = Pkg.dir("Dyno")
    dynare_exe = joinpath(root_dir, "bin", "dynare_m")

    out = try
        readstring(`$dynare_exe $modfile json jsonstdout`)
    catch y
        # TODO: find a way to capture output
        try run(`$dynare_exe $modfile json`) end
        error("Modfile parsing failed")
    end

    regex =  r"\{(.*)\}"s
    m = match(regex, out)
    model_data = JSON.parse(m.match)
    return model_data
end

function import_data(model_data)

    # symbols_str = model_data["symbols"]

    symbols = OrderedDict()
    for k_str in ["endogenous", "exogenous", "parameters"]
        if k_str in keys(model_data)
            symbols[Symbol(k_str)] = [Symbol(strip(x["name"])) for x in model_data[k_str]]
        end
    end
    #
    equations = [parse_equation(ee["equation"]) for ee in model_data["model"]]
    #
    #list all declaration statements
    decl_stmts = []
    for stmt in model_data["statements"]
        if stmt["statementName"] == "param_init"
            push!(decl_stmts, stmt)
        elseif stmt["statementName"] == "init_val"
            push!(decl_stmts, stmt["vals"]...)
        end
    end

    cdict = Dict{Symbol, Union{Symbol, Expr, Float64, Int64}}() #Symbol(k)=>parse(String(v)) for (k,v) in  model_data["calibration"])
    for stmt in decl_stmts
        val = stmt["value"]
        if typeof(val) <: String
            cdict[Symbol(stmt["name"])] = parse(val)
        else
            cdict[Symbol(stmt["name"])] = val
        end
    end

    calibration = Dolang.solve_triangular_system(cdict)

    #
    p = length(symbols[:exogenous])
    sigma = zeros(p,p)
    for i in 1:p
        # This is really rediculous ! (to avoid degenerate distribution)
        sigma[i,i]=1e-10
    end

    exogenous = MultivariateNormal(zeros(size(sigma,1)), sigma)

    model = Model(symbols, equations, calibration, exogenous)
    return model

end

function import_model(filename)
    model_data = read_modfile(filename)
    model = import_data(model_data)
    return model
end
