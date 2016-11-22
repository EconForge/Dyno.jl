import Dolang

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
        calibration_grouped[vg] =  [calibration[v] for v in symbols[vg]]
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
    println(code_2)

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

function import_data(model_data)

    symbols_str = model_data["symbols"]

    symbols = OrderedDict()
    for k_str in ["endogenous", "exogenous", "parameters"]
        if k_str in keys(symbols_str)
            symbols[Symbol(k_str)] = [Symbol(x) for x in symbols_str[k_str]]
        end
    end

    equations = [parse_equation(ee) for ee in model_data["equations"]]

    cdict = Dict{Symbol, Union{Expr,Float64, Int64}}(Symbol(k)=>parse(String(v)) for (k,v) in  model_data["calibration"])
    calibration = triangular_system(cdict)

    p = length(model_data["exogenous"]["sigma"])
    sigma = zeros(p,p)
    for i=1:size(sigma,1)
        for j=1:size(sigma,2)
            sigma[i,j] = eval(parse(model_data["exogenous"]["sigma"][i,j]))
        end
    end

    exogenous = MultivariateNormal(zeros(size(sigma,1)), sigma)
    model = Model(symbols, equations, calibration, exogenous)
    return model

end


function import_modfile(filename)
    model_data = modfile_parser(filename)
    model = import_data(model_data)
    return model
end
