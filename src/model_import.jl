
function read_modfile(modfile)

    
    dynare_exe = joinpath(CondaPkg.bindirs()[1], "dynare-preprocessor")

    # root_dir = Pkg.dir("Dyno")
    # dynare_exe = joinpath(root_dir, "deps", "usr", "bin", "dynare_m")

    out = try
        read(`$dynare_exe $modfile json=parse jsonstdout`, String)
    catch y
        # TODO: find a way to capture output
        error("Modfile parsing failed")
    end

    regex =  r"\{(.*)\}"s
    m = match(regex, out)
    model_data = JSON.parse(m.match)
    return model_data
end

function parse_equation(eq::Dict{String, Any}; variables=Symbol[])
    lhs = Dolang.parse_string(eq["lhs"]; variables=variables)
    rhs = Dolang.parse_string(eq["rhs"]; variables=variables)
    return :($rhs-$lhs)
end

struct Model{ID}
    symbols:: OrderedDict{Symbol,Array{Symbol,1}}
    equations:: Array{Expr,1}
    calibration:: Dict{Symbol, Float64}
    calibration_grouped:: Dict{Symbol, Array{Float64,1}}
    exogenous # type to be defined
end

function Model(symbols, equations, calibration, exogenous; print_code=false)

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

    # steady-state equations

    ss_equations = OrderedDict(
        (Symbol("eq_$i")=>Dolang.stringify(Dolang.steady_state(eq)) for (i,eq) 
        in enumerate(equations)
        )...
    )

    ss_arguments = OrderedDict(
        :y => [Dolang.stringify((v,0)) for v in sym_vars],
        :e => [Dolang.stringify((v,0)) for v in sym_exo],
        :p => [Dolang.stringify(e) for e in sym_parms]
    )

    d_equations = OrderedDict(
        (Symbol("eq_$i")=>Dolang.stringify(eq) for (i,eq) 
        in enumerate(equations)
        )...
    )
    d_arguments = OrderedDict(
        :y_2=>[Dolang.stringify(v,1) for v in sym_vars ],
        :y_1=>[Dolang.stringify(v,0) for v in sym_vars ],
        :y_0=>[Dolang.stringify(v,-1) for v in sym_vars ],
        :e=>[Dolang.stringify(v,0) for v in sym_exo ],
        :p=>[Dolang.stringify(e) for e in sym_parms]
    )

    ff_ss = Dolang.FunctionFactory(ss_equations, ss_arguments, Dict(), :f_s)
    ff_d = Dolang.FunctionFactory(d_equations, d_arguments, Dict(), :f_d)

    
    # code_1 = make_function(ss_equations, ss_args, p_args, funname=:f_s, orders=[0,1])
    # code_1 = Dolang.gen_generated_gufun(ff_ss; funname= :static)

    ID = gensym()

    code_1 = Dolang.gen_kernel2(ff_ss, 0; funname= :static, dispatch=Model{ID})

    if print_code
        println("****************")
        println("Static functions")
        println("****************")
        println(code_1)
    end
    eval(code_1)

    code_2 = Dolang.gen_kernel2(ff_d, 0; funname= :dynamic, dispatch=Model{ID})

    if print_code
        println("*****************")
        println("Dynamic functions")
        println("*****************")
        println(code_2)
    end
    eval(code_2)



    Model{ID}(symbols, equations, calibration, calibration_grouped, exogenous)

end



function import_data(model_data; print_code=true)

    # symbols_str = model_data["symbols"]

    modfile = model_data["modfile"]

    symbols = OrderedDict()
    for k_str in ["endogenous", "exogenous", "parameters"]
        symbols[Symbol(k_str)] = [Symbol(strip(x["name"])) for x in modfile[k_str]]
    end

    variables = cat(symbols[:endogenous], symbols[:exogenous]; dims=1)
    #
    equations = [parse_equation(ee; variables=variables) for ee in modfile["model"]]

    variables = [symbols[:endogenous]; symbols[:exogenous]]
    # equations = [sanitize(eq, variables) for eq in equations]

    #
    #list all declaration statements
    decl_stmts = []
    for stmt in modfile["statements"]
        if stmt["statementName"] == "param_init"
            push!(decl_stmts, stmt)
        elseif stmt["statementName"] == "init_val"
            push!(decl_stmts, stmt["vals"]...)
        end
    end

    # cdict = Dict{Symbol, Union{Symbol, Expr, Float64, Int64}}() #Symbol(k)=>parse(String(v)) for (k,v) in  model_data["calibration"])
    # for stmt in decl_stmts
    #     val = stmt["value"]
    #     if typeof(val) <: String
    #         cdict[Symbol(stmt["name"])] = Dolang.parse_string(val)
    #     else
    #         cdict[Symbol(stmt["name"])] = val
    #     end
    # end
    # calibration = Dolang.solve_triangular_system(cdict)

    calibration = Dict{Symbol, Float64}()
    for stmt in decl_stmts
        val_s = stmt["value"]
        sym = Symbol(stmt["name"])
        val = Dolang.eval_string(val_s, calibration)
        calibration[sym] = val
    end

    # return calibration
    #
    p = length(symbols[:exogenous])
    sigma = zeros(p,p)
    for i in 1:p
        # This is really ridiculous ! (to avoid degenerate distribution)
        sigma[i,i]=1e-20
    end
    for stmt in modfile["statements"]
        if stmt["statementName"] == "shocks"
            for stmt_var in stmt["variance"]
                exo_name = Symbol(stmt_var["name"])
                exo_i = findfirst(isequal(exo_name),symbols[:exogenous])
                val_s = string(stmt_var["variance"])
                val = Dolang.eval_string(val_s, calibration)
                sigma[exo_i, exo_i] = val
            end
            for stmt_var in stmt["stderr"]
                exo_name = Symbol(stmt_var["name"])
                # @show symbols
                # @show exo_name
                exo_i = findfirst(isequal(exo_name),symbols[:exogenous])
                val_s = string(stmt_var["stderr"])
                
                val = Dolang.eval_string(val_s, calibration)
                sigma[exo_i, exo_i] = val
            end


            for stmt_cvar in stmt["covariance"]
                exo_name1 = Symbol(stmt_cvar["name"])
                exo_name2 = Symbol(stmt_cvar["name2"])
                exo_i1 = findfirst(isequal(exo_name1), symbols[:exogenous])
                exo_i2 = findfirst(isequal(exo_name2), symbols[:exogenous])
                val_s = string(stmt_cvar["covariance"])
                val = Dolang.eval_string(val_s, calibration)
                sigma[exo_i1, exo_i2] = val
                sigma[exo_i2, exo_i1] = val
            end

            #     throw("Unsupported 'shocks' statement (for now).")
            # end
        end
    end



    exogenous = MultivariateNormal(zeros(size(sigma,1)), sigma)

    model = Model(symbols, equations, calibration, exogenous; print_code=print_code)
    return model

end

function import_modfile(filename; print_code=false)
    model_data = read_modfile(filename)
    model = import_data(model_data; print_code=print_code)
    return model
end


# function Model(filename; print_code=false)
#     return import_modfile(filename; print_code=print_code)
# end
