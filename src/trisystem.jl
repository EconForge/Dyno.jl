
function triangular_system(dict::Dict{Symbol,Union{Expr,Float64, Int64}})
    context = Dict() # context
    solutions = Dict{Symbol, Float64}()
    finished = false
    N = length(dict)
    n = 0
#     return
    while ~(  (finished) || n>N )
        done_smthg = false
        n += 1
        for k in keys(dict)
            if ~(in(k,keys(solutions)))
                expr = dict[k]
                try
                    sol = eval( :(let $([:($x=$y) for (x, y) in solutions]...); $expr end) )
                    sol = Float64(sol)
                    solutions[k] = sol
                    context[k] = sol
                    done_smthg = true
                catch
                    0
                end
            end
        end
        if done_smthg==false
            finished=true
        end
    end
    if length(solutions)<length(dict)
        error("Not a triangular system")
#         throw(Exception())
    else
        return solutions
    end
end

function string_dict_to_symbol_dict(dict::Dict{String, String})
    res = Dict{Symbol, Union{Expr,Float64}}()
    for k in keys(dict)
        v = dict[k]
        res[Symbol(k)] = parse(v)
    end
    return res
end

function triangular_system(dict::Dict{String, String})
    return triangular_system(string_dict_to_symbol_dict(dict))
end
