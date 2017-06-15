# This duplicates some functionality from Dolo.jl. Should be moved to Dolang.

sanitize(s, dynvars::Array{Symbol,1}) = s
sanitize(s::Symbol, dynvars::Array{Symbol,1})= s in dynvars ? :($s(0)) : s
function sanitize(eq::Expr, dynvars::Array{Symbol,1})
    if eq.head == :(=)
        return sanitize( :( $(eq.args[1])==$( eq.args[2] ) ), dynvars )
    else
        if eq.head == :call
            if eq.args[1] in dynvars
                return eq
            else
                return Expr(eq.head, [sanitize(e, dynvars) for e in eq.args]...)
            end
        else
            head = eq.head
            args = eq.args
            return Expr(eq.head, [sanitize(e, dynvars) for e in args]...)
        end
    end
end
