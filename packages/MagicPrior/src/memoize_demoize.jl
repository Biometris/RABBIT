# Copyright (c) 2012-2013: Simon Kornblith, Stefan Karpinski, Kevin Squire,
# and other contributors:
#
# https://github.com/simonster/Memoize.jl/contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

# Chaozhi added funtion dememoize() on 12/5/2018, according to the
# https://github.com/JuliaCollections/Memoize.jl/issues/14

macro memoize(args...)
    if length(args) == 1
        dicttype = :(IdDict)
        ex = args[1]
    elseif length(args) == 2
        (dicttype, ex) = args
    else
        error("Memoize accepts at most two arguments")
    end
    # a return type declaration of Any is a No-op because everything is <: Any
    rettype = Any
    # if the return type is provided we need to strip it out and put it back later
    if ex.args[1].head == :(::)
        rettype = ex.args[1].args[2]
        ex.args[1] = ex.args[1].args[1]
    end
    # error handling for expressions that are not method definitions
    if !isa(ex,Expr) || (ex.head != :function && ex.head != Symbol("=")) ||
       isempty(ex.args) || ex.args[1].head != :call || isempty(ex.args[1].args)
        error("@memoize must be applied to a method definition")
    end
    f = ex.args[1].args[1]
    ex.args[1].args[1] = u = Symbol("##",f,"_unmemoized")

    args = ex.args[1].args[2:end]

    # Extract keywords from AST
    kws = Any[]
    vals = copy(args)
    if length(vals) > 0 && isa(vals[1], Expr) && vals[1].head == :parameters
        kws = popfirst!(vals).args
    end

    # Set up arguments for tuple to encode keywords
    tup = Array{Any}(undef,length(kws)+length(vals))
    i = 1
    for val in vals
        tup[i] = if isa(val, Expr)
                if val.head == :... || val.head == :kw
                    val.args[1]
                elseif val.head == :(::)
                    val
                else
                    error("@memoize did not understand method syntax $val")
                end
            else
                val
            end
        i += 1
    end

    for kw in kws
        if isa(kw, Expr) && (kw.head == :kw || kw.head == :...)
            tup[i] = kw.args[1]
        else
            error("@memoize did not understand method syntax")
        end
        i += 1
    end

    # Set up identity arguments to pass to unmemoized function
    identargs = Array{Any}(undef,(length(kws) > 0)+length(vals))
    i = (length(kws) > 0) + 1
    for val in vals
        if isa(val, Expr)
            if val.head == :kw
                val = val.args[1]
            end
            if isa(val, Expr) && val.head == :(::)
                val = val.args[1]
            end
        end
        identargs[i] = val
        i += 1
    end
    if length(kws) > 0
        identkws = map(kws) do kw
            if kw.head == :kw
                key = kw.args[1]
                if isa(key, Expr) && key.head == :(::)
                    key = key.args[1]
                end
                Expr(:kw, key, key)
            else
                kw
            end
        end
        identargs[1] = Expr(:parameters, identkws...)
    end

    fcachename = Symbol("##",f,"_memoized_cache")
    mod = @__MODULE__
    fcache = isdefined(mod, fcachename) ?
             getfield(mod, fcachename) :
             Core.eval(mod, :(const $fcachename = ($dicttype)()))

    if length(kws) == 0
        lookup = :($fcache[($(tup...),)]::Core.Compiler.return_type($u, typeof(($(identargs...),))))
    else
        lookup = :($fcache[($(tup...),)])
    end

    esc(quote
        $ex
        empty!($fcache)
        $f($(args...),)::$rettype =
            haskey($fcache, ($(tup...),)) ? $lookup :
            ($fcache[($(tup...),)] = $u($(identargs...),))
    end)

end

function dememoize()
  mod =  @__MODULE__
  all_caches = filter( x -> startswith(string(x), "##") &&
    endswith(string(x),"_memoized_cache"), names(mod, all=true))
  # println(all_caches)
  for cur_cache_symbol in all_caches
    empty!(getfield(mod, cur_cache_symbol))
  end
end
