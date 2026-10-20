# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module MyApp

import HiGHS
import MathOptInterface as MOI

# We should fix this upstream in MathOptInterface.jl
function Base.show(io::IO, x::MOI.VariableIndex)
    return print(io, "MOI.VariableIndex(", x.value, ")")
end

function @main(args::Vector{String})::Cint
    capacity = 10.0
    profit = [5.0, 3.0, 2.0, 7.0, 4.0]
    weight = [2.0, 8.0, 4.0, 2.0, 5.0]
    model = HiGHS.Optimizer()
    silent = "--silent" in args
    MOI.set(model, MOI.Silent(), silent)
    @assert MOI.get(model, MOI.Silent()) == silent
    param = MOI.RawOptimizerAttribute("presolve")
    try
        MOI.set(model, param, length(args))
        error("Expected SetAttributeNotAllowed")
    catch err
        @assert err isa MOI.SetAttributeNotAllowed{MOI.RawOptimizerAttribute}
        @assert occursin("Expected a value of type `String`.", err.message)
    end
    x = MOI.add_variables(model, 5)
    MOI.add_constraint.(model, x, MOI.ZeroOne())
    MOI.add_constraint(model, weight' * x, MOI.LessThan(capacity))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
    f = profit' * x
    MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
    MOI.optimize!(model)
    @assert MOI.get(model, MOI.TerminationStatus()) == MOI.OPTIMAL
    return 0
end

end  # MyApp
