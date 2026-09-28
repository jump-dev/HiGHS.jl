# Copyright (c) 2019 Mathieu Besançon, Oscar Dowson, and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module MyApp

import HiGHS
import MathOptInterface as MOI

function @main(args::Vector{String})::Cint
    capacity = 10.0
    profit = [5.0, 3.0, 2.0, 7.0, 4.0]
    weight = [2.0, 8.0, 4.0, 2.0, 5.0]
    model = HiGHS.Optimizer()
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
