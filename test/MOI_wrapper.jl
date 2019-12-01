import MathOptInterface
const MOI = MathOptInterface

@testset "MOI wrapper" begin
    o = HiGHS.Optimizer()
    x1 = MOI.add_variable(o)
    x2 = MOI.add_constrained_variable(o, MOI.Interval(0, 1))
end
