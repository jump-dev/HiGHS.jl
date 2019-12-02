import MathOptInterface
const MOI = MathOptInterface

@testset "MOI variable count and empty" begin
    o = HiGHS.Optimizer()
    x1 = MOI.add_variable(o)
    x2 = MOI.add_constrained_variable(o, MOI.Interval(0, 1))
    @test MOI.get(o, MOI.NumberOfVariables()) == 2
    MOI.empty!(o)
    @test MOI.get(o, MOI.NumberOfVariables()) == 0
end
