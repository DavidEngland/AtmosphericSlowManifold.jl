using Symbolics

@testset "Discovery parametric IR" begin
    lib = build_feature_library([:u], [:Ri], 1)
    @test length(lib.features) == 3
    @test length(lib.evaluators) == 3

    sample = Dict{Symbol, Float64}(:u => 2.0, :d1_u => -0.5, :Ri => 3.0)
    vals = [f(sample) for f in lib.evaluators]
    @test vals[1] == 2.0
    @test vals[2] == -0.5
    @test vals[3] == 3.0

    term1 = OperatorTerm(1.5, [BasisOperator(StateVariable(:u), 2.0)])
    term2 = OperatorTerm(0.5, [BasisOperator(DiagnosticVariable(:Ri), 1.0), BasisOperator(SpatialDerivative(:u, 1), 1.0)])

    @variables u Ri d1_u
    var_map = Dict(:u => u, :Ri => Ri, :d1_u => d1_u)
    model = DiscoveredModel(:K_m, [term1, term2], 0.1, 2)
    expr = to_mtk_expression(model, var_map)

    exv = Symbolics.substitute(expr, Dict(u => 2.0, Ri => 3.0, d1_u => -0.5))
    @test isapprox(float(Symbolics.value(exv)), 5.25; atol = 1e-10)

    constraints = AbstractPhysicalConstraint[
        PositivityConstraint(StateVariable(:u)),
        MonotonicityConstraint(DiagnosticVariable(:Ri), SpatialDerivative(:u, 1)),
        EnergyConstraint(DiagnosticVariable(:Ri)),
    ]
    grid = [
        2.0 -0.5 3.0
        1.5 -0.2 2.5
        1.0 -0.1 2.0
    ]
    cm = assemble_constraint_matrix(constraints, lib.features, grid)
    @test size(cm.A_ineq, 2) == length(lib.features)
    @test size(cm.A_ineq, 1) == length(constraints)
    @test length(cm.b_ineq) == length(constraints)
end
