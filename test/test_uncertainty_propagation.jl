using Test
using AtmosphericSlowManifold

@testset "Uncertainty propagation" begin
    terms = OperatorTerm{Float64}[
        OperatorTerm{Float64}(0.0, BasisOperator[BasisOperator(StateVariable(:z), 1.0)]),
        OperatorTerm{Float64}(0.0, BasisOperator[]),
    ]
    model = DiscoveredModel{Float64}(:K_m, terms, 0.0, 2)

    z = [0.0, 1.0, 2.0]
    draws_mat = [
        1.0 2.0 3.0;
        0.0 1.0 2.0
    ]

    out = evaluate_profile_uncertainty(model, draws_mat, z)
    @test size(out.profiles) == (3, 3)
    @test out.profiles[:, 1] ≈ [0.0, 1.0, 2.0]
    @test out.profiles[:, 2] ≈ [1.0, 3.0, 5.0]
    @test out.profiles[:, 3] ≈ [2.0, 5.0, 8.0]
    @test out.median ≈ [1.0, 3.0, 5.0]
    @test all(out.lower .<= out.median)
    @test all(out.median .<= out.upper)

    draws_dict = Dict{Symbol, Any}(
        :c_1 => [1.0, 2.0, 3.0],
        :c_2 => [0.0, 1.0, 2.0],
    )
    out_dict = evaluate_profile_uncertainty(model, draws_dict, z)
    @test out_dict.profiles ≈ out.profiles
    @test out_dict.median ≈ out.median

    bad_model = DiscoveredModel{Float64}(
        :K_m,
        OperatorTerm{Float64}[OperatorTerm{Float64}(0.0, BasisOperator[BasisOperator(DiagnosticVariable(:Ri), 1.0)])],
        0.0,
        1,
    )
    @test_throws ArgumentError evaluate_profile_uncertainty(bad_model, [1.0 2.0 3.0], z)

    out_with_feature = evaluate_profile_uncertainty(
        bad_model,
        [1.0 2.0 3.0],
        z;
        feature_profiles = Dict{Symbol, Any}(:Ri => [0.1, 0.2, 0.3]),
    )
    @test size(out_with_feature.profiles) == (3, 3)

    zero_model = DiscoveredModel{Float64}(:K_m, OperatorTerm{Float64}[], 0.0, 0)
    zero_out = evaluate_profile_uncertainty(zero_model, zeros(0, 5), z)
    @test all(zero_out.median .== 0.0)
    @test all(zero_out.lower .== 0.0)
    @test all(zero_out.upper .== 0.0)
end
