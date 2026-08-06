@testset "Gegenbauer transform" begin
    obs = ObservationTable(
        Dict(
            :z => [2.0, 10.0, 20.0, 40.0, 80.0],
            :u => [1.0, 1.6, 2.1, 2.8, 3.2],
            :v => [0.1, 0.2, 0.25, 0.35, 0.4],
            :theta => [288.8, 289.0, 289.3, 289.9, 290.6],
            :q => [0.0042, 0.0040, 0.0037, 0.0033, 0.0030],
            :u_star => [0.31, 0.31, 0.31, 0.31, 0.31],
        ),
        Dict(
            :z => "m",
            :u => "m s^-1",
            :v => "m s^-1",
            :theta => "K",
            :q => "kg kg^-1",
            :u_star => "m s^-1",
        ),
    )

    proj = project_to_gegenbauer(obs; n_modes = 6, lambda = 0.75)
    @test proj.n_modes == 6
    @test :u in proj.variables
    @test length(proj.coefficients[:u]) == 6
    @test all(isfinite, proj.coefficients[:theta])

    obs_small = ObservationTable(
        Dict(
            :z => [2.0, 10.0],
            :u => [1.0, 1.2],
        ),
        Dict(
            :z => "m",
            :u => "m s^-1",
        ),
    )
    @test_throws ArgumentError project_to_gegenbauer(obs_small)
end
