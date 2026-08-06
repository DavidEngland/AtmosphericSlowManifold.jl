using DataFrames
using JSON3
using NCDatasets

@testset "Export utilities" begin
    mktempdir() do d
        # DataFrame -> CSV
        df = DataFrame(z = [0.0, 10.0], u = [1.0, 2.0])
        csv_path = joinpath(d, "table.csv")
        out_csv = export_to_csv(csv_path, df)
        @test out_csv == csv_path
        @test isfile(csv_path)

        # ObservationTable -> CSV
        obs = ObservationTable(
            Dict(
                :z => [2.0, 20.0],
                :u => [1.0, 2.0],
                :v => [0.1, 0.2],
                :theta => [289.0, 290.0],
                :q => [0.0, 0.0],
                :u_star => [0.3, 0.3],
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
        obs_csv_path = joinpath(d, "obs.csv")
        out_obs_csv = export_to_csv(obs_csv_path, obs)
        @test out_obs_csv == obs_csv_path
        @test isfile(obs_csv_path)

        # Modal budget -> CSV
        budget = ModalBudgetDiagnostic(
            [1.0, 2.0],
            [0.1, 0.2],
            [0.05, 0.15],
            [0.85, 1.65],
        )
        budget_csv_path = joinpath(d, "budget.csv")
        out_budget_csv = export_to_csv(budget_csv_path, budget, [0.0, 10.0])
        @test out_budget_csv == budget_csv_path
        @test isfile(budget_csv_path)

        # Discovered model -> JSON
        terms = OperatorTerm{Float64}[
            OperatorTerm{Float64}(0.2, BasisOperator[BasisOperator(StateVariable(:u), 1.0)]),
            OperatorTerm{Float64}(-0.01, BasisOperator[BasisOperator(DiagnosticVariable(:Ri), 1.0)]),
        ]
        model = DiscoveredModel{Float64}(:K_m, terms, 0.01, 2)
        diagnostics = Dict{Symbol, Any}(:status => :ok, :samples => 100)
        json_path = joinpath(d, "model.json")
        out_json = export_to_json(json_path, model, diagnostics)
        @test out_json == json_path
        @test isfile(json_path)

        parsed = JSON3.read(read(json_path, String))
        @test String(parsed["target_variable"]) == "K_m"
        @test Int(parsed["num_terms"]) == 2

        # NetCDF trajectory export
        z_grid = [0.0, 10.0, 20.0]
        t_grid = [0.0, 5.0]
        u_matrix = [1.0 1.1; 2.0 2.1; 3.0 3.2]
        nc_path = joinpath(d, "traj.nc")
        out_nc = export_to_netcdf(nc_path, z_grid, t_grid, u_matrix, "u")
        @test out_nc == nc_path
        @test isfile(nc_path)

        ds = NCDataset(nc_path)
        try
            @test haskey(ds, "z")
            @test haskey(ds, "t")
            @test haskey(ds, "u")
            @test vec(ds["z"][:]) == z_grid
            @test vec(ds["t"][:]) == t_grid
            @test size(ds["u"][:, :]) == (3, 2)
        finally
            close(ds)
        end
    end
end
