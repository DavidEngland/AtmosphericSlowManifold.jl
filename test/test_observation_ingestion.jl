using NCDatasets

@testset "Observation ingestion" begin
    mktempdir() do d
        csv_ok = joinpath(d, "tower_ok.csv")
        open(csv_ok, "w") do io
            write(io, "z_m,u_ms,v_ms,theta_k,q_kgkg,u_star_ms\n")
            write(io, "2.0,1.5,0.2,289.1,0.0041,0.35\n")
            write(io, "10.0,2.0,0.4,289.6,0.0039,0.36\n")
            write(io, "40.0,3.5,0.7,290.4,0.0031,0.40\n")
        end

        obs = read_tower_csv(csv_ok)
        @test haskey(obs.columns, :z)
        @test haskey(obs.columns, :u)
        @test obs.units[:theta] == "K"
        @test length(obs.columns[:z]) == 3

        csv_bad = joinpath(d, "tower_bad_units.csv")
        open(csv_bad, "w") do io
            write(io, "z,u,v,theta,q,u_star\n")
            write(io, "2.0,1.0,0.2,289.0,0.004,0.30\n")
        end
        @test_throws ArgumentError read_tower_csv(csv_bad)

        nc = joinpath(d, "tower_profile.nc")
        ds = NCDataset(nc, "c")
        try
            defDim(ds, "n", 3)
            zvar = defVar(ds, "z", Float64, ("n",))
            uvar = defVar(ds, "u", Float64, ("n",))
            vvar = defVar(ds, "v", Float64, ("n",))
            thvar = defVar(ds, "theta", Float64, ("n",))
            qvar = defVar(ds, "q", Float64, ("n",))
            usvar = defVar(ds, "u_star", Float64, ("n",))

            zvar[:] = [2.0, 10.0, 40.0]
            uvar[:] = [1.5, 2.0, 3.5]
            vvar[:] = [0.2, 0.4, 0.7]
            thvar[:] = [289.1, 289.6, 290.4]
            qvar[:] = [0.0041, 0.0039, 0.0031]
            usvar[:] = [0.35, 0.36, 0.40]
        finally
            close(ds)
        end

        obs_nc = read_tower_netcdf(nc)
        @test length(obs_nc.columns[:z]) == 3
        @test obs_nc.units[:q] == "kg kg^-1"
    end
end
