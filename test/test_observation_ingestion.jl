# test/test_observation_ingestion.jl
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

@testset "Observation path utilities and mapped ingestion" begin
    data_dir = resolve_sibling_data_dir(; must_exist = false)
    @test endswith(data_dir, joinpath("SpectralBL-Analytics", "data"))

    mktempdir() do d
        mkpath(joinpath(d, "nested"))
        csv_path = joinpath(d, "profiles.csv")
        nc_path = joinpath(d, "nested", "profiles.nc")
        txt_path = joinpath(d, "notes.txt")

        open(csv_path, "w") do io
            write(io, "height,u_velocity,v_velocity,potential_temp\n")
            write(io, "2.0,1.0,0.1,289.0\n")
            write(io, "20.0,2.0,0.2,290.0\n")
        end

        open(txt_path, "w") do io
            write(io, "ignore me\n")
        end

        ds = NCDataset(nc_path, "c")
        try
            defDim(ds, "n", 1)
            zvar = defVar(ds, "z", Float64, ("n",))
            uvar = defVar(ds, "u", Float64, ("n",))
            vvar = defVar(ds, "v", Float64, ("n",))
            tvar = defVar(ds, "theta", Float64, ("n",))
            zvar[:] = [5.0]
            uvar[:] = [1.2]
            vvar[:] = [0.0]
            tvar[:] = [289.5]
        finally
            close(ds)
        end

        files = find_data_files(d; extensions = [".csv", ".nc"], recursive = true)
        @test length(files) == 2
        @test csv_path in files
        @test nc_path in files

        obs = read_observation_data(
            csv_path;
            z_col = :height,
            u_col = :u_velocity,
            v_col = :v_velocity,
            temp_col = :potential_temp,
        )
        @test obs isa ObservationTable
        @test obs.columns[:z] == [2.0, 20.0]
        @test obs.columns[:u] == [1.0, 2.0]
        @test obs.columns[:q] == [0.0, 0.0]
        @test obs.columns[:u_star] == [0.3, 0.3]
    end
end

@testset "Observation mapped ingestion with surface flux aliases" begin
    mktempdir() do d
        csv_flux = joinpath(d, "profiles_flux.csv")
        open(csv_flux, "w") do io
            write(io, "height,u_velocity,v_velocity,potential_temp,ustar,hs\n")
            write(io, "2.0,1.0,0.1,289.0,0.30,-20.0\n")
            write(io, "20.0,2.0,0.2,290.0,0.35,-25.0\n")
        end

        obs = read_observation_data(
            csv_flux;
            z_col = :height,
            u_col = :u_velocity,
            v_col = :v_velocity,
            temp_col = :potential_temp,
            auto_surface_flux_aliases = true,
            include_derived_obukhov = true,
        )

        @test haskey(obs.columns, :sensible_heat_flux)
        @test haskey(obs.columns, :L_obukhov)
        @test obs.units[:sensible_heat_flux] == "W m^-2"
        @test obs.units[:L_obukhov] == "m"
        @test all(isfinite, obs.columns[:L_obukhov])
        @test length(obs.columns[:L_obukhov]) == 2
    end
end
