using Test
using AtmosphericSlowManifold
using Symbolics

@testset "LaTeX exporter" begin
    terms = OperatorTerm{Float64}[
        OperatorTerm{Float64}(0.25, BasisOperator[BasisOperator(SpatialDerivative(:u, 2), 1.0)]),
        OperatorTerm{Float64}(1.5, BasisOperator[BasisOperator(StateVariable(:z), 1.0), BasisOperator(SpatialDerivative(:u, 1), 1.0)]),
        OperatorTerm{Float64}(-0.3, BasisOperator[BasisOperator(DiagnosticVariable(:Ri), 1.0)]),
    ]
    model = DiscoveredModel{Float64}(:u, terms, 0.0123, 3)

    eqn = to_latex(model)
    @test occursin("\\frac{\\partial u}{\\partial t}", eqn)
    @test occursin("\\frac{\\partial^{2} u}{\\partial z^{2}}", eqn)
    @test occursin("\\mathrm{Ri}", eqn)

    ttable = latex_term_table(model; r2 = 0.97)
    @test occursin("Term & Basis & Coefficient", ttable)
    @test occursin("R^{2} = 0.97", ttable)
    @test occursin("Residual norm = 0.0123", ttable)

    models = Dict(
        :CASES99 => model,
        :GABLS1 => DiscoveredModel{Float64}(:u, terms[1:2], 0.02, 2),
    )
    stable = latex_site_summary_table(models; r2_by_site = Dict(:CASES99 => 0.94, :GABLS1 => 0.91))
    @test occursin("Site & Active terms", stable)
    @test occursin("0.94", stable)
    @test occursin("0.91", stable)

    multieqn = to_latex([model, DiscoveredModel{Float64}(:theta, terms[1:1], 0.05, 1)])
    @test occursin("\\begin{align}", multieqn)
    @test occursin("\\partial \\theta", multieqn)

    temp_path = joinpath(mktempdir(), "latex-export.tex")
    written_path = write_latex(temp_path, eqn)
    @test written_path == temp_path
    @test isfile(temp_path)
    @test occursin("\\frac{\\partial u}{\\partial t}", read(temp_path, String))

    @variables z
    closure = WSINDyClosure(Num(0.1) + Num(z), Num(0.2), Num(0.05))
    mom = to_latex(closure; field = :momentum)
    heat = to_latex(closure; field = :heat)
    flux = to_latex(closure; field = :flux)
    @test occursin("0.1", mom)
    @test occursin("0.2", heat)
    @test occursin("0.05", flux)
    @test_throws ArgumentError to_latex(closure; field = :badfield)
end
