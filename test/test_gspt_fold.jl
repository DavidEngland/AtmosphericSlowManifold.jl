@testset "Fold diagnostics" begin
    ms = ManifoldState()
    fold = FoldConstraint(1 - 0.35 * ms.r^2 * ms.chi - 0.12 * ms.pi_g)

    dchi = fold_transversality(fold, ms.chi)
    @test dchi !== nothing
end
