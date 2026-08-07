@testset "Canard detection" begin
    g_fn(y, z) = [y[1]^2 + z[1]]

    # Saddle case: nonzero eigenvalues have opposite signs.
    f_saddle(y, z) = [y[1] + z[2], z[1] + 2z[2]]
    saddle = classify_folded_singularity([0.0], [0.0, 0.0], g_fn, f_saddle)
    @test saddle.classification == :folded_saddle
    @test abs(saddle.det_fast_jacobian) <= 1e-8

    # Node case: nonzero eigenvalues are real and same sign.
    f_node(y, z) = [y[1] + 2z[1], y[1] + z[2]]
    node = classify_folded_singularity([0.0], [0.0, 0.0], g_fn, f_node)
    @test node.classification == :folded_node

    # Focus case: nonzero eigenvalues are complex.
    f_focus(y, z) = [y[1] + z[2], y[1] - z[2]]
    focus = classify_folded_singularity([0.0], [0.0, 0.0], g_fn, f_focus)
    @test focus.classification == :folded_focus

    # Matrix overload should also classify a focus.
    @test classify_folded_singularity([0.0 -1.0; 1.0 0.0]) == :folded_focus

    t = collect(0.0:0.1:1.0)
    states = [
        [0.4, 0.0],
        [0.2, 0.0],
        [-0.5, 0.0],
        [-0.7, 0.0],
        [-0.8, 0.0],
        [-0.6, 0.0],
        [-0.4, 0.0],
        [0.1, 0.0],
        [0.3, 0.0],
        [0.4, 0.0],
        [0.5, 0.0],
    ]

    segments = detect_canard_trajectories(
        t,
        states,
        g_fn;
        ny = 1,
        fold_threshold = 1e-10,
        persistence_threshold = 0.3,
    )
    @test length(segments) == 1
    @test segments[1].duration >= 0.3

    # Solution-like object overload.
    sol_like = (t = t, u = states)
    segments2 = detect_canard_trajectories(
        sol_like,
        g_fn;
        ny = 1,
        fold_threshold = 1e-10,
        persistence_threshold = 0.3,
    )
    @test length(segments2) == 1
end
