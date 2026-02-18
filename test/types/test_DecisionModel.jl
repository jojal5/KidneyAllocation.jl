@testset "DecisionModel.jl" begin

    import KidneyAllocation: fit_threshold_f1, fit_threshold_prevalence

    @testset "fit_threshold_f1" begin

        @test_throws AssertionError fit_threshold_f1([true, false, false], [0.8, 0.7])
        @test_throws ArgumentError fit_threshold_f1([true, false, false], [0.8, 0.7, 1.1])

        # No ties
        gt = [true, true, false, true, false]
        p = [0.7, 0.6, 0.5, 0.4, 0.3]
        @test fit_threshold_f1(gt, p) ≈ 0.4

        # with ties
        gt = [true, true, false, true, false]
        p = [0.7, 0.6, 0.5, 0.5, 0.4]
        @test fit_threshold_f1(gt, p) ≈ 0.5

    end

    @testset "fit_threshold_prevalence" begin

        @test_throws AssertionError fit_threshold_prevalence([true, false, false], [0.8, 0.7])
        @test_throws ArgumentError fit_threshold_prevalence([true, false, false], [0.8, 0.7, 1.1])

        # No ties
        gt = [true, true, false, true, false]
        p = [0.7, 0.6, 0.5, 0.4, 0.3]
        @test fit_threshold_prevalence(gt, p) ≈ 0.5

        # with ties
        gt = [true, true, false, true, false]
        p = [0.7, 0.7, 0.5, 0.5, 0.4]
        @test fit_threshold_prevalence(gt, p) ≈ 0.5

    end



end