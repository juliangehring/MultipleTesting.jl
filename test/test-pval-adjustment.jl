## Test pvalue adjustment methods ##
module Test_pval_adjustment

using MultipleTesting
using Test


@testset "p-Value adjustment" begin

    pi0 = 0.4

    # p-values without ties
    pval1 = [0.0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.75, 1.0]
    ref1 = Dict(
        Bonferroni         => [0.0, 0.001, 0.01, 0.1, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0],
        Holm               => [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0],
        Hochberg           => [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0],
        Hommel             => [0.0, 9e-4, 8e-3, 0.07, 0.3, 0.5, 0.8, 1.0, 1.0, 1.0],
        BenjaminiHochberg  => [0.0, 0.0005, 0.003333333, 0.025, 0.1, 0.166666667, 0.285714286, 0.5, 0.833333333, 1.0],
        BenjaminiYekutieli => [0.0, 0.001464484, 0.009763228, 0.073224206, 0.292896825, 0.488161376, 0.836848073, 1.0, 1.0, 1.0],
        BenjaminiLiu       => [0.0, 0.0008096761, 0.0063776447, 0.0475542565, 0.1589448656, 0.2047550000, 0.2361600000, 0.2361600000, 0.2361600000, 0.2361600000],
        Sidak              => [0.0, 0.0009995501, 0.0099551198, 0.0956179250, 0.4012630608, 0.6513215599, 0.8926258176, 0.9939533824, 0.9999990463, 1.0000000000],
        ForwardStop        => [0.0, 0.0000500025, 0.0003668351, 0.0027877103, 0.0124888271, 0.0279674419, 0.0558497432, 0.1127217283, 0.2542297986, 1.0],
        BarberCandes       => [0.2857142857, 0.2857142857, 0.2857142857, 0.2857142857, 0.2857142857, 0.2857142857, 0.2857142857, 0.375, 1.0, 1.0]
    )

    # p-values with ties
    pval2 = [0.0001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.1, 0.4, 0.75, 1.0]
    ref2 = Dict(
        Bonferroni         => [0.001, 0.001, 0.01, 0.1, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0],
        Holm               => [0.001, 0.001, 0.008, 0.07, 0.3, 0.5, 0.5, 1.0, 1.0, 1.0],
        Hochberg           => [9e-4, 9e-4, 8e-3, 0.07, 0.3, 0.4, 0.4, 1.0, 1.0, 1.0],
        Hommel             => [9e-4, 9e-4, 8e-3, 0.07, 0.25, 0.4, 0.4, 1.0, 1.0, 1.0],
        BenjaminiHochberg  => [0.0005, 0.0005, 0.003333333, 0.025, 0.1, 0.142857143, 0.142857143, 0.5, 0.833333333, 1.0],
        BenjaminiYekutieli => [0.001464484, 0.001464484, 0.009763228, 0.073224206, 0.292896825, 0.418424036, 0.418424036, 1.0, 1.0, 1.0],
        BenjaminiLiu       => [0.0009995501, 0.0009995501, 0.0063776447, 0.0475542565, 0.1589448656, 0.2047550000, 0.2047550000, 0.2352000000, 0.2352000000, 0.2352000000],
        Sidak              => [0.0009995501, 0.0009995501, 0.0099551198, 0.0956179250, 0.4012630608, 0.6513215599, 0.6513215599, 0.9939533824, 0.9999990463, 1.0000000000],
        ForwardStop        => [0.0001000050, 0.0001000050, 0.0004001701, 0.0028127115, 0.0125088281, 0.0279841094, 0.0390378817, 0.0980113495, 0.2411539063, 1.0],
        BarberCandes       => [0.2857142857, 0.2857142857, 0.2857142857, 0.2857142857, 0.2857142857, 0.2857142857, 0.2857142857, 0.375, 1.0, 1.0]
    )

    # smallest p-values from larger set
    k3 = 6
    n3 = length(pval1)
    pval3 = pval1[1:k3]
    pval3pad = [pval3; fill(1.0, n3-k3)]  # non-observed p-values equal to 1
    methods3 = [Bonferroni, Holm, Hochberg, Hommel, BenjaminiHochberg,
                BenjaminiYekutieli, BenjaminiLiu, Sidak, ForwardStop]

    @testset "pvalue adjustment $method" for method in keys(ref1)

        @test method <: PValueAdjustment
        @test typeof(method()) <: PValueAdjustment

        @test_throws MethodError method(0.1)

        ## no valid p-values as input
        @test_throws DomainError adjust([-1.0, 0.5], method())
        @test_throws DomainError adjust([0.5, 1.5], method())

        ## any single p-value is returned unchanged
        # this test is not valid for ForwardStop or Barber-Candès
        if !(method in [ForwardStop, BarberCandes])
            pval = rand(1)
            @test adjust(pval, method()) == pval
        end

        if (method == BarberCandes)
            pval = rand(1)
            @test adjust(pval, method()) == fill(1.0, size(pval))
        end

        ## compare with reference values
        @test isapprox( adjust(pval1, method()), ref1[method], atol = 1e-9 )
        @test isapprox( adjust(PValues(pval1), method()), ref1[method], atol = 1e-9 )

        # unsorted inputs
        for i in 1:10
            ord = MultipleTesting.unorder(pval1)
            @test isapprox( adjust(pval1[ord], method()), ref1[method][ord], atol = 1e-9 )
        end


        ## compare with reference values having ties
        @test isapprox( adjust(pval2, method()), ref2[method], atol = 1e-9 )
        @test isapprox( adjust(PValues(pval2), method()), ref2[method], atol = 1e-9 )

        # unsorted inputs
        # this test is not valid for ForwardStop
        if method != ForwardStop
            for i in 1:10
                ord = MultipleTesting.unorder(ref2[method])
                @test isapprox( adjust(pval2[ord], method()), ref2[method][ord], atol = 1e-9 ) # FIXME
            end
        end

        ## sorting order does not play a role
        for i in 1:10
            pval4 = sort(rand(10)) # all under H0
            ord = MultipleTesting.unorder(pval4)
            @test adjust(pval4[ord], method()) == adjust(pval4, method())[ord]
        end


        ## total number of tests explicitly specified
        if method in methods3
            # all p-values present: same reference values as for pval1
            @test adjust(pval3, length(pval3), method()) == adjust(pval3, method())
            @test adjust(PValues(pval3), length(pval3), method()) == adjust(pval3, method())

            # k smallest p-values present, total number n known
            @test adjust(pval3, n3, method()) == adjust(pval3pad, n3, method())[1:k3]

            # unsorted inputs
            for i in 1:10
                ord = MultipleTesting.unorder(pval3)
                @test adjust(pval3[ord], n3, method()) == (adjust(pval3pad, n3, method())[1:k3])[ord]
            end

            # k > n not allowed: test for any n in [1,k-1]
            @test_throws ArgumentError adjust(pval3, rand(1:k3-1), method())
        end

    end


    @testset "argument checking for number of tests" begin

        @test_throws ArgumentError MultipleTesting.check_number_tests(3, 2)
        @test MultipleTesting.check_number_tests(2, 2) == nothing
        @test MultipleTesting.check_number_tests(2, 3) == nothing

    end


    @testset "BarberCandès #2:" begin
        for k = 1:5
            pv = rand(BetaUniformMixtureModel(0.5, 0.5, 7.0), 40)
            @test isapprox(adjust(pv, BarberCandes()),
                           MultipleTesting.barber_candes_brute_force(pv), atol=1e-9)
        end
    end

    @testset "BarberCandès: All p-values < 0.5 (#87)" begin
        for pv in ([0.05, 0.1, 0.3], [0.01, 0.17, 0.25, 0.37, 0.47])
            n = length(pv)
            @test isapprox( adjust(PValues(pv), BarberCandes()), fill(1/n, n) )
        end
    end

    @testset "Step-up/down" begin

        x        = [1.0, 2.0, 2.0, 5.0, 4.0, 5.0, 7.0]
        ref_up   = [1.0, 2.0, 2.0, 4.0, 4.0, 5.0, 7.0]
        ref_down = [1.0, 2.0, 2.0, 5.0, 5.0, 5.0, 7.0]

        # step-up
        y = copy(x)
        MultipleTesting.stepup!(y)
        @test y == ref_up

        # step-down
        y = copy(x)
        MultipleTesting.stepdown!(y)
        @test y == ref_down

    end

end

end
