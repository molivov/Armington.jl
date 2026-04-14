using Test
using Armington

@testset "Armington.jl" begin

    # ================================================================
    # Construction
    # ================================================================
    @testset "Construction" begin
        @testset "Basic construction" begin
            leaf = CESLeaf(:A)
            @test leaf.name == :A

            leaf_str = CESLeaf("B")
            @test leaf_str.name == :B

            node = CESNode(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test node.σ == 2.0
            @test node.α == (0.5, 0.5)
        end

        @testset "Type promotion" begin
            # Int σ with Float64 α → Float64
            node = CESNode(2, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test node.σ isa Float64

            # All Int → Int (if someone really wants that)
            node_int = CESNode(2, (1, 1), (CESLeaf(:A), CESLeaf(:B)))
            @test node_int.σ isa Int

            # BigFloat
            node_big = CESNode(big"2.0", (big"0.5", big"0.5"), (CESLeaf(:A), CESLeaf(:B)))
            @test node_big.σ isa BigFloat
            @test eltype(node_big.α) == BigFloat
        end

        @testset "Vector/iterable convenience" begin
            node = CESNode(2.0, [0.5, 0.5], [CESLeaf(:A), CESLeaf(:B)])
            @test node.α == (0.5, 0.5)
        end

        @testset "Anonymous leaves" begin
            leaf = CESLeaf()
            @test leaf.name == Symbol()

            node = CESNode(2.0, (0.5, 0.5), (CESLeaf(), CESLeaf()))
            @test aggregate(node, [4.0, 9.0]) ≈ 6.25
        end

        @testset "Leaf-free construction" begin
            # Two-arg form: leaves created automatically
            node = CESNode(2.0, (0.5, 0.5))
            @test aggregate(node, [4.0, 9.0]) ≈ 6.25

            # Works with vectors too
            node_v = CESNode(2.0, [0.5, 0.5])
            @test aggregate(node_v, [4.0, 9.0]) ≈ 6.25

            # Three inputs
            node3 = CESNode(3.0, (0.4, 0.35, 0.25))
            @test aggregate(node3, [2.0, 3.0, 4.0]) > 0
        end

        @testset "Leaf-free nesting" begin
            # Inner node with anonymous leaves, outer with named + nested
            inner = CESNode(4.0, (0.6, 0.4))
            tree = CESNode(1.5, (0.7, 0.3), (CESLeaf(:dom), inner))
            @test aggregate(tree, [1.0, 1.0, 1.0]) > 0
        end

        @testset "Invalid construction" begin
            @test_throws ArgumentError CESNode(-1.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test_throws ArgumentError CESNode(2.0, (0.5, -0.1), (CESLeaf(:A), CESLeaf(:B)))
            @test_throws ArgumentError CESNode(2.0, (0.5,), (CESLeaf(:A), CESLeaf(:B)))
        end

        @testset "Limiting σ values accepted" begin
            leon = CESNode(0.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test leon.σ == 0.0
            lin = CESNode(Inf, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test isinf(lin.σ)
        end
    end

    # ================================================================
    # Tree introspection
    # ================================================================
    @testset "Tree introspection" begin
        inner = CESNode(3.0, (0.6, 0.4), (CESLeaf(:US), CESLeaf(:EU)))
        tree = CESNode(1.5, (0.7, 0.3), (CESLeaf(:dom), inner))

        @test leaf_names(tree) == [:dom, :US, :EU]
        @test leaf_names(CESLeaf(:solo)) == [:solo]

        sub = CESNode(5.0, (0.5, 0.5), (CESLeaf(:v1), CESLeaf(:v2)))
        mid = CESNode(3.0, (0.4, 0.6), (sub, CESLeaf(:C)))
        top = CESNode(1.5, (0.5, 0.5), (CESLeaf(:A), mid))
        @test leaf_names(top) == [:A, :v1, :v2, :C]
    end

    # ================================================================
    # General CES
    # ================================================================
    @testset "General CES" begin
        @testset "Equal weights, equal inputs" begin
            tree = CESNode(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test aggregate(tree, [1.0, 1.0]) ≈ 1.0
        end

        @testset "Known analytic value (σ = 2)" begin
            tree = CESNode(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            # σ=2 → ρ=1/2, Q = [0.5*4^(1/2) + 0.5*9^(1/2)]^2 = 6.25
            @test aggregate(tree, [4.0, 9.0]) ≈ 6.25
        end

        @testset "Three inputs" begin
            tree = CESNode(3.0, (0.4, 0.35, 0.25),
                (CESLeaf(:A), CESLeaf(:B), CESLeaf(:C)))
            Q = aggregate(tree, [2.0, 3.0, 4.0])
            @test Q > 0
        end
    end

    # ================================================================
    # Limiting cases
    # ================================================================
    @testset "Leontief (σ = 0)" begin
        tree = CESNode(0.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
        @test aggregate(tree, [1.0, 100.0]) ≈ 2.0
        @test aggregate(tree, [10.0, 3.0]) ≈ 6.0

        # Asymmetric weights
        tree2 = CESNode(0.0, (0.25, 0.75), (CESLeaf(:A), CESLeaf(:B)))
        # min(2/0.25, 6/0.75) = min(8, 8) = 8
        @test aggregate(tree2, [2.0, 6.0]) ≈ 8.0
    end

    @testset "Cobb-Douglas (σ = 1)" begin
        tree = CESNode(1.0, (0.3, 0.7), (CESLeaf(:A), CESLeaf(:B)))
        x = [2.0, 3.0]
        expected = (2.0 / 0.3)^0.3 * (3.0 / 0.7)^0.7
        @test aggregate(tree, x) ≈ expected atol = 1e-12

        # Shares should work as exponents
        tree_eq = CESNode(1.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
        @test aggregate(tree_eq, [4.0, 4.0]) ≈ (4.0 / 0.5)^0.5 * (4.0 / 0.5)^0.5
    end

    @testset "Linear / perfect substitutes (σ = Inf)" begin
        tree = CESNode(Inf, (0.6, 0.4), (CESLeaf(:A), CESLeaf(:B)))
        @test aggregate(tree, [3.0, 7.0]) ≈ 4.6

        # Only one input matters
        @test aggregate(tree, [10.0, 0.0]) ≈ 6.0
    end

    # ================================================================
    # Nesting
    # ================================================================
    @testset "Two-level nesting" begin
        inner = CESNode(3.0, (0.6, 0.4), (CESLeaf(:US), CESLeaf(:EU)))
        tree = CESNode(1.5, (0.7, 0.3), (CESLeaf(:dom), inner))
        Q = aggregate(tree, [1.0, 1.0, 1.0])
        @test Q > 0
    end

    @testset "Three-level nesting" begin
        sub = CESNode(5.0, (0.5, 0.5), (CESLeaf(:v1), CESLeaf(:v2)))
        mid = CESNode(3.0, (0.4, 0.6), (sub, CESLeaf(:C)))
        top = CESNode(1.5, (0.5, 0.5), (CESLeaf(:A), mid))
        Q = aggregate(top, [1.0, 1.1, 0.95, 1.3])
        @test Q > 0
    end

    @testset "Mixed limiting cases in nesting" begin
        # Leontief at top, high σ at bottom
        inner = CESNode(10.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
        tree = CESNode(0.0, (0.5, 0.5), (CESLeaf(:C), inner))
        Q = aggregate(tree, [4.0, 3.0, 5.0])
        @test Q > 0

        # Linear at top, Cobb-Douglas at bottom
        inner_cd = CESNode(1.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
        tree_lin = CESNode(Inf, (0.6, 0.4), (inner_cd, CESLeaf(:C)))
        Q2 = aggregate(tree_lin, [2.0, 3.0, 5.0])
        @test Q2 > 0
    end

    # ================================================================
    # Homogeneity: Q(λx) = λQ(x)
    # ================================================================
    @testset "Homogeneity" begin
        λ = 3.7

        @testset "General CES" begin
            tree = CESNode(2.0, (0.6, 0.4), (CESLeaf(:A), CESLeaf(:B)))
            x = [1.5, 2.3]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x) atol = 1e-10
        end

        @testset "Cobb-Douglas" begin
            tree = CESNode(1.0, (0.4, 0.6), (CESLeaf(:A), CESLeaf(:B)))
            x = [2.0, 5.0]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x) atol = 1e-10
        end

        @testset "Leontief" begin
            tree = CESNode(0.0, (0.3, 0.7), (CESLeaf(:A), CESLeaf(:B)))
            x = [2.0, 5.0]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x)
        end

        @testset "Linear" begin
            tree = CESNode(Inf, (0.6, 0.4), (CESLeaf(:A), CESLeaf(:B)))
            x = [3.0, 7.0]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x)
        end

        @testset "Nested" begin
            inner = CESNode(4.0, (0.5, 0.5), (CESLeaf(:x), CESLeaf(:y)))
            tree = CESNode(1.5, (0.7, 0.3), (CESLeaf(:A), inner))
            x = [1.0, 2.0, 1.5]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x) atol = 1e-10
        end
    end

    # ================================================================
    # Monotonicity: increasing any input increases Q
    # ================================================================
    @testset "Monotonicity" begin
        tree = CESNode(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
        x = [3.0, 5.0]
        Q_base = aggregate(tree, x)
        @test aggregate(tree, [3.1, 5.0]) > Q_base
        @test aggregate(tree, [3.0, 5.1]) > Q_base
    end

    # ================================================================
    # CESNode_ρ
    # ================================================================
    @testset "CESNode_ρ" begin
        @testset "Equivalence with CESNode" begin
            tree_rho = CESNode_ρ(0.5, (0.6, 0.4), (CESLeaf(:A), CESLeaf(:B)))
            tree_sig = CESNode(2.0, (0.6, 0.4), (CESLeaf(:A), CESLeaf(:B)))
            x = [3.0, 7.0]
            @test aggregate(tree_rho, x) ≈ aggregate(tree_sig, x)
        end

        @testset "Cobb-Douglas at ρ = 0" begin
            tree = CESNode_ρ(0.0, (0.3, 0.7), (CESLeaf(:A), CESLeaf(:B)))
            x = [2.0, 3.0]
            expected = (2.0 / 0.3)^0.3 * (3.0 / 0.7)^0.7
            @test aggregate(tree, x) ≈ expected atol = 1e-12
        end

        @testset "Negative ρ gives complements" begin
            tree = CESNode_ρ(-1.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test tree.σ ≈ 0.5
        end

        @testset "Invalid ρ ≥ 1" begin
            @test_throws ArgumentError CESNode_ρ(1.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test_throws ArgumentError CESNode_ρ(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
        end
    end

    # ================================================================
    # Dimension checks
    # ================================================================
    @testset "Dimension mismatch" begin
        tree = CESNode(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
        @test_throws DimensionMismatch aggregate(tree, [1.0])
        @test_throws DimensionMismatch aggregate(tree, [1.0, 2.0, 3.0])
    end

    # ================================================================
    # Numeric type propagation
    # ================================================================
    @testset "Type propagation" begin
        @testset "Float32 inputs" begin
            tree = CESNode(2.0f0, (0.5f0, 0.5f0), (CESLeaf(:A), CESLeaf(:B)))
            Q = aggregate(tree, Float32[4.0, 9.0])
            @test Q isa Float32
            @test Q ≈ 6.25f0
        end

        @testset "BigFloat inputs" begin
            tree = CESNode(big"2.0", (big"0.5", big"0.5"), (CESLeaf(:A), CESLeaf(:B)))
            Q = aggregate(tree, [big"4.0", big"9.0"])
            @test Q isa BigFloat
            @test Q ≈ big"6.25"
        end
    end

    # ================================================================
    # Log-sum-exp method
    # ================================================================
    @testset "LSE method" begin
        @testset "Agrees with standard for general σ" begin
            tree = CESNode(2.0, (0.6, 0.4), (CESLeaf(:A), CESLeaf(:B)))
            x = [4.0, 9.0]
            @test aggregate(tree, x; method = :lse) ≈ aggregate(tree, x) atol = 1e-12
        end

        @testset "Agrees with standard for nested tree" begin
            inner = CESNode(4.0, (0.5, 0.5), (CESLeaf(:x), CESLeaf(:y)))
            tree = CESNode(1.5, (0.7, 0.3), (CESLeaf(:A), inner))
            x = [1.0, 2.0, 1.5]
            @test aggregate(tree, x; method = :lse) ≈ aggregate(tree, x) atol = 1e-10
        end

        @testset "Cobb-Douglas (σ = 1)" begin
            tree = CESNode(1.0, (0.3, 0.7), (CESLeaf(:A), CESLeaf(:B)))
            x = [2.0, 3.0]
            expected = (2.0 / 0.3)^0.3 * (3.0 / 0.7)^0.7
            @test aggregate(tree, x; method = :lse) ≈ expected atol = 1e-12
        end

        @testset "Leontief (σ = 0)" begin
            tree = CESNode(0.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test aggregate(tree, [1.0, 100.0]; method = :lse) ≈ 2.0
        end

        @testset "Linear (σ = Inf)" begin
            tree = CESNode(Inf, (0.6, 0.4), (CESLeaf(:A), CESLeaf(:B)))
            @test aggregate(tree, [3.0, 7.0]; method = :lse) ≈ 4.6
        end

        @testset "Homogeneity" begin
            tree = CESNode(2.0, (0.6, 0.4), (CESLeaf(:A), CESLeaf(:B)))
            x = [1.5, 2.3]
            λ = 3.7
            @test aggregate(tree, λ .* x; method = :lse) ≈ λ * aggregate(tree, x; method = :lse) atol = 1e-10
        end

        @testset "Near σ = 1: LSE is smoother" begin
            x = [2.0, 3.0]
            α = (0.4, 0.6)
            # Perturbation wider than √eps so one side falls outside the σ ≈ 1 branch
            σ_lo = 1.0 - 1e-4
            σ_hi = 1.0 + 1e-4

            tree_lo = CESNode(σ_lo, α, (CESLeaf(:A), CESLeaf(:B)))
            tree_hi = CESNode(σ_hi, α, (CESLeaf(:A), CESLeaf(:B)))

            lse_gap = abs(aggregate(tree_hi, x; method = :lse) - aggregate(tree_lo, x; method = :lse))
            std_gap = abs(aggregate(tree_hi, x) - aggregate(tree_lo, x))

            @test lse_gap ≤ std_gap + 1e-10
        end
    end

    # ================================================================
    # Invalid method keyword
    # ================================================================
    @testset "Invalid method" begin
        tree = CESNode(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
        @test_throws ArgumentError aggregate(tree, [1.0, 2.0]; method = :bogus)
    end

end
