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

            node = CES(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test node.σ == 2.0
            @test node.α == (0.5, 0.5)
        end

        @testset "Type promotion" begin
            # Int σ with Float64 α → Float64
            node = CES(2, (0.5, 0.5), (:A, :B))
            @test node.σ isa Float64

            # All Int → Int (if someone really wants that)
            node_int = CES(2, (1, 1), (:A, :B))
            @test node_int.σ isa Int

            # BigFloat
            node_big = CES(big"2.0", (big"0.5", big"0.5"), (:A, :B))
            @test node_big.σ isa BigFloat
            @test eltype(node_big.α) == BigFloat
        end

        @testset "Heterogeneous α promotion" begin
            # Mixed Float64/Int α — promotes to Float64
            node = CES(2.0, (0.5, 1), (:A, :B))
            @test node.σ isa Float64
            @test node.α isa NTuple{2, Float64}
            @test node.α == (0.5, 1.0)

            # Int σ with mixed-numeric α
            node2 = CES(2, (0.5, 1), (:A, :B))
            @test node2.σ isa Float64
            @test node2.α === (0.5, 1.0)

            # All-Int α with Float σ — already covered indirectly but worth pinning
            node3 = CES(2.0, (1, 1), (:A, :B))
            @test node3.σ isa Float64
            @test node3.α === (1.0, 1.0)

            # Three-element heterogeneous
            node4 = CES(2.0, (1, 0.5, big"0.25"), (:A, :B, :C))
            @test node4.σ isa BigFloat
            @test eltype(node4.α) == BigFloat
        end

        @testset "Vector/iterable convenience" begin
            node = CES(2.0, [0.5, 0.5], [:A, :B])
            @test node.α == (0.5, 0.5)
        end

        @testset "Anonymous leaves" begin
            leaf = CESLeaf()
            @test leaf.name == Symbol()

            node = CES(2.0, (0.5, 0.5), (CESLeaf(), CESLeaf()))
            @test aggregate(node, [4.0, 9.0]) ≈ 6.25
        end

        @testset "Leaf-free construction" begin
            # Two-arg form: leaves created automatically
            node = CES(2.0, (0.5, 0.5))
            @test aggregate(node, [4.0, 9.0]) ≈ 6.25

            # Works with vectors too
            node_v = CES(2.0, [0.5, 0.5])
            @test aggregate(node_v, [4.0, 9.0]) ≈ 6.25

            # Three inputs
            node3 = CES(3.0, (0.4, 0.35, 0.25))
            @test aggregate(node3, [2.0, 3.0, 4.0]) > 0
        end

        @testset "Leaf-free nesting" begin
            # Inner node with anonymous leaves, outer with named + nested
            inner = CES(4.0, (0.6, 0.4))
            tree = CES(1.5, (0.7, 0.3), (:dom, inner))
            @test aggregate(tree, [1.0, 1.0, 1.0]) > 0
        end

        @testset "Symbol shorthand for leaves" begin
            # Tuple of symbols is normalized to tuple of named CESLeaf
            node = CES(2.0, (0.5, 0.5), (:steel, :aluminum))
            @test node.children[1] isa CESLeaf
            @test node.children[1].name === :steel
            @test node.children[2].name === :aluminum

            # Equivalent to explicit construction
            explicit = CES(2.0, (0.5, 0.5), (CESLeaf(:steel), CESLeaf(:aluminum)))
            @test node.children == explicit.children
            @test typeof(node) == typeof(explicit)

            # Mixed: symbol + CESLeaf + nested CES
            inner = CES(3.0, (0.6, 0.4))
            mixed = CES(1.5, (0.5, 0.3, 0.2), (:dom, CESLeaf(:imp), inner))
            @test mixed.children[1].name === :dom
            @test mixed.children[2].name === :imp
            @test mixed.children[3] === inner

            # Strings also work
            from_str = CES(2.0, (0.5, 0.5), ("a", "b"))
            @test from_str.children[1].name === :a

            # Aggregation works through normalized children
            @test aggregate(node, [4.0, 9.0]) ≈ aggregate(explicit, [4.0, 9.0])

            # Leaf names recoverable via leaf_names
            @test leaf_names(node) == [:steel, :aluminum]

            # Invalid child type rejected at construction
            @test_throws ArgumentError CES(2.0, (0.5, 0.5), (1, 2))
        end

        @testset "CESNode back-compat alias" begin
            # CESNode is retained as an alias for CES
            @test CESNode === CES
            @test CESNode_ρ === CES_ρ

            # Old-style construction still works
            old = CESNode(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            via_ces = CES(2.0, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test typeof(old) == typeof(via_ces)
            @test aggregate(old, [4.0, 9.0]) ≈ aggregate(via_ces, [4.0, 9.0])

            # ρ-constructor alias
            old_ρ = CESNode_ρ(0.5, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            via_ces_ρ = CES_ρ(0.5, (0.5, 0.5), (CESLeaf(:A), CESLeaf(:B)))
            @test old_ρ.σ == via_ces_ρ.σ
        end

        @testset "Invalid construction" begin
            @test_throws ArgumentError CES(-1.0, (0.5, 0.5), (:A, :B))
            @test_throws ArgumentError CES(2.0, (0.5, -0.1), (:A, :B))
            @test_throws ArgumentError CES(2.0, (0.5,), (:A, :B))
        end

        @testset "Limiting σ values accepted" begin
            leon = CES(0.0, (0.5, 0.5), (:A, :B))
            @test leon.σ == 0.0
            lin = CES(Inf, (0.5, 0.5), (:A, :B))
            @test isinf(lin.σ)
        end
    end

    # ================================================================
    # Tree introspection
    # ================================================================
    @testset "Tree introspection" begin
        inner = CES(3.0, (0.6, 0.4), (:US, :EU))
        tree = CES(1.5, (0.7, 0.3), (:dom, inner))

        @test leaf_names(tree) == [:dom, :US, :EU]
        @test leaf_names(CESLeaf(:solo)) == [:solo]

        sub = CES(5.0, (0.5, 0.5), (:v1, :v2))
        mid = CES(3.0, (0.4, 0.6), (sub, :C))
        top = CES(1.5, (0.5, 0.5), (:A, mid))
        @test leaf_names(top) == [:A, :v1, :v2, :C]
    end

    # ================================================================
    # Display
    # ================================================================
    @testset "Display" begin
        # Compact show
        tree = CES(1.5, (0.7, 0.3), (:dom,
            CES(4.0, (0.6, 0.4), (:US, :EU))))
        str = sprint(show, tree)
        @test occursin("σ=1.5", str)
        @test occursin("2 children", str)
        @test !occursin("└", str)

        # Leaf display
        @test sprint(show, CESLeaf(:A)) == "CESLeaf(:A)"
        @test sprint(show, CESLeaf()) == "CESLeaf()"

        # Named node compact show
        named = CES(2.0, (0.5, 0.5); name=:test)
        @test occursin("test: CES", sprint(show, named))

        # Unnamed compact show
        @test startswith(sprint(show, CES(2.0, (0.5, 0.5))), "CES")

        # show_tree: full tree display
        tree_str = sprint(show_tree, tree)
        @test occursin("σ=1.5", tree_str)
        @test occursin("σ=4.0", tree_str)
        @test occursin("dom", tree_str)
        @test occursin("US", tree_str)
        @test occursin("EU", tree_str)
        @test occursin("└", tree_str)

        # show_tree: anonymous leaves show bullet
        anon_tree = CES(2.0, (0.5, 0.5))
        anon_str = sprint(show_tree, anon_tree)
        @test occursin("•", anon_str)

        # show_tree: named nodes
        named_tree = CES(1.5, (0.7, 0.3), (:dom,
            CES(4.0, (0.6, 0.4), (:US, :EU); name=:imports));
            name=:total)
        named_str = sprint(show_tree, named_tree)
        @test occursin("total: CES", named_str)
        @test occursin("imports: CES", named_str)
    end

    # ================================================================
    # General CES
    # ================================================================
    @testset "General CES" begin
        @testset "Equal weights, equal inputs" begin
            tree = CES(2.0, (0.5, 0.5), (:A, :B))
            @test aggregate(tree, [1.0, 1.0]) ≈ 1.0
        end

        @testset "Known analytic value (σ = 2)" begin
            tree = CES(2.0, (0.5, 0.5), (:A, :B))
            # σ=2 → ρ=1/2, Q = [0.5*4^(1/2) + 0.5*9^(1/2)]^2 = 6.25
            @test aggregate(tree, [4.0, 9.0]) ≈ 6.25
        end

        @testset "Three inputs" begin
            tree = CES(3.0, (0.4, 0.35, 0.25),
                (:A, :B, :C))
            Q = aggregate(tree, [2.0, 3.0, 4.0])
            @test Q > 0
        end
    end

    # ================================================================
    # Limiting cases
    # ================================================================
    @testset "Leontief (σ = 0)" begin
        tree = CES(0.0, (0.5, 0.5), (:A, :B))
        @test aggregate(tree, [1.0, 100.0]) ≈ 2.0
        @test aggregate(tree, [10.0, 3.0]) ≈ 6.0

        # Asymmetric weights
        tree2 = CES(0.0, (0.25, 0.75), (:A, :B))
        # min(2/0.25, 6/0.75) = min(8, 8) = 8
        @test aggregate(tree2, [2.0, 6.0]) ≈ 8.0
    end

    @testset "Cobb-Douglas (σ = 1)" begin
        tree = CES(1.0, (0.3, 0.7), (:A, :B))
        x = [2.0, 3.0]
        expected = (2.0 / 0.3)^0.3 * (3.0 / 0.7)^0.7
        @test aggregate(tree, x) ≈ expected atol = 1e-12

        # Shares should work as exponents
        tree_eq = CES(1.0, (0.5, 0.5), (:A, :B))
        @test aggregate(tree_eq, [4.0, 4.0]) ≈ (4.0 / 0.5)^0.5 * (4.0 / 0.5)^0.5
    end

    @testset "Linear / perfect substitutes (σ = Inf)" begin
        tree = CES(Inf, (0.6, 0.4), (:A, :B))
        @test aggregate(tree, [3.0, 7.0]) ≈ 4.6

        # Only one input matters
        @test aggregate(tree, [10.0, 0.0]) ≈ 6.0
    end

    # ================================================================
    # Nesting
    # ================================================================
    @testset "Two-level nesting" begin
        inner = CES(3.0, (0.6, 0.4), (:US, :EU))
        tree = CES(1.5, (0.7, 0.3), (:dom, inner))
        Q = aggregate(tree, [1.0, 1.0, 1.0])
        @test Q > 0
    end

    @testset "Three-level nesting" begin
        sub = CES(5.0, (0.5, 0.5), (:v1, :v2))
        mid = CES(3.0, (0.4, 0.6), (sub, :C))
        top = CES(1.5, (0.5, 0.5), (:A, mid))
        Q = aggregate(top, [1.0, 1.1, 0.95, 1.3])
        @test Q > 0
    end

    @testset "Mixed limiting cases in nesting" begin
        # Leontief at top, high σ at bottom
        inner = CES(10.0, (0.5, 0.5), (:A, :B))
        tree = CES(0.0, (0.5, 0.5), (:C, inner))
        Q = aggregate(tree, [4.0, 3.0, 5.0])
        @test Q > 0

        # Linear at top, Cobb-Douglas at bottom
        inner_cd = CES(1.0, (0.5, 0.5), (:A, :B))
        tree_lin = CES(Inf, (0.6, 0.4), (inner_cd, :C))
        Q2 = aggregate(tree_lin, [2.0, 3.0, 5.0])
        @test Q2 > 0
    end

    # ================================================================
    # Homogeneity: Q(λx) = λQ(x)
    # ================================================================
    @testset "Homogeneity" begin
        λ = 3.7

        @testset "General CES" begin
            tree = CES(2.0, (0.6, 0.4), (:A, :B))
            x = [1.5, 2.3]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x) atol = 1e-10
        end

        @testset "Cobb-Douglas" begin
            tree = CES(1.0, (0.4, 0.6), (:A, :B))
            x = [2.0, 5.0]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x) atol = 1e-10
        end

        @testset "Leontief" begin
            tree = CES(0.0, (0.3, 0.7), (:A, :B))
            x = [2.0, 5.0]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x)
        end

        @testset "Linear" begin
            tree = CES(Inf, (0.6, 0.4), (:A, :B))
            x = [3.0, 7.0]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x)
        end

        @testset "Nested" begin
            inner = CES(4.0, (0.5, 0.5), (:x, :y))
            tree = CES(1.5, (0.7, 0.3), (:A, inner))
            x = [1.0, 2.0, 1.5]
            @test aggregate(tree, λ .* x) ≈ λ * aggregate(tree, x) atol = 1e-10
        end
    end

    # ================================================================
    # Monotonicity: increasing any input increases Q
    # ================================================================
    @testset "Monotonicity" begin
        tree = CES(2.0, (0.5, 0.5), (:A, :B))
        x = [3.0, 5.0]
        Q_base = aggregate(tree, x)
        @test aggregate(tree, [3.1, 5.0]) > Q_base
        @test aggregate(tree, [3.0, 5.1]) > Q_base
    end

    # ================================================================
    # CES_ρ
    # ================================================================
    @testset "CES_ρ" begin
        @testset "Equivalence with CES" begin
            tree_rho = CES_ρ(0.5, (0.6, 0.4), (:A, :B))
            tree_sig = CES(2.0, (0.6, 0.4), (:A, :B))
            x = [3.0, 7.0]
            @test aggregate(tree_rho, x) ≈ aggregate(tree_sig, x)
        end

        @testset "Cobb-Douglas at ρ = 0" begin
            tree = CES_ρ(0.0, (0.3, 0.7), (:A, :B))
            x = [2.0, 3.0]
            expected = (2.0 / 0.3)^0.3 * (3.0 / 0.7)^0.7
            @test aggregate(tree, x) ≈ expected atol = 1e-12
        end

        @testset "Negative ρ gives complements" begin
            tree = CES_ρ(-1.0, (0.5, 0.5), (:A, :B))
            @test tree.σ ≈ 0.5
        end

        @testset "Invalid ρ ≥ 1" begin
            @test_throws ArgumentError CES_ρ(1.0, (0.5, 0.5), (:A, :B))
            @test_throws ArgumentError CES_ρ(2.0, (0.5, 0.5), (:A, :B))
        end
    end

    # ================================================================
    # Dimension checks
    # ================================================================
    @testset "Dimension mismatch" begin
        tree = CES(2.0, (0.5, 0.5), (:A, :B))
        @test_throws DimensionMismatch aggregate(tree, [1.0])
        @test_throws DimensionMismatch aggregate(tree, [1.0, 2.0, 3.0])
    end

    # ================================================================
    # Numeric type propagation
    # ================================================================
    @testset "Type propagation" begin
        @testset "Float32 inputs" begin
            tree = CES(2.0f0, (0.5f0, 0.5f0), (:A, :B))
            Q = aggregate(tree, Float32[4.0, 9.0])
            @test Q isa Float32
            @test Q ≈ 6.25f0
        end

        @testset "BigFloat inputs" begin
            tree = CES(big"2.0", (big"0.5", big"0.5"), (:A, :B))
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
            tree = CES(2.0, (0.6, 0.4), (:A, :B))
            x = [4.0, 9.0]
            @test aggregate(tree, x; method = :lse) ≈ aggregate(tree, x) atol = 1e-12
        end

        @testset "Agrees with standard for nested tree" begin
            inner = CES(4.0, (0.5, 0.5), (:x, :y))
            tree = CES(1.5, (0.7, 0.3), (:A, inner))
            x = [1.0, 2.0, 1.5]
            @test aggregate(tree, x; method = :lse) ≈ aggregate(tree, x) atol = 1e-10
        end

        @testset "Cobb-Douglas (σ = 1)" begin
            tree = CES(1.0, (0.3, 0.7), (:A, :B))
            x = [2.0, 3.0]
            expected = (2.0 / 0.3)^0.3 * (3.0 / 0.7)^0.7
            @test aggregate(tree, x; method = :lse) ≈ expected atol = 1e-12
        end

        @testset "Leontief (σ = 0)" begin
            tree = CES(0.0, (0.5, 0.5), (:A, :B))
            @test aggregate(tree, [1.0, 100.0]; method = :lse) ≈ 2.0
        end

        @testset "Linear (σ = Inf)" begin
            tree = CES(Inf, (0.6, 0.4), (:A, :B))
            @test aggregate(tree, [3.0, 7.0]; method = :lse) ≈ 4.6
        end

        @testset "Homogeneity" begin
            tree = CES(2.0, (0.6, 0.4), (:A, :B))
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

            tree_lo = CES(σ_lo, α, (:A, :B))
            tree_hi = CES(σ_hi, α, (:A, :B))

            lse_gap = abs(aggregate(tree_hi, x; method = :lse) - aggregate(tree_lo, x; method = :lse))
            std_gap = abs(aggregate(tree_hi, x) - aggregate(tree_lo, x))

            @test lse_gap ≤ std_gap + 1e-10
        end
    end

    # ================================================================
    # Invalid method keyword
    # ================================================================
    @testset "Invalid method" begin
        tree = CES(2.0, (0.5, 0.5), (:A, :B))
        @test_throws ArgumentError aggregate(tree, [1.0, 2.0]; method = :bogus)
    end

    # ================================================================
    # Type stability
    # ================================================================
    @testset "Type stability" begin

        # ── Flat trees ─────────────────────────────────────────
    
        @testset "Flat tree — standard" begin
            tree = CES(2.0, (0.4, 0.6))
            x = [3.0, 5.0]
            @test @inferred(aggregate(tree, x)) isa Float64
        end
    
        @testset "Flat tree — lse" begin
            tree = CES(2.0, (0.4, 0.6))
            x = [3.0, 5.0]
            @test @inferred(aggregate(tree, x; method = :lse)) isa Float64
        end
    
        # ── Nested tree ────────────────────────────────────────
    
        @testset "Nested tree — standard" begin
            inner = CES(4.0, (0.3, 0.7), (:A, :B))
            tree = CES(1.5, (0.5, 0.5), (:C, inner))
            x = [1.0, 2.0, 3.0]
            @test @inferred(aggregate(tree, x)) isa Float64
        end
    
        @testset "Nested tree — lse" begin
            inner = CES(4.0, (0.3, 0.7), (:A, :B))
            tree = CES(1.5, (0.5, 0.5), (:C, inner))
            x = [1.0, 2.0, 3.0]
            @test @inferred(aggregate(tree, x; method = :lse)) isa Float64
        end
    
        # ── Three levels deep ─────────────────────────────────
    
        @testset "Three-level nesting" begin
            bottom = CES(3.0, (0.6, 0.4), (:a, :b))
            mid    = CES(2.0, (0.5, 0.5), (bottom, :c))
            top    = CES(1.5, (0.7, 0.3), (:d, mid))
            x = [1.0, 2.0, 3.0, 4.0]
            @test @inferred(aggregate(top, x)) isa Float64
        end
    
        # ── Limiting cases ────────────────────────────────────
    
        @testset "Leontief (σ = 0)" begin
            tree = CES(0.0, (0.4, 0.6))
            x = [3.0, 5.0]
            @test @inferred(aggregate(tree, x)) isa Float64
        end
    
        @testset "Cobb-Douglas (σ = 1)" begin
            tree = CES(1.0, (0.4, 0.6))
            x = [3.0, 5.0]
            @test @inferred(aggregate(tree, x)) isa Float64
        end
    
        @testset "Linear (σ = Inf)" begin
            tree = CES(Inf, (0.4, 0.6))
            x = [3.0, 5.0]
            @test @inferred(aggregate(tree, x)) isa Float64
        end
    
        # ── Limiting cases nested ─────────────────────────────
    
        @testset "Leontief nested in general CES" begin
            inner = CES(0.0, (0.5, 0.5), (:A, :B))
            tree = CES(2.0, (0.6, 0.4), (:C, inner))
            x = [1.0, 2.0, 3.0]
            @test @inferred(aggregate(tree, x)) isa Float64
        end
    
        # ── 3-ary node ────────────────────────────────────────
    
        @testset "3-ary node" begin
            tree = CES(2.5, (0.3, 0.3, 0.4))
            x = [1.0, 2.0, 3.0]
            @test @inferred(aggregate(tree, x)) isa Float64
        end
    
        # ── Kernel-level stability ────────────────────────────
    
        @testset "Kernels" begin
            α = (0.4, 0.6)
            xv = (3.0, 5.0)
            @test @inferred(Armington._ces(2.0, α, xv))         isa Float64
            @test @inferred(Armington._ces(0.0, α, xv))         isa Float64
            @test @inferred(Armington._ces(1.0, α, xv))         isa Float64
            @test @inferred(Armington._ces(Inf, α, xv))         isa Float64
            @test @inferred(Armington._ces_lse(2.0, α, xv))     isa Float64
            @test @inferred(Armington._ces_lse(1.0, α, xv))     isa Float64
            @test @inferred(Armington._leontief(α, xv))         isa Float64
            @test @inferred(Armington._linear(α, xv))           isa Float64
            @test @inferred(Armington._cobb_douglas(α, xv))     isa Float64
        end
    
        # ── Internal recursion ────────────────────────────────
    
        @testset "_compute and _collect_children" begin
            inner = CES(4.0, (0.3, 0.7), (:A, :B))
            tree = CES(1.5, (0.5, 0.5), (:C, inner))
            x = [1.0, 2.0, 3.0]
    
            # _compute on a leaf returns (Float64, Int)
            @test @inferred(Armington._compute(CESLeaf(:A), x, 1, :standard)) isa Tuple{Float64, Int}
    
            # _compute on a node returns (Float64, Int)
            @test @inferred(Armington._compute(tree, x, 1, :standard)) isa Tuple{Float64, Int}
    
            # _collect_children on empty tuple
            @test @inferred(Armington._collect_children((), x, 1, :standard)) isa Tuple{Tuple{}, Int}
    
            # _collect_children on the tree's children
            @test @inferred(Armington._collect_children(tree.children, x, 1, :standard)) isa Tuple{<:Tuple, Int}
        end
    
        # ── BigFloat propagation ──────────────────────────────
    
        @testset "BigFloat stability" begin
            tree = CES(big"2.0", (big"0.5", big"0.5"))
            x = [big"4.0", big"9.0"]
            @test @inferred(aggregate(tree, x)) isa BigFloat
        end

    end
end
