module Armington

export CESLeaf, CESNode, aggregate, leaf_names, CESNode_ρ, show_tree

# ─────────────────────────────────────────────────────────────
# Types
# ─────────────────────────────────────────────────────────────

"""
    CESLeaf(name::Symbol)

Terminal node in a CES nesting tree. Represents a single input.
The `name` field is for labeling only; it does not affect computation.
"""
struct CESLeaf
    name::Symbol
end

CESLeaf(name::AbstractString) = CESLeaf(Symbol(name))
CESLeaf() = CESLeaf(Symbol())

"""
    CESNode(σ, α, children)

Interior node in a CES nesting tree.

# Arguments
- `σ`: Elasticity of substitution (σ ≥ 0). Special cases:
        σ = 0    → Leontief (perfect complements)
        σ = 1    → Cobb-Douglas
        σ = Inf  → Linear (perfect substitutes)
- `α`: Distribution parameters. A `Tuple` of length `N`, one per child.
        These enter the CES directly: U = [Σ αᵢ xᵢ^ρ]^(1/ρ) with ρ = (σ-1)/σ.
- `children`: `Tuple` of `CESNode` or `CESLeaf`, one per distribution parameter.

If `children` is omitted, anonymous leaves are created automatically
(one per element of α):

    CESNode(2.0, (0.5, 0.5))  # equivalent to CESNode(2.0, (0.5, 0.5), (CESLeaf(), CESLeaf()))

The type parameter `T` is the numeric type of σ and α. This is inferred
automatically and allows the tree parameters to carry AD dual numbers,
arbitrary-precision types, etc.

The nesting structure is encoded in the type, so Julia compiles specialized code
for each tree topology.
"""
struct CESNode{T<:Real, N, C<:Tuple}
    name::Symbol
    σ::T
    α::NTuple{N, T}
    children::C

    function CESNode(σ::T, α::NTuple{N, T}, children::C; name::Symbol = Symbol()) where {T<:Real, N, C<:Tuple}
        length(children) == N || throw(ArgumentError(
            "Length of α ($(N)) must equal number of children ($(length(children)))."
        ))
        σ >= 0 || throw(ArgumentError("Elasticity σ must be non-negative, got $σ."))
        all(a -> a > 0, α) || throw(ArgumentError("All distribution parameters α must be positive."))
        new{T, N, C}(name, σ, α, children)
    end
end

# Promote σ and α to a common type
function CESNode(σ::Real, α::NTuple{N, Real}, children::Tuple; name::Symbol = Symbol()) where {N}
    T = promote_type(typeof(σ), eltype(α))
    CESNode(convert(T, σ), convert(NTuple{N, T}, α), children; name)
end

# Convenience: accept vectors / non-tuple iterables
function CESNode(σ::Real, α, children; name::Symbol = Symbol())
    CESNode(σ, Tuple(α), Tuple(children); name)
end

# Leaf-free: create anonymous leaves from α length
function CESNode(σ::Real, α::Union{Tuple, AbstractVector}; name::Symbol = Symbol())
    leaves = ntuple(_ -> CESLeaf(), length(α))
    CESNode(σ, α, leaves; name)
end

const CESTree = Union{CESLeaf, CESNode}

"""
    CESNode_ρ(ρ, α, children)

Construct a `CESNode` parametrized by ρ = (σ-1)/σ instead of σ.

The mapping is σ = 1/(1-ρ), with ρ ∈ (-∞, 1):
- ρ → -∞  ⟹  σ → 0   (Leontief)
- ρ = 0   ⟹  σ = 1   (Cobb-Douglas)
- ρ → 1   ⟹  σ → ∞   (Linear)

The node stores σ internally; this is purely a convenience for construction.
"""
function CESNode_ρ(ρ::Real, α, children; name::Symbol = Symbol())
    ρ < 1 || throw(ArgumentError("ρ must be < 1, got $ρ."))
    σ = 1 / (1 - ρ)
    CESNode(σ, α, children; name)
end

# ─────────────────────────────────────────────────────────────
# Tree introspection
# ─────────────────────────────────────────────────────────────

"""
    leaf_names(tree) -> Vector{Symbol}

Return leaf names in depth-first order. This is the order in which
a flat input vector is expected by `aggregate`.
"""
leaf_names(leaf::CESLeaf) = [leaf.name]
leaf_names(node::CESNode) = mapreduce(leaf_names, vcat, node.children)

"""Number of leaves in the tree."""
_nleaves(::CESLeaf) = 1
_nleaves(node::CESNode) = sum(_nleaves, node.children)

# ─────────────────────────────────────────────────────────────
# Display
# ─────────────────────────────────────────────────────────────

function Base.show(io::IO, leaf::CESLeaf)
    print(io, "CESLeaf(:", leaf.name, ")")
end

function Base.show(io::IO, node::CESNode{T,N}) where {T,N}
    label = node.name == Symbol() ? "" : string(node.name, ": ")
    print(io, label, "CESNode(σ=", node.σ, ", α=", node.α, ", ", N, " children)")
end

"""
    show_tree(tree)
    show_tree(io, tree)

Print the full nesting structure of a CES tree.
"""
show_tree(tree::CESTree) = show_tree(stdout, tree)
show_tree(io::IO, tree::CESTree) = _show_tree(io, tree, 0, "")

function _show_tree(io::IO, leaf::CESLeaf, depth::Int, prefix::String)
    label = leaf.name == Symbol() ? "•" : string(leaf.name)
    println(io, "  "^depth, prefix, label)
end

function _show_tree(io::IO, node::CESNode{T,N}, depth::Int, prefix::String) where {T,N}
    label = node.name == Symbol() ? "" : string(node.name, ": ")
    println(io, "  "^depth, prefix, label, "CESNode(σ=", node.σ, ", α=", node.α, ")")
    for i in 1:N
        _show_tree(io, node.children[i], depth + 1, "└ ")
    end
end

# ─────────────────────────────────────────────────────────────
# Core: aggregate
# ─────────────────────────────────────────────────────────────

"""
    aggregate(tree, x; method = :standard)

Recursively compute the CES aggregate for `tree`, given a flat vector `x`
whose entries correspond to leaves in depth-first order.

For an interior node with elasticity σ and distribution parameters α:

- σ = 0:   U = min(xᵢ / αᵢ)               (Leontief)
- σ = 1:   U = Π (xᵢ / αᵢ)^αᵢ             (Cobb-Douglas)
- σ = Inf: U = Σ αᵢ xᵢ                     (Linear)
- else:    U = [Σ αᵢ xᵢ^ρ]^(1/ρ)  where ρ = (σ-1)/σ

# Keyword arguments
- `method`: `:standard` (default) uses the direct CES formula with explicit
  branches for σ = 0, σ ≈ 1, and σ = Inf. `:lse` uses a log-sum-exp
  formulation that is numerically smoother near σ = 1, at the cost of
  additional transcendental function evaluations. Use `:lse` when σ is
  a parameter being estimated and may pass through 1.

Returns a scalar whose type is determined by promotion of the tree
parameter type and the element type of `x`.
"""
function aggregate(tree::CESTree, x::AbstractVector; method::Symbol = :standard)
    length(x) == _nleaves(tree) || throw(DimensionMismatch(
        "Expected $(_nleaves(tree)) inputs, got $(length(x))."
    ))
    method ∈ (:standard, :lse) || throw(ArgumentError(
        "method must be :standard or :lse, got :$method."
    ))
    _aggregate(tree, x, Ref(1), method)
end

# ── Internal recursion ─────────────────────────────────────

function _aggregate(leaf::CESLeaf, x::AbstractVector, idx::Ref{Int}, method::Symbol)
    val = x[idx[]]
    idx[] += 1
    return val
end

function _aggregate(node::CESNode{T, N}, x::AbstractVector, idx::Ref{Int}, method::Symbol) where {T, N}
    child_vals = ntuple(N) do i
        _aggregate(node.children[i], x, idx, method)
    end
    kernel = method === :lse ? _ces_lse : _ces
    return kernel(node.σ, node.α, child_vals)
end

# ── Regimes ────────────────────────────────────────────────

_leontief(α::NTuple{N}, x::NTuple{N}) where {N} = minimum(i -> x[i] / α[i], 1:N)
_linear(α::NTuple{N}, x::NTuple{N}) where {N} = sum(i -> α[i] * x[i], 1:N)
_cobb_douglas(α::NTuple{N}, x::NTuple{N}) where {N} = prod(i -> (x[i] / α[i])^α[i], 1:N)

# ── Kernels ────────────────────────────────────────────────

"""
Standard CES kernel with explicit branches for limiting cases.
"""
function _ces(σ, α::NTuple{N}, x::NTuple{N}) where {N}
    if σ == 0
        return _leontief(α, x)
    elseif isinf(σ)
        return _linear(α, x)
    elseif σ ≈ 1
        return _cobb_douglas(α, x)
    else
        ρ = (σ - one(σ)) / σ
        s = sum(i -> α[i] * x[i]^ρ, 1:N)
        return s^(one(s) / ρ)
    end
end

"""
Log-sum-exp CES kernel. Numerically smoother near σ = 1 (ρ ≈ 0).

Uses the identity U = exp((1/ρ) · logsumexp(log αᵢ + ρ log xᵢ)) with
a stable logsumexp implementation. Falls back to exact Cobb-Douglas
at σ = 1 and retains explicit branches for σ = 0 and σ = Inf.
"""
function _ces_lse(σ, α::NTuple{N}, x::NTuple{N}) where {N}
    if σ == 0
        return _leontief(α, x)
    elseif isinf(σ)
        return _linear(α, x)
    elseif σ == 1
        return _cobb_douglas(α, x)
    else
        ρ = (σ - one(σ)) / σ
        z = ntuple(i -> log(α[i]) + ρ * log(x[i]), N)
        z_max = maximum(z)
        lse = z_max + log(sum(i -> exp(z[i] - z_max), 1:N))
        return exp(lse / ρ)
    end
end

end # module
