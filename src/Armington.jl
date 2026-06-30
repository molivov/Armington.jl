module Armington

export CESLeaf, CES, CESNode, aggregate, leaf_names, CES_ρ, CESNode_ρ, show_tree

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
    CES(σ, α, children)

Interior node in a CES nesting tree.

# Arguments
- `σ`: Elasticity of substitution (σ ≥ 0). Special cases:
        σ = 0    → Leontief
        σ = 1    → Cobb-Douglas
        σ = Inf  → Linear
- `α`: Distribution parameters. A `Tuple` of length `N`, one per child.
        These enter the CES directly: U = (Σ αᵢ xᵢ^ρ)^(1/ρ) with ρ = (σ-1)/σ.
- `children`: `Tuple` of `CES`, `CESLeaf`, `Symbol`, or `AbstractString` (one per
        distribution parameter). Symbols and strings are converted to named
        `CESLeaf` instances, so `CES(σ, (0.5, 0.5), (:steel, :aluminum))` is
        shorthand for `CES(σ, (0.5, 0.5), (CESLeaf(:steel), CESLeaf(:aluminum)))`.

If `children` is omitted, anonymous leaves are created automatically
(one per element of α):

    CES(2.0, (0.5, 0.5))  # equivalent to CES(2.0, (0.5, 0.5), (CESLeaf(), CESLeaf()))

"""
struct CES{T<:Real, N, C<:Tuple}
    name::Symbol
    σ::T
    α::NTuple{N, T}
    children::C

    function CES(σ::T, α::NTuple{N, T}, children::Tuple; name::Symbol = Symbol()) where {T<:Real, N}
        normalized = map(_to_leaf, children)
        length(normalized) == N || throw(ArgumentError(
            "Length of α ($(N)) must equal number of children ($(length(normalized)))."
        ))
        σ >= 0 || throw(ArgumentError("Elasticity σ must be non-negative, got $σ."))
        all(a -> a > 0, α) || throw(ArgumentError("All distribution parameters α must be positive."))
        C = typeof(normalized)
        new{T, N, C}(name, σ, α, normalized)
    end
end

# Promote σ and α to a common type
function CES(σ::Real, α::NTuple{N, Real}, children::Tuple; name::Symbol = Symbol()) where {N}
    T = promote_type(typeof(σ), typeof.(α)...)
    CES(convert(T, σ), convert(NTuple{N, T}, α), children; name)
end

# Convenience: accept vectors / non-tuple iterables
function CES(σ::Real, α, children; name::Symbol = Symbol())
    CES(σ, Tuple(α), Tuple(children); name)
end

# Leaf-free: create anonymous leaves from α length
function CES(σ::Real, α::Union{Tuple, AbstractVector}; name::Symbol = Symbol())
    leaves = ntuple(_ -> CESLeaf(), length(α))
    CES(σ, α, leaves; name)
end

# Normalize a child input to a tree element. Symbols and strings are
# promoted to named leaves; existing leaves and nodes pass through.
# Anything else throws with an informative message; this is what lets
# `CES(σ, α, (:steel, :aluminum))` work as shorthand for
# `CES(σ, α, (CESLeaf(:steel), CESLeaf(:aluminum)))`.
_to_leaf(s::Symbol) = CESLeaf(s)
_to_leaf(s::AbstractString) = CESLeaf(s)
_to_leaf(l::CESLeaf) = l
_to_leaf(n::CES) = n
_to_leaf(x) = throw(ArgumentError(
    "Children must be CESLeaf, CES, Symbol, or AbstractString, got $(typeof(x))."
))

const CESTree = Union{CESLeaf, CES}

"""
    CESNode

Alias for [`CES`](@ref). Retained for compatibility with code written
against earlier versions of the package; `CES` is the canonical name.
"""
const CESNode = CES

"""
    CES_ρ(ρ, α, children)

Construct a `CES` parametrized by ρ = (σ-1)/σ instead of σ.

The mapping is σ = 1/(1-ρ), with ρ ∈ (-∞, 1):
- ρ → -∞  ⟹  σ → 0   (Leontief)
- ρ = 0   ⟹  σ = 1   (Cobb-Douglas)
- ρ → 1   ⟹  σ → ∞   (Linear)

The node stores σ internally; this is purely a convenience for construction.
"""
function CES_ρ(ρ::Real, α, children; name::Symbol = Symbol())
    ρ < 1 || throw(ArgumentError("ρ must be < 1, got $ρ."))
    σ = 1 / (1 - ρ)
    CES(σ, α, children; name)
end

"""
    CESNode_ρ

Alias for [`CES_ρ`](@ref). Retained for compatibility with code written
against earlier versions of the package; `CES_ρ` is the canonical name.
"""
const CESNode_ρ = CES_ρ

# ─────────────────────────────────────────────────────────────
# Tree introspection
# ─────────────────────────────────────────────────────────────

"""
    leaf_names(tree) -> Vector{Symbol}

Return leaf names in depth-first order. This is the order in which
a flat input vector is expected by `aggregate`.
"""
leaf_names(leaf::CESLeaf) = [leaf.name]
leaf_names(node::CES) = mapreduce(leaf_names, vcat, node.children)

"""Number of leaves in the tree."""
_nleaves(::CESLeaf) = 1
_nleaves(node::CES) = sum(_nleaves, node.children)

# ─────────────────────────────────────────────────────────────
# Display
# ─────────────────────────────────────────────────────────────

function Base.show(io::IO, leaf::CESLeaf)
    if leaf.name == Symbol()
        print(io, "CESLeaf()")
    else
        print(io, "CESLeaf(:", leaf.name, ")")
    end
end

function Base.show(io::IO, node::CES{T,N}) where {T,N}
    label = node.name == Symbol() ? "" : string(node.name, ": ")
    print(io, label, "CES(σ=", node.σ, ", α=", node.α, ", ", N, " children)")
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

function _show_tree(io::IO, node::CES{T,N}, depth::Int, prefix::String) where {T,N}
    label = node.name == Symbol() ? "" : string(node.name, ": ")
    println(io, "  "^depth, prefix, label, "CES(σ=", node.σ, ", α=", node.α, ")")
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
whose entries correspond to leaves in depth-first order. This is the
order from `show_tree`.

For an interior node with elasticity σ and distribution parameters α:

- σ = 0:   U = min(xᵢ / αᵢ)               (Leontief)
- σ = 1:   U = Π (xᵢ / αᵢ)^αᵢ             (Cobb-Douglas)
- σ = Inf: U = Σ αᵢ xᵢ                    (Linear)
- else:    U = (Σ αᵢ xᵢ^ρ)^(1/ρ)  where ρ = (σ-1)/σ

# Keyword arguments
- `method`: `:standard` (default) uses the direct CES formula with explicit
  branches for σ = 0, σ ≈ 1, and σ = Inf. `:lse` uses a log-sum-exp
  formulation that is numerically smoother near σ = 1, at the cost of
  additional transcendental function evaluations. Use `:lse` when σ is
  a parameter being estimated and may pass through 1.

"""
function aggregate(tree::CESTree, x::AbstractVector; method::Symbol = :standard)
    length(x) == _nleaves(tree) || throw(DimensionMismatch(
        "Expected $(_nleaves(tree)) inputs, got $(length(x))."
    ))
    method ∈ (:standard, :lse) || throw(ArgumentError(
        "method must be :standard or :lse, got :$method."
    ))
    val, _ = _compute(tree, x, 1, method)
    return val
end

# ── Internal recursion ─────────────────────────────────────
#
# Every function returns (value, next_index). The index is a plain
# Int passed by value.

function _compute(::CESLeaf, x::AbstractVector, idx::Int, ::Symbol)
    return x[idx], idx + 1
end

function _compute(node::CES{T,N}, x::AbstractVector, idx::Int, method::Symbol) where {T,N}
    child_vals, idx = _collect_children(node.children, x, idx, method)
    kernel = method === :lse ? _ces_lse : _ces
    return kernel(node.σ, node.α, child_vals), idx
end

# ── Recursive tuple peeling ───────────────────────────────
#
# Processes children left to right. Each call peels the first
# element with children[1] and shrinks the rest with Base.tail.

function _collect_children(::Tuple{}, ::AbstractVector, idx::Int, ::Symbol)
    return (), idx
end

function _collect_children(children::Tuple, x::AbstractVector, idx::Int, method::Symbol)
    val, idx = _compute(children[1], x, idx, method)
    rest, idx = _collect_children(Base.tail(children), x, idx, method)
    return (val, rest...), idx
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
