# Armington.jl

Recursive CES (constant elasticity of substitution) aggregation over arbitrary nesting trees.

Build nested CES structures, also called Armington aggregators (utility functions, production functions) and evaluate them in a single call. The nesting topology is encoded in Julia's type system, so each tree compiles to specialized, allocation-free code.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/molivov/Armington.jl")
```

## Quick start

```julia
using Armington

# Two-level Armington: domestic vs. imports (US, EU)
imports = CESNode(4.0, (0.6, 0.4), (CESLeaf(:US), CESLeaf(:EU)))
tree = CESNode(1.5, (0.7, 0.3), (CESLeaf(:dom), imports))

leaf_names(tree)           # [:dom, :US, :EU]
aggregate(tree, [10, 5, 8]) # CES composite quantity
```

## API

### Types

**`CESLeaf(name::Symbol)`**: terminal node (a single good/input). `CESLeaf()` creates an anonymous leaf.

**`CESNode(σ, α, children)`**: interior CES node.
- `σ`: elasticity of substitution (σ ≥ 0, including `Inf`)
- `α`: distribution parameters (tuple, one per child)
- `children`: tuple of `CESNode` or `CESLeaf`

If `children` is omitted, anonymous leaves are created from the length of `α`:

```julia
CESNode(2.0, (0.5, 0.5))  # two anonymous leaves
```

### Functions

**`aggregate(tree, x; method = :standard)`**: recursively compute the CES composite given a flat input vector `x` in depth-first leaf order. `method = :lse` uses a log-sum-exp formulation that is numerically smoother near σ = 1, useful when estimating σ.

**`leaf_names(tree)`**: return leaf names in the order expected by `aggregate`.

**`CESNode_ρ(ρ, α, children)`**: construct a `CESNode` using ρ = (σ−1)/σ instead of σ. Useful when you think in terms of the substitution parameter directly.

### Limiting cases

| σ | ρ | Regime | Formula |
|---|---|--------|---------|
| 0 | −∞ | Leontief | Q = min(xᵢ / αᵢ) |
| 1 | 0 | Cobb-Douglas | Q = Π (xᵢ / αᵢ)^αᵢ |
| ∞ | 1 | Linear | Q = Σ αᵢ xᵢ |

Pass `σ = 0`, `σ = 1` or `σ = Inf` for exact limiting cases. The general CES formula is used for all other values.

## Three-level example

```julia
# Varieties within each origin
us = CESNode(8.0, (0.5, 0.5), (CESLeaf(:US_east), CESLeaf(:US_west)))
eu = CESNode(8.0, (0.6, 0.4), (CESLeaf(:EU_north), CESLeaf(:EU_south)))

# Origins within imports
imports = CESNode(4.0, (0.55, 0.45), (us, eu))

# Domestic vs. import composite
tree = CESNode(1.5, (0.7, 0.3), (CESLeaf(:domestic), imports))

aggregate(tree, [20.0, 5.0, 4.0, 6.0, 3.0])
```

```julia
show_tree(tree)

total: CESNode(σ=1.5, α=(0.7, 0.3))
  └ domestic
  └ imports: CESNode(σ=4.0, α=(0.55, 0.45))
    └ us: CESNode(σ=8.0, α=(0.5, 0.5))
      └ US_east
      └ US_west
    └ eu: CESNode(σ=8.0, α=(0.6, 0.4))
      └ EU_north
      └ EU_south
```
