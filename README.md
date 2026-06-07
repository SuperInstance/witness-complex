# witness-complex

> **Witness complex — approximate the topology of massive datasets using a small set of landmark points**

[![crates.io](https://img.shields.io/crates/v/witness-complex.svg)](https://crates.io/crates/witness-complex)
[![docs.rs](https://docs.rs/witness-complex/badge.svg)](https://docs.rs/witness-complex)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## What is the Witness Complex?

The Vietoris-Rips complex is the workhorse of topological data analysis — but it scales terribly. For n points, computing the full VR complex requires O(n³) distance computations for just 2-simplices, and the complex size grows combinatorially.

The **witness complex** (de Silva & Carlsson, 2004) solves this scalability problem. Instead of using all points, you select a small set of **landmark** points. The remaining points (the "witnesses") vote for which simplices should exist based on proximity to the landmarks. A witness "testifies" that a simplex exists when all its vertices are among the witness's nearest landmarks.

This reduces the complex from O(nᵏ) simplices (for k-dimensional VR) to O(mᵏ) where m ≪ n is the number of landmarks — a dramatic compression that preserves the essential topology.

## Why Does This Matter?

The witness complex makes topological data analysis practical for large datasets:

- **Scalability**: Analyze 100K+ point clouds by selecting 100–1000 landmarks instead of building the full complex
- **Noise robustness**: Witnesses "average out" noise by voting — spurious features don't get enough votes
- **Memory efficiency**: Store only the landmark complex, not the full pairwise distance matrix
- **Approximation quality**: Under mild conditions, the witness complex recovers the same homology as the full complex

Real-world applications:
- **Single-cell genomics**: 100K+ cells → 500 landmarks → tractable topological analysis of differentiation trajectories
- **Point cloud processing**: LiDAR scans with millions of points reduced to thousands of landmarks
- **Image collections**: Analyze the topology of large image datasets via landmark embeddings
- **Network analysis**: Study the shape of social networks, citation graphs, and neural connectivity

## Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                  Witness Complex Pipeline                     │
│                                                              │
│  Point Cloud (n points)                                      │
│  ┌──────────────────────┐                                    │
│  │  · · ·  ·  · · · ·  │     Landmark Selection              │
│  │  · · · · · · · · ·  │──── (random / maxmin / k-means++)   │
│  │  · · · · · · · · ·  │     ──────────────────────────────  │
│  │  · · ★ · · · · · ·  │              │                      │
│  │  · · · · · ★ · · ·  │              ▼                      │
│  │  · · · · · · · · ·  │     Landmarks (m ≪ n)               │
│  │  · ★ · · · · · · ·  │     ┌───┐                           │
│  │  · · · · · · · ★ ·  │     │ ★ │ ← Selected points         │
│  └──────────────────────┘     └───┘                           │
│                                      │                        │
│                    Witness Voting     ▼                        │
│                    ┌──────────────────────────┐               │
│  Weak Witness:     │ For each non-landmark:   │               │
│  Point votes for   │   find k nearest landmarks│               │
│  σ if all vertices │   vote for simplices from │               │
│  are among k       │   nearest landmark set    │               │
│  nearest landmarks │                           │               │
│                    └──────────┬───────────────┘               │
│                               ▼                               │
│                    ┌──────────────────────┐                    │
│                    │  Witness Complex     │                    │
│                    │  (on landmark set)   │                    │
│                    │  Vertices: m         │                    │
│                    │  Simplices: O(mᵏ)    │                    │
│                    └──────────────────────┘                    │
└──────────────────────────────────────────────────────────────┘
```

## Quick Start

```rust
use witness_complex::{LandmarkSelector, WitnessComplex};

// A dataset of 10 points in 1D
let points = vec![
    vec![0.0], vec![1.0], vec![2.0], vec![3.0], vec![4.0],
    vec![5.0], vec![6.0], vec![7.0], vec![8.0], vec![9.0],
];

// Select 4 landmarks using maxmin (spreads landmarks evenly)
let landmarks = LandmarkSelector::maxmin(&points, 4);
println!("Landmarks at indices: {:?}", landmarks);

// Build a weak witness complex (up to dimension 2)
let complex = WitnessComplex::build_weak(&points, &landmarks, 2);

println!("Vertices: {}", complex.num_vertices());
println!("Total simplices: {}", complex.num_simplices());
println!("Edges: {:?}", complex.simplices_of_dim(1).len());
```

### Landmark Selection Strategies

```rust
// Random selection (fast, reproducible with seed)
let lm = LandmarkSelector::random(&points, 5, 42);

// MaxMin: greedily maximizes minimum distance between landmarks
// Best for coverage — landmarks are spread as far apart as possible
let lm = LandmarkSelector::maxmin(&points, 5);

// K-means++: weighted random sampling proportional to distance²
// Good balance of speed and coverage
let lm = LandmarkSelector::kmeans_pp(&points, 5, 0);
```

### Strong Witness Complex

```rust
// Strong witness: stricter criterion — witness must certify the simplex
// by being closer to it than any other simplex of the same dimension
let strong = WitnessComplex::build_strong(&points, &landmarks, 2);
println!("Strong complex simplices: {}", strong.num_simplices());
```

### Querying the Complex

```rust
// Check if a specific simplex exists
if complex.has_simplex(&[landmarks[0], landmarks[1]]) {
    println!("Edge exists between landmarks {} and {}", landmarks[0], landmarks[1]);
}

// Get all simplices of a given dimension
let triangles = complex.simplices_of_dim(2);
println!("Found {} triangles", triangles.len());
```

## API Reference

### LandmarkSelector

| Method | Returns | Description |
|--------|---------|-------------|
| `LandmarkSelector::random(points, k, seed)` | `Vec<usize>` | Random selection with deterministic seed |
| `LandmarkSelector::maxmin(points, k)` | `Vec<usize>` | Greedy farthest-point sampling |
| `LandmarkSelector::kmeans_pp(points, k, seed)` | `Vec<usize>` | K-means++ initialization |

### Witness (internal)

| Method | Returns | Description |
|--------|---------|-------------|
| `Witness::new(points, idx, landmarks)` | `Witness` | Build witness with sorted landmark distances |
| `witness.nearest_landmarks(k)` | `Vec<usize>` | K nearest landmark indices |

### WitnessComplex

| Method | Returns | Description |
|--------|---------|-------------|
| `WitnessComplex::build_weak(points, landmarks, max_dim)` | `WitnessComplex` | Weak witness complex |
| `WitnessComplex::build_strong(points, landmarks, max_dim)` | `WitnessComplex` | Strong witness complex |
| `complex.num_vertices()` | `usize` | Number of landmark vertices |
| `complex.num_simplices()` | `usize` | Total simplices in the complex |
| `complex.simplices_of_dim(d)` | `Vec<&Vec<usize>>` | All d-dimensional simplices |
| `complex.has_simplex(vertices)` | `bool` | Check if simplex exists |

## Mathematical Background

### Weak Witness

A point w ∈ X is a **weak witness** for a simplex σ = {l₀, ..., lₖ} if the k+1 vertices of σ are among the k+1 nearest landmarks of w. This means:

```
w weakly witnesses σ ⟺ {l₀, ..., lₖ} ⊆ k-nearest-landmarks(w)
```

Weak witnesses are permissive — many witnesses can vote for the same simplex, providing robustness against noise.

### Strong Witness

A point w is a **strong witness** for σ if the (k+2)-th nearest landmark of w is strictly farther than the (k+1)-th. This ensures w "belongs" to σ more than to any competing simplex:

```
w strongly witnesses σ ⟺ d(w, σ) < d(w, σ') for all other σ' of same dimension
```

Strong witnesses produce smaller, more conservative complexes with higher topological fidelity.

### Landmark Selection

- **Random**: O(k) time, no distance computations needed. Good for quick exploration.
- **MaxMin**: O(n·k) time. Greedily selects points that maximize the minimum inter-landmark distance. Guarantees good coverage.
- **K-means++**: O(n·k) time. Samples proportional to squared distance from nearest existing landmark. Balances exploration and exploitation.

### Approximation Guarantees

Under the hypothesis that the landmark set is an ε-net of the point cloud (every point is within distance ε of some landmark), the witness complex approximates the Čech complex at scale 2ε, recovering the correct homology up to dimension (max_dim − 1).

## Installation

```bash
cargo add witness-complex
```

Or add to your `Cargo.toml`:

```toml
[dependencies]
witness-complex = "0.1.0"
```

## Related Crates

- [`cech-complex`](https://github.com/SuperInstance/cech-complex) — Čech complex (exact topology via ball intersection)
- [`mapper-graph`](https://github.com/SuperInstance/mapper-graph) — Mapper algorithm for topological summaries
- [`persistence-landscape`](https://github.com/SuperInstance/persistence-landscape) — Persistence landscapes for statistical TDA
- [`betti-curve`](https://github.com/SuperInstance/betti-curve) — Betti curves and Euler characteristic curves

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
