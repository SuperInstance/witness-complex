//! # witness-complex
//!
//! **Witness complex** construction for landmark-based topological data analysis.
//!
//! Given a set of *landmark* points selected from a point cloud, the witness
//! complex is a simplicial complex built using the remaining *witness* points
//! that "vote" for simplices based on proximity to landmarks.
//!
//! # Example
//!
//! ```
//! use witness_complex::{LandmarkSelector, WitnessComplex};
//!
//! let points = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
//! let landmarks = LandmarkSelector::random(&points, 3, 42);
//! let complex = WitnessComplex::build_weak(&points, &landmarks, 2);
//! assert!(complex.num_vertices() > 0);
//! ```

use std::collections::HashSet;

// ---------------------------------------------------------------------------
// LandmarkSelector
// ---------------------------------------------------------------------------

/// Strategies for selecting landmark points from a dataset.
pub struct LandmarkSelector;

impl LandmarkSelector {
    /// Select `k` landmarks uniformly at random using seed `seed`.
    pub fn random(points: &[Vec<f64>], k: usize, seed: u64) -> Vec<usize> {
        let n = points.len();
        if k >= n {
            return (0..n).collect();
        }
        // Simple LCG pseudo-random number generator
        let mut rng = seed.wrapping_add(1);
        let mut selected = HashSet::new();
        while selected.len() < k {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let idx = (rng >> 33) as usize % n;
            selected.insert(idx);
        }
        let mut result: Vec<usize> = selected.into_iter().collect();
        result.sort();
        result
    }

    /// MaxMin selection: greedily pick landmarks that maximise the minimum
    /// distance to already-chosen landmarks.
    pub fn maxmin(points: &[Vec<f64>], k: usize) -> Vec<usize> {
        let n = points.len();
        if k >= n {
            return (0..n).collect();
        }
        if k == 0 {
            return vec![];
        }

        // Start with the first point
        let mut landmarks = vec![0usize];
        let mut min_dists: Vec<f64> = vec![f64::INFINITY; n];

        #[allow(clippy::needless_range_loop)]
        for j in 0..n {
            min_dists[j] = euclidean(&points[0], &points[j]);
        }

        while landmarks.len() < k {
            // Find point with largest minimum distance
            let (next, _) = min_dists
                .iter()
                .enumerate()
                .filter(|(i, _)| !landmarks.contains(i))
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, &d)| (i, d))
                .unwrap_or((0, 0.0));

            landmarks.push(next);
            for j in 0..n {
                let d = euclidean(&points[next], &points[j]);
                if d < min_dists[j] {
                    min_dists[j] = d;
                }
            }
        }

        landmarks.sort();
        landmarks
    }

    /// K-means++ landmark initialisation: pick first landmark randomly (seed),
    /// then sample proportional to squared distance to nearest landmark.
    pub fn kmeans_pp(points: &[Vec<f64>], k: usize, seed: u64) -> Vec<usize> {
        let n = points.len();
        if k >= n {
            return (0..n).collect();
        }
        if k == 0 {
            return vec![];
        }

        let mut rng = seed.wrapping_add(1);
        let next_random = |rng: &mut u64| -> f64 {
            *rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (*rng >> 33) as f64 / (1u64 << 31) as f64
        };

        // First landmark from seed
        let first = (seed as usize) % n;
        let mut landmarks = vec![first];
        let mut min_dists: Vec<f64> = vec![0.0; n];
        for j in 0..n {
            min_dists[j] = euclidean(&points[first], &points[j]).powi(2);
        }

        while landmarks.len() < k {
            let total: f64 = min_dists.iter().enumerate()
                .filter(|(i, _)| !landmarks.contains(i))
                .map(|(_, d)| *d)
                .sum();

            if total <= 0.0 {
                // All remaining points coincide with landmarks
                break;
            }

            let r = next_random(&mut rng) * total;
            let mut cumsum = 0.0;
            let mut next = 0usize;
            #[allow(clippy::needless_range_loop)]
            for j in 0..n {
                if landmarks.contains(&j) {
                    continue;
                }
                cumsum += min_dists[j];
                if cumsum >= r {
                    next = j;
                    break;
                }
            }

            landmarks.push(next);
            for j in 0..n {
                let d = euclidean(&points[next], &points[j]).powi(2);
                if d < min_dists[j] {
                    min_dists[j] = d;
                }
            }
        }

        landmarks.sort();
        landmarks
    }
}

/// Euclidean distance between two points.
fn euclidean(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

// ---------------------------------------------------------------------------
// Witness
// ---------------------------------------------------------------------------

/// A witness point with its nearest landmark distances.
#[derive(Debug, Clone)]
pub struct Witness {
    /// Index of the witness in the original point cloud.
    pub index: usize,
    /// Distances to all landmarks, sorted ascending.
    pub landmark_dists: Vec<(usize, f64)>,
}

impl Witness {
    /// Create a witness from a point and the set of landmarks.
    pub fn new(points: &[Vec<f64>], point_idx: usize, landmarks: &[usize]) -> Self {
        let mut dists: Vec<(usize, f64)> = landmarks
            .iter()
            .map(|&li| (li, euclidean(&points[point_idx], &points[li])))
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        Self {
            index: point_idx,
            landmark_dists: dists,
        }
    }

    /// The k nearest landmarks (landmark indices).
    pub fn nearest_landmarks(&self, k: usize) -> Vec<usize> {
        self.landmark_dists.iter().take(k).map(|(li, _)| *li).collect()
    }
}

// ---------------------------------------------------------------------------
// WeakWitness
// ---------------------------------------------------------------------------

/// A weak witness votes for a simplex σ if the witness's nearest landmarks
/// include all vertices of σ.
pub struct WeakWitness;

impl WeakWitness {
    /// Check if witness `w` votes for the simplex given by `landmark_indices`.
    pub fn votes_for(w: &Witness, landmark_indices: &[usize], nu: usize) -> bool {
        // A weak witness votes for σ if all vertices of σ are among the nu
        // nearest landmarks of w.
        let nearest: HashSet<usize> = w.nearest_landmarks(nu).into_iter().collect();
        landmark_indices.iter().all(|li| nearest.contains(li))
    }
}

// ---------------------------------------------------------------------------
// StrongWitness
// ---------------------------------------------------------------------------

/// A strong witness certifies a simplex if the witness is closer to the
/// affine hull of the simplex vertices than to any landmark not in the simplex.
pub struct StrongWitness;

impl StrongWitness {
    /// Simplified strong witness: the witness certifies the simplex if all
    /// vertices of the simplex are among its `k` nearest landmarks and the
    /// (k+1)-th nearest landmark is strictly farther.
    pub fn certifies(w: &Witness, landmark_indices: &[usize]) -> bool {
        let k = landmark_indices.len();
        if k > w.landmark_dists.len() {
            return false;
        }
        let nearest_k: HashSet<usize> = w.nearest_landmarks(k).into_iter().collect();
        if !landmark_indices.iter().all(|li| nearest_k.contains(li)) {
            return false;
        }
        // Check that (k+1)-th nearest is strictly farther
        if let Some((_, d_next)) = w.landmark_dists.get(k) {
            if let Some((_, d_k)) = w.landmark_dists.get(k - 1) {
                return *d_next > *d_k + 1e-12;
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// WitnessComplex
// ---------------------------------------------------------------------------

/// A simplicial complex built from witness relations.
#[derive(Debug, Clone)]
pub struct WitnessComplex {
    /// Vertices (landmark indices).
    vertices: Vec<usize>,
    /// Simplices stored as sorted Vec<usize> of landmark indices.
    simplices: HashSet<Vec<usize>>,
}

impl WitnessComplex {
    /// Build a weak witness complex of dimension up to `max_dim`.
    ///
    /// A witness votes for a k-simplex if all k+1 landmarks are among the
    /// witness's (k+1) nearest landmarks.
    pub fn build_weak(
        points: &[Vec<f64>],
        landmarks: &[usize],
        max_dim: usize,
    ) -> Self {
        let landmark_set: HashSet<usize> = landmarks.iter().copied().collect();
        let witnesses: Vec<Witness> = (0..points.len())
            .filter(|i| !landmark_set.contains(i))
            .map(|i| Witness::new(points, i, landmarks))
            .collect();

        let mut simplices: HashSet<Vec<usize>> = HashSet::new();

        // Add 0-simplices (vertices)
        for &li in landmarks {
            simplices.insert(vec![li]);
        }

        // For each witness, add simplices from its nearest landmarks
        for w in &witnesses {
            for dim in 1..=max_dim {
                if dim + 1 > w.landmark_dists.len() {
                    break;
                }
                let nearest: Vec<usize> = w.nearest_landmarks(dim + 1);
                let mut sigma = nearest;
                sigma.sort();
                simplices.insert(sigma);
            }
        }

        Self {
            vertices: landmarks.to_vec(),
            simplices,
        }
    }

    /// Build a strong witness complex.
    pub fn build_strong(
        points: &[Vec<f64>],
        landmarks: &[usize],
        max_dim: usize,
    ) -> Self {
        let landmark_set: HashSet<usize> = landmarks.iter().copied().collect();
        let witnesses: Vec<Witness> = (0..points.len())
            .filter(|i| !landmark_set.contains(i))
            .map(|i| Witness::new(points, i, landmarks))
            .collect();

        let mut simplices: HashSet<Vec<usize>> = HashSet::new();
        for &li in landmarks {
            simplices.insert(vec![li]);
        }

        // Generate candidate simplices and certify them
        for w in &witnesses {
            let nearest = w.nearest_landmarks(max_dim + 1);
            for dim in 1..=max_dim {
                if dim + 1 > nearest.len() {
                    break;
                }
                let mut sigma: Vec<usize> = nearest[..=dim].to_vec();
                sigma.sort();
                if StrongWitness::certifies(w, &sigma) {
                    simplices.insert(sigma);
                }
            }
        }

        Self {
            vertices: landmarks.to_vec(),
            simplices,
        }
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Total number of simplices.
    pub fn num_simplices(&self) -> usize {
        self.simplices.len()
    }

    /// Get all simplices of a given dimension.
    pub fn simplices_of_dim(&self, dim: usize) -> Vec<&Vec<usize>> {
        self.simplices
            .iter()
            .filter(|s| s.len() == dim + 1)
            .collect()
    }

    /// Check if a specific simplex exists.
    pub fn has_simplex(&self, vertices: &[usize]) -> bool {
        let mut v = vertices.to_vec();
        v.sort();
        self.simplices.contains(&v)
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_points() -> Vec<Vec<f64>> {
        vec![
            vec![0.0], vec![1.0], vec![2.0], vec![3.0], vec![4.0],
            vec![5.0], vec![6.0], vec![7.0], vec![8.0], vec![9.0],
        ]
    }

    #[test]
    fn test_random_landmarks() {
        let pts = sample_points();
        let lm = LandmarkSelector::random(&pts, 3, 42);
        assert_eq!(lm.len(), 3);
        assert!(lm.iter().all(|&i| i < pts.len()));
    }

    #[test]
    fn test_random_landmarks_deterministic() {
        let pts = sample_points();
        let l1 = LandmarkSelector::random(&pts, 3, 42);
        let l2 = LandmarkSelector::random(&pts, 3, 42);
        assert_eq!(l1, l2);
    }

    #[test]
    fn test_maxmin_landmarks() {
        let pts = sample_points();
        let lm = LandmarkSelector::maxmin(&pts, 3);
        assert_eq!(lm.len(), 3);
        // For 10 equally spaced points, maxmin should spread them
        assert!(lm[2] - lm[0] >= 5);
    }

    #[test]
    fn test_kmeans_pp_landmarks() {
        let pts = sample_points();
        let lm = LandmarkSelector::kmeans_pp(&pts, 4, 7);
        assert_eq!(lm.len(), 4);
    }

    #[test]
    fn test_witness_nearest() {
        let pts = sample_points();
        let landmarks = vec![0usize, 5, 9];
        let w = Witness::new(&pts, 2, &landmarks);
        let nearest = w.nearest_landmarks(1);
        assert_eq!(nearest, vec![0]); // point 2 is closest to landmark 0
    }

    #[test]
    fn test_witness_dists_sorted() {
        let pts = sample_points();
        let landmarks = vec![0usize, 5, 9];
        let w = Witness::new(&pts, 7, &landmarks);
        for i in 1..w.landmark_dists.len() {
            assert!(w.landmark_dists[i].1 >= w.landmark_dists[i - 1].1);
        }
    }

    #[test]
    fn test_weak_witness_votes() {
        let pts = sample_points();
        let landmarks = vec![0usize, 5];
        let w = Witness::new(&pts, 1, &landmarks);
        assert!(WeakWitness::votes_for(&w, &[0], 1));
    }

    #[test]
    fn test_weak_complex_has_vertices() {
        let pts = sample_points();
        let lm = vec![0usize, 5, 9];
        let wc = WitnessComplex::build_weak(&pts, &lm, 1);
        assert!(wc.has_simplex(&[0]));
        assert!(wc.has_simplex(&[5]));
        assert!(wc.has_simplex(&[9]));
    }

    #[test]
    fn test_weak_complex_edges() {
        let pts = vec![
            vec![0.0, 0.0], vec![0.1, 0.0], vec![1.0, 0.0], vec![1.1, 0.0],
        ];
        let lm = vec![0usize, 2];
        let wc = WitnessComplex::build_weak(&pts, &lm, 1);
        // Points 1 and 3 are witnesses; 1 is close to both 0 and 2
        assert!(wc.num_simplices() >= 3); // 2 vertices + at least 1 edge
    }

    #[test]
    fn test_strong_complex() {
        let pts = sample_points();
        let lm = vec![0usize, 5, 9];
        let wc = WitnessComplex::build_strong(&pts, &lm, 1);
        assert_eq!(wc.num_vertices(), 3);
    }

    #[test]
    fn test_simplices_of_dim() {
        let pts = sample_points();
        let lm = vec![0usize, 3, 6, 9];
        let wc = WitnessComplex::build_weak(&pts, &lm, 2);
        let dim0 = wc.simplices_of_dim(0);
        assert_eq!(dim0.len(), 4);
    }

    #[test]
    fn test_more_landmarks_than_points() {
        let pts = vec![vec![0.0], vec![1.0]];
        let lm = LandmarkSelector::random(&pts, 5, 1);
        assert_eq!(lm.len(), 2);
    }
}
