# witness-complex

> **Witness complex for landmark-based topological approximation**

[![crates.io](https://img.shields.io/crates/v/witness-complex.svg)](https://crates.io/crates/witness-complex)
[![license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

The witness complex approximates the topology of a large point cloud using a small set of landmarks. Instead of computing the full Vietoris-Rips complex (expensive), witnesses 'vote' for simplices based on their proximity to landmarks.

## Algorithm

1. **Select landmarks** from the point cloud (random, maxmin, or k-means++)
2. **Witness voting**: Each non-landmark point acts as a witness
3. **Weak witness**: Point votes for simplex σ if all vertices of σ are among its k nearest landmarks
4. **Strong witness**: Stricter condition requiring the point to be closer to σ than any other simplex
5. **Build complex**: Add simplices that receive enough witness votes

Much more scalable than VR complex for large datasets.

## Installation

```toml
[dependencies]
witness-complex = "0.1.0"
```

## License

MIT © [SuperInstance](https://github.com/SuperInstance)

---

*Part of the [Exocortex](https://github.com/SuperInstance/exocortex) project — persistent cognitive substrate for multi-agent systems.*
