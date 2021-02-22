//! Sampling of a Poisson disk distribution in multiple dimensions.
//!
//! The Poisson disk distribution produces samples of which no two samples are too close
//! to each other. This results in a more uniform distribution than from pure random sampling.
//!
//! This library is an implementation of the algorithm introduced by Robert Bridson \[1\]
//! which is O(N) for producing N samples. That is, the sampling time increases linearly
//! with the number of produced samples. For two-dimensional sampling, the sampling time
//! increases with the area and for three-dimensional sampling with the volume.
//!
//! # Examples
//! ## Three dimensions
//! ```rust
//! use poisson_diskus::bridson;
//!
//! let box_size = [3.0, 5.0, 7.0];
//! let rmin = 0.5;
//! let num_attempts = 30;
//! let use_pbc = true;
//!
//! let coords = bridson(&box_size, rmin, num_attempts, use_pbc).unwrap();
//! ```
//!
//! ## Larger number of dimensions
//! ```rust
//! use poisson_diskus::bridson;
//!
//! let box_size = [3.0, 5.0, 3.0, 2.0, 1.0];
//! let rmin = 1.0;
//! let num_attempts = 30;
//! let use_pbc = true;
//!
//! let coords = bridson(&box_size, rmin, num_attempts, use_pbc).unwrap();
//!
//! for coord in coords {
//!     assert_eq!(coord.len(), box_size.len());
//! }
//! ```
//!
//! # Citations
//! \[1\] Bridson, R. (2007). Fast Poisson disk sampling in arbitrary dimensions. SIGGRAPH sketches, 10, 1.
//!
mod bridson;
mod coord;
mod error;
mod grid;
mod sample;

pub use bridson::{bridson, bridson_rng};
pub use error::Error;
