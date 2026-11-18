//! Gravitas - Gravity field modeling for astrodynamics
//!
//! Rust reimplementation of the gravitas library using EGM96 model,
//! parallelized with rayon for multi-threaded computation.

pub mod egm96;
pub mod gravity;

pub mod python;

pub use egm96::EGM96;
pub use gravity::earth_acceleration;
