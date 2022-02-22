use ndarray::prelude::*;
use num_traits::{Float, FromPrimitive, NumAssignOps};

#[allow(dead_code)]
pub type Point<T> = Array<T, Ix1>;
pub type Vector<T> = Array<T, Ix1>;
pub type Matrix<T> = Array<T, Ix2>;

pub trait MyNum: Clone + std::fmt::Display + Float + FromPrimitive + NumAssignOps {}
impl<S> MyNum for S where S: Clone + std::fmt::Display + Float + FromPrimitive + NumAssignOps {}
