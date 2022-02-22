#[cfg(test)]
use super::ELEMENT_NUM;
#[cfg(test)]
use crate::bem::{bem_calc, calc_conjugate_values};
#[cfg(test)]
use crate::error::MyLinalgError;
#[cfg(test)]
use crate::subprograms::{circ_integral_trapez, exact_u_normal_dv, gauss, CircleCurve};
#[cfg(test)]
use ndarray::prelude::*;

#[test]
fn gauss_solve_sample() {
    let a = array![[1., 3., 1.], [1., 1., -1.], [3., 11., 5.]];
    let b = array![9., 1., 35.];
    let sol = gauss(&a, &b).unwrap();
    assert_eq!(sol, array![5., 0., 4.])
}

#[test]
fn circular_integral_sample() {
    #[allow(unused_variables)]
    let func = |x: usize| 1.;
    let curve = CircleCurve::new(array![0.0, 0.0], 1.0, ELEMENT_NUM as u32);
    let integral_val = circ_integral_trapez(func, &curve);
    let exact_val = 2. * std::f64::consts::PI;
    assert!((integral_val - exact_val).abs() < 1.0e-9)
}

#[test]
fn dot_test() {
    let v = array![2., 3.];
    let diff: f64 = v.dot(&v) - 13.;
    assert!(diff.abs() < 1.0e-10)
}

#[test]
fn bem_solve_accuracy() -> Result<(), MyLinalgError> {
    let curve = CircleCurve::new(array![0.0, 0.0], 1.0, ELEMENT_NUM as u32);
    let div_num = ELEMENT_NUM;
    let vectors = calc_conjugate_values(div_num, &curve);
    let q_vec = vectors.1?;
    let mut exact_q_vec = Array1::<f64>::zeros(div_num);
    // 厳密なq_vecを計算
    for (i, el) in exact_q_vec.iter_mut().enumerate() {
        *el = exact_u_normal_dv(&curve.points[i]);
    }
    let mut error = 0.0;
    for el in (q_vec - exact_q_vec).iter() {
        error += el * el;
    }
    println!("{}", error);
    Ok(())
}

#[test]
fn interior_point_calc_test() {
    let curve = CircleCurve::new(array![0.0, 0.0], 1.0, ELEMENT_NUM as u32);
    let p = array![0.5, 0.0];
    let val = bem_calc(ELEMENT_NUM, &p, &curve).unwrap();
    println!("{}", val);
}
