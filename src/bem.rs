use super::type_and_trait::Point;
use crate::error::MyLinalgError;
use crate::subprograms::{
    circ_integral_trapez, exact_u, fund_gamma, fund_gamma_normal_dv, gauss, CircleCurve,
};
use ndarray::prelude::*;

/// 曲線構造体データからDirichlet条件に基づいて共役な境界量を求める.
/// 用いた条件のベクトルとともに返す.
pub fn calc_conjugate_values(
    div_num: usize,
    curve: &CircleCurve,
) -> (Point<f64>, Result<Point<f64>, MyLinalgError>) {
    let mut u_mat = Array2::<f64>::zeros((div_num, div_num));
    let mut w_mat = Array2::<f64>::zeros((div_num, div_num));
    let mut u_vec = Array1::<f64>::zeros(div_num);
    // U, W行列の要素を計算.
    // iter_mut()は列方向の順番に要素への参照をつくるイテレータ
    for (i, el) in u_mat.iter_mut().enumerate() {
        let m = i / div_num;
        let n = i % div_num;
        *el = curve.u_components(m, n);
    }
    for (i, el) in w_mat.iter_mut().enumerate() {
        let m = i / div_num;
        let n = i % div_num;
        *el = curve.w_components(m, n);
    }
    // Dirichlet境界条件としてuベクトルの値を計算.
    for (i, el) in u_vec.iter_mut().enumerate() {
        *el = exact_u(&curve.points[i]);
    }
    let b_vec = w_mat.dot(&u_vec);
    (u_vec, gauss(&u_mat, &b_vec))
}

pub fn bem_calc(
    div_num: usize,
    calc_point: &Point<f64>,
    curve: &CircleCurve,
) -> Result<f64, MyLinalgError> {
    // let mut u_mat = Array2::<f64>::zeros((div_num, div_num));
    // let mut w_mat = Array2::<f64>::zeros((div_num, div_num));
    // let mut u_vec = Array1::<f64>::zeros(div_num);
    // // U, W行列の要素を計算.
    // // iter_mut()は列方向の順番に要素への参照をつくるイテレータ
    // for (i, el) in u_mat.iter_mut().enumerate() {
    //     let m = i / div_num;
    //     let n = i % div_num;
    //     *el = curve.u_components(m, n);
    // }
    // for (i, el) in w_mat.iter_mut().enumerate() {
    //     let m = i / div_num;
    //     let n = i % div_num;
    //     *el = curve.w_components(m, n);
    // }
    // // Dirichlet境界条件としてuベクトルの値を計算.
    // for (i, el) in u_vec.iter_mut().enumerate() {
    //     *el = exact_u(&curve.points[i]);
    // }
    // let b_vec = w_mat.dot(&u_vec);
    // let q_vec = gauss(&u_mat, &b_vec)?;
    let vectors = calc_conjugate_values(div_num, curve);
    let u_vec = vectors.0;
    let q_vec = vectors.1?;

    // 計算用のクロージャを作成（calc_pointはxに束縛）
    let x = calc_point;
    let func = |i: usize| {
        fund_gamma(x, &curve.points[i]) * q_vec[i]
            - fund_gamma_normal_dv(x, &curve.points[i]) * u_vec[i]
    };
    let result = circ_integral_trapez(func, curve);
    Ok(result)
}
