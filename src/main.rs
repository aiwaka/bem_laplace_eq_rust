mod bem;
mod error;
mod output;
mod subprograms;
mod test;
mod type_and_trait;

use bem::bem_calc;
use ndarray::prelude::*;
use type_and_trait::Point;

// 境界要素数. 精度にかかわる.
const ELEMENT_NUM: usize = 32;

fn main() {
    const CALC_N: i32 = 8;
    let curve = subprograms::CircleCurve::new(array![0.0, 0.0], 1.0, ELEMENT_NUM as u32);
    let curve_center = &curve.center;
    let curve_radius = curve.radius;
    // 点列コンテナ. 点, 計算値の順に格納
    let mut points_list = Vec::<(Point<f64>, f64)>::with_capacity(32 * 32);
    // 値を計算する点列をつくる
    for i in 0..=2 * CALC_N {
        for j in 0..=2 * CALC_N {
            let x1 = curve_radius * ((i - CALC_N) as f64) / (CALC_N as f64);
            let x2 = curve_radius * ((j - CALC_N) as f64) / (CALC_N as f64);
            if x1.powf(2.) + x2.powf(2.) >= curve_radius.powf(2.) {
                continue;
            }
            let new_point = array![x1, x2] + curve_center.clone();
            points_list.push((new_point, 0.));
        }
    }
    // 点を一つずつ持ってきてその点での値を計算.
    for p_data in points_list.iter_mut() {
        let point = &p_data.0;
        let res = bem_calc(ELEMENT_NUM, point, &curve);
        match res {
            Ok(val) => {
                // Okなら取り出した値を格納
                (*p_data).1 = val;
            }
            Err(e) => println!("{:?}", e),
        };
    }
    // 出力
    match output::output_data(&points_list) {
        Ok(_) => (),
        Err(e) => println!("{}", &e),
    };
}
