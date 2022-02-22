use crate::error::MyLinalgError;
use crate::type_and_trait::{Matrix, Point, Vector};
use ndarray::prelude::*;
use num::Float;

pub struct CircleCurve {
    pub center: Point<f64>,
    pub radius: f64,
    pub div_num: u32,
    pub points: Vec<Point<f64>>,
}
impl CircleCurve {
    /// 円周曲線を作成し, div_num個に分割する点列も生成する.
    pub fn new(center: Point<f64>, radius: f64, div_num: u32) -> Self {
        let mut curve = CircleCurve {
            center,
            radius,
            div_num,
            points: Vec::<Point<f64>>::with_capacity(div_num as usize),
        };
        curve.make_points();
        curve
    }

    /// カーブを指定して点列（境界要素）を作成する.
    /// これを使って積分も行う.
    /// vecを初期化したのちに格納する.
    /// カーブは円に限る.
    fn make_points(&mut self) {
        // TAUは2*piと等価
        let tau: f64 = std::f64::consts::TAU;
        self.points.clear();
        for k in 0..self.div_num {
            let theta = tau * (k as f64) / (self.div_num as f64);
            let x1 = theta.cos();
            let x2 = theta.sin();
            // unit_p: 単位円上の点
            let unit_p = array![x1, x2];
            let p = self.radius * unit_p + self.center.clone();
            self.points.push(p);
        }
    }

    /// 区間番号（0-index）m, nに対応する要素の値を返す.
    pub fn u_components(&self, m: usize, n: usize) -> f64 {
        let tau: f64 = std::f64::consts::TAU;
        let mid_index_2 = (m + 1_usize) % (self.div_num as usize);
        let y_index = (n + 1_usize) % (self.div_num as usize);
        let mid = (self.points[m].clone() + self.points[mid_index_2].clone()) / 2.;
        let p_x = self.points[n].clone();
        let p_y = self.points[y_index].clone();
        let comp_val = ComponentValues::new(&mid, &p_x, &p_y);

        if m == n {
            (1. - Float::ln(comp_val.h / 2.)) * comp_val.h / tau
        } else {
            let l_x1 = comp_val.l_x1;
            let l_x2 = comp_val.l_x2;
            let l_y1 = comp_val.l_y1;
            let r1 = comp_val.r1;
            let r2 = comp_val.r2;
            let h = comp_val.h;
            let theta = comp_val.theta;
            (l_x2 * Float::ln(r2) - l_x1 * Float::ln(r1) + h - l_y1 * theta) / tau
        }
    }

    pub fn w_components(&self, m: usize, n: usize) -> f64 {
        let tau: f64 = std::f64::consts::TAU;
        let mid_index_2 = (m + 1_usize) % (self.div_num as usize);
        let y_index = (n + 1_usize) % (self.div_num as usize);
        let mid = (self.points[m].clone() + self.points[mid_index_2].clone()) / 2.;
        let p_x = self.points[n].clone();
        let p_y = self.points[y_index].clone();
        if m == n {
            return 0.5_f64;
        }
        let comp_val = ComponentValues::new(&mid, &p_x, &p_y);
        comp_val.theta / tau
    }
}

/// 境界要素の計算に用いる値
#[derive(Default)]
pub struct ComponentValues {
    pub l_x1: f64,
    pub l_x2: f64,
    pub l_y1: f64,
    pub l_y2: f64,
    pub r1: f64,
    pub r2: f64,
    pub h: f64,
    pub theta: f64,
}
impl ComponentValues {
    /// mid: m番目の区間の中点.
    /// x, y: n番目の区間の端点.
    fn new(mid: &Point<f64>, x: &Point<f64>, y: &Point<f64>) -> Self {
        let norm = |v: &Point<f64>| Float::sqrt(v.dot(v));
        let h = norm(&(x - y));
        let t_vec = (y - x) / h;
        let n_vec = array![x[1] - y[1], y[0] - x[0]] / h;
        let l_x1 = (mid - x).dot(&t_vec);
        let l_x2 = (mid - y).dot(&t_vec);
        let l_y1 = (mid - x).dot(&n_vec);
        let l_y2 = (mid - y).dot(&n_vec);
        let r1 = norm(&(mid - x));
        let r2 = norm(&(mid - y));
        let theta = l_y2.atan2(l_x2) - l_y1.atan2(l_x1);
        ComponentValues {
            l_x1,
            l_x2,
            l_y1,
            l_y2,
            r1,
            r2,
            h,
            theta,
        }
    }
}

/// 消去法の係数行列と定数ベクトルのi,j行を入れ替える.
fn _gauss_swap(a: &mut ArrayViewMut2<f64>, b: &mut ArrayViewMut1<f64>, i: usize, j: usize) {
    let row_i = a.slice(s![i, ..]).to_owned();
    let row_j = a.slice(s![j, ..]).to_owned();
    a.slice_mut(s![i, ..]).assign(&row_j);
    a.slice_mut(s![j, ..]).assign(&row_i);
    let b_own = b.to_owned();
    b.slice_mut(s![i]).assign(&b_own.slice(s![j]));
    b.slice_mut(s![j]).assign(&b_own.slice(s![i]));
}

pub fn gauss(a: &Matrix<f64>, b: &Vector<f64>) -> Result<Vector<f64>, MyLinalgError> {
    const EPSILON: f64 = 1.0e-10;
    let mut acopy = a.clone();
    let mut bcopy = b.clone();
    let size = a.nrows();
    // 正方行列でなければErrを返す.
    if size != a.ncols() {
        return Err(MyLinalgError::NonSquareMatrix);
    }
    for k in 0..size {
        // pivoting
        let mut p = k;
        let mut pmax = acopy[[k, k]].abs();
        for i in (k + 1)..size {
            let el_abs = acopy[[i, k]].abs();
            if el_abs > pmax {
                p = i;
                pmax = el_abs;
            }
        }
        // 最大ピボットが小さければsingularとしてErrを返す
        if pmax < EPSILON {
            println!("{}, {}", k, pmax);
            return Err(MyLinalgError::SingularMatrix);
        }
        if p != k {
            _gauss_swap(&mut acopy.view_mut(), &mut bcopy.view_mut(), k, p);
        }

        // 前進消去
        for i in (k + 1)..size {
            let ratio = acopy[[i, k]] / acopy[[k, k]];
            acopy[[i, k]] = 0.;
            for j in (k + 1)..size {
                let val = acopy[[k, j]];
                acopy[[i, j]] -= val * ratio;
            }
            let val = bcopy[[k]];
            bcopy[[i]] -= val * ratio;
        }
        // 後退代入
        for i in (0..size).rev() {
            for j in (i + 1)..size {
                let val = acopy[[i, j]] * bcopy[[j]];
                bcopy[[i]] -= val;
                acopy[[i, j]] = 0.;
            }
            let a_ii = acopy[[i, i]];
            bcopy[[i]] /= a_ii;
            acopy[[i, i]] = 1.;
        }
    }
    Ok(bcopy.to_owned())
}

/// 台形則を用いて円周上の周回積分を行う.
/// funcはcurve上の境界要素の添字を受け取りスカラー値を返す関数.
/// 最初と最後の点の関数値は必ず一致することを用いてアルゴリズムを簡略化している.
pub fn circ_integral_trapez<F1>(func: F1, curve: &CircleCurve) -> f64
where
    F1: Fn(usize) -> f64,
{
    let h = std::f64::consts::TAU * curve.radius / (curve.div_num as f64);

    let mut result = 0.0;
    for i in 0..(curve.div_num as usize) {
        result += func(i)
    }
    result *= h;
    result
}

pub fn exact_u(x: &Point<f64>) -> f64 {
    x[0].powf(3.) - 3. * x[0] * x[1].powf(2.)
}

/// 経路は円に限定しているので勾配は単純に計算する.
#[cfg(test)]
pub fn exact_u_normal_dv(x: &Point<f64>) -> f64 {
    let theta = x[1].atan2(x[0]);
    3. * (x[0].powf(2.) - x[1].powf(2.)) * theta.cos() - 6. * x[0] * x[1] * theta.sin()
}

pub fn fund_gamma(x: &Point<f64>, y: &Point<f64>) -> f64 {
    let diff_vec = x - y;
    let norm2 = |v: Point<f64>| Float::sqrt(v.dot(&v));
    -1. * Float::ln(norm2(diff_vec)) / std::f64::consts::TAU
}

pub fn fund_gamma_normal_dv(x: &Point<f64>, y: &Point<f64>) -> f64 {
    let theta = x[1].atan2(x[0]);
    let normal_vec = array![theta.cos(), theta.sin()];
    (x - y).dot(&normal_vec) / std::f64::consts::TAU / (x - y).dot(&(x - y))
}
