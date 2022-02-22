# 境界要素法でラプラス方程式を解く（Rust）

タイトルの通り.
![\Delta u=0](https://latex.codecogs.com/svg.image?\Delta&space;u=0)という問題を
![\bar{u}=x^3-3xy^2](https://latex.codecogs.com/svg.image?\bar{u}=x^3-3xy^2)
という Dirichlet 境界条件のもとで解く. 境界は半径 1, 原点中心の単位円（一応, 円の形に限り境界は任意に指定できる）.

## 実行

1. `git clone git@github.com:littleIkawa/bem_laplace_eq_rust.git ***`（`***`は適当なディレクトリ名）
1. `cd ***`
1. `cargo build`
1. `cargo run`
1. `plot.txt`が生成されるので gnuplot でプロットする.
