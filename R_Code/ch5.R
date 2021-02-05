
# Chapter 5 多変数のメトロポリス法 ---------------------------------------------------

# 5章で利用するパッケージ
library(tidyverse)


# 5.1 多変数のガウス分布 ------------------------------------------------------------

### 図5.1 多変数ガウス分布 -----

# 作用(負の対数尤度)を指定
fn_S <- function(x, y) {
  (x^2 + y^2 + x * y) * 0.5
}

# 確率分布を指定
fn_P <- function(x, y) {
  exp(-fn_S(x, y))
}

# 一様分布の範囲(ステップ幅)を指定
c <- 4

# 真の分布を計算
vec <- seq(-c, c, 0.01) # 描画用の点
true_df <- tibble(
  x = rep(vec, times = length(vec)), 
  y = rep(vec, each = length(vec)), 
  P_xy = fn_P(x, y) # 確率密度
)

# 真の分布をプロット
ggplot(true_df, aes(x = x, y = y, z = P_xy, color = ..level..)) + 
  geom_contour() +
  labs(title = "Mixture Gaussian Distribution")


# 繰り返し回数を指定
K <- 10000

# 初期値を指定
x <- 0
y <- 0

# Main loop
trace_x <- rep(0, K)
trace_y <- rep(0, K)
for(k in 1:K) {
  # 更新候補を計算
  dash_x <- x + runif(n = 1, min = -c, max = c)
  dash_y <- y + runif(n = 1, min = -c, max = c)
  
  # テスト値を生成
  r <- runif(n = 1, min = 0, max = 1)
  
  # メトロポリステストにより更新
  if(r < exp(fn_S(x, y) - fn_S(dash_x, dash_y))) {
    x <- dash_x
    y <- dash_y
  }
  
  # 更新値を記録
  trace_x[k] <- x
  trace_y[k] <- y
}

# 更新値のデータフレームを作成
update_df <- tibble(
  k = 1:K, 
  x = trace_x, 
  y = trace_y
)

# 散布図を作成
ggplot() + 
  geom_point(data = update_df, aes(x = x, y = y)) + # サンプルの散布図
  geom_contour(data = true_df, aes(x = x, y = y, z = P_xy, color = ..level..)) + # 真の確率密度
  labs(title = "Metropolis Method", 
       subtitle = paste0("K=", K))


### 図5.3-4 :1つ -----

# 関数を指定
fn_s_physics <- function(x) {
  (2 + tanh(x)) / 3
}
fn_s_baseball <- function(y) {
  #  データ数を取得
  n <- length(y)
  
  # 受け皿を初期化
  P_y <- rep(0, n)
  
  # 出力を計算
  for(i in 1:n) {
    y_i <- y[i]
    if(y_i <= 2) {
      P_y_i <- 0
    } else if(y_i > 2) {
      P_y_i <- y_i^2 * 0.5
    }
    P_y[i] <- P_y_i
  }
  
  # 出力
  P_y
}

# 期待値を計算
expected_df <- update_df %>% 
  dplyr::mutate(
    s_x = fn_s_physics(x), 
    s_y = fn_s_baseball(y)
  ) %>% 
  dplyr::mutate(
    E_s_x = cumsum(s_x) / k, 
    E_s_y = cumsum(s_y) / k
  )

# 推移をプロット
ggplot(expected_df, aes(x = k)) + 
  geom_line(aes(y = E_s_x), color = "orange") + 
  geom_hline(aes(yintercept = 2 / 3), linetype = "dashed", color = "orange") + 
  geom_line(aes(y = E_s_y), color = "#00A986") + 
  geom_hline(aes(yintercept = 0.1305), linetype = "dashed", color = "#00A986") + 
  labs(y = "E[f]")


### 図5.3-4 :複数 -----

# 繰り返し回数を指定
K <- 10000

# シミュレーション回数を指定
n_simu <- 100

# Main loop
trace_df <- tibble()
for(k in 1:K) {
  # 前ステップの変数を取得
  if(k > 1) {
    x <- tmp_df[["new_x"]]
    y <- tmp_df[["new_y"]]
  } else if(k == 1) { # 初回のみ
    # 初期値を指定
    x <- rep(0, n_simu)
    y <- rep(0, n_simu)
  }
  
  # 全てのシミュレーションにおけるk回目の更新
  tmp_df <- tibble(
    simulation = as.factor(1:n_simu), 
    k = k, 
    old_x = x, 
    old_y = y, 
    S_xy = fn_S(old_x, old_y), 
    delta_x = runif(n = n_simu, min = -c, max = c), 
    delta_y = runif(n = n_simu, min = -c, max = c), 
    S_dash_xy = fn_S(old_x + delta_x, old_y + delta_y), 
    r = runif(n = n_simu, min = 0, max = 1)
  ) %>% 
    # メトロポリステスト
    dplyr::mutate(
      test = dplyr::if_else(
        r < exp(S_xy - S_dash_xy), true = 1, false = 0
      )
    ) %>% 
    # 値を更新
    dplyr::mutate(
      new_x = old_x + test * delta_x, 
      new_y = old_y + test * delta_y
    ) 
  
  # 結果を結合
  trace_df <- rbind(trace_df, tmp_df)
  
  # n回ごとに途中経過を表示
  if(k %% 25 == 0) {
    print(paste0(
      "k=", k, " (", round(k / K * 100, 1), "%)"
    ))
  }
}

# 期待値を計算
expected_df <- tibble()
for(s in 1:n_simu) {
  # s番目のシミュレーションを計算
  tmp_df <- trace_df %>% 
    dplyr::filter(simulation == s) %>% 
    dplyr::select(
      simulation, k, x = new_x, y = new_y
    ) %>% 
    # 関数の計算
    dplyr::mutate(
      s_x = fn_s_physics(x), 
      s_y = fn_s_baseball(y)
    ) %>% 
    # 期待値を計算
    dplyr::mutate(
      E_s_x = cumsum(s_x) / k, 
      E_s_y = cumsum(s_y) / k
    )
  
  # 計算結果を結合
  expected_df <- rbind(expected_df, tmp_df)
}

# 図5.3をプロット
ggplot(expected_df, aes(x = k, y = E_s_x, color = simulation)) + 
  geom_line(alpha = 0.5) + 
  geom_hline(aes(yintercept = 2 / 3), linetype = "dashed", color = "orange") + 
  theme(legend.position = "none") + 
  labs(title = "Metropolis Method", 
       subtitle = paste0("simulation:", n_simu), 
       y = "E[s_phy(x)]")

# 図5.4をプロット
ggplot(expected_df, aes(x = k, y = E_s_y, color = simulation)) + 
  geom_line(alpha = 0.5) + 
  geom_hline(aes(yintercept = 0.1305), linetype = "dashed", color = "#00A986") + 
  theme(legend.position = "none") + 
  labs(title = "Metropolis Method", 
       subtitle = paste0("simulation:", n_simu), 
       y = "E[s_bb(y)]")

