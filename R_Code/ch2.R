
# Chapter 2 そもそもモンテカルロ法とは -------------------------------------------------

# 2.1 そもそも乱数とは ------------------------------------------------------------

# 利用パッケージ
library(tidyverse)
library(gganimate)

# 2.1.2 ガウス乱数(正規乱数) -------------------------------------------------------

### 図2.1 ガウス関数

# 作図用のデータフレームを作成
gaussian_function_df <- tibble(
  x = seq(-4, 4, 0.01), 
  f_x = exp(-x^2)
)

# 作図
ggplot(gaussian_function_df, aes(x = x, y = f_x)) + 
  geom_line() + 
  labs(title = expression(f(x) == exp(-x^2)), 
       y = "f(x)", x = "x")


### ガウス分布:(平均0,標準偏差1/√2)

# 作図用のデータフレームを作成
gaussian_distribution_df <- tibble(
  x = seq(-4, 4, 0.01), 
  P_x = exp(-x^2) / pi
)

# 作図
ggplot(gaussian_distribution_df, aes(x = x, y = P_x)) + 
  geom_line() + 
  labs(title = expression(P(x) == exp(-x^2) / sqrt(pi)), 
       y = "density", x = "x")


### 図2.2 ガウス分布:(平均mu,標準偏差sigma)

# パラメータを指定
mu <- 0
sigma <- 2

# 作図用のデータフレームを作成
gaussian_distribution_df <- tibble(
  x = seq(-10, 10, 0.01), 
  C_N = sqrt(2 * pi) * sigma, # 規格化係数
  P_x = exp(-(x - mu)^2 / (2 * sigma^2)) / C_N, 
  #P_x = dnorm(x = x, mean = mu, sd = sigma)
)

# 作図
ggplot(gaussian_distribution_df, aes(x = x, y = P_x)) + 
  geom_line() + 
  labs(title = "Gaussian Distribution", 
       subtitle = paste0("mu=", mu, ", sigma=", sigma),
       y = "density", x = "x")


### 図2.3 中心極限定理

# 誤差の種類を指定
K <- 100

# 作図用のデータフレームを作成
gaussian_distribution_df <- tibble(
  x = seq(-1.5, 1.5, 0.01), 
  density = dnorm(x = x, mean = 0, sd = 1)
)

# 誤差を生成
sum_error <- 0 # 変数を初期化
for(k in 1:K) {
  # 一様乱数を生成
  k_error <- runif(n = 1000000, min = -0.5, max = 0.5)
  
  # 誤差を加算
  sum_error <- sum_error + k_error
}

# 作図用のデータフレームを作成
error_df = tibble(
  error = sum_error
)

# for()無しver
error_df <- (runif(n = 1000000 * K, min = -0.5, max = 0.5)) %>% 
  matrix(ncol = K) %>% 
  rowSums() %>% 
  tibble(error = .)

# ヒストグラムを作成
ggplot(error_df, aes(x = error / sqrt(K))) + 
  geom_histogram(binwidth = 0.1) + 
  labs(subtitle = (paste0("K=", K)), 
       x = expression(x / K^0.5))


### アニメーション

# 誤差の種類を指定
K <- 100

# 作図用のデータフレームを作成
error_animation_df <- tibble() # 変数を初期化
for(k in 1:K) {
  # k種類の誤差を加算
  tmp_df <- (runif(n = 10000 * k, min = -0.5, max = 0.5) / sqrt(k)) %>% 
    matrix(ncol = k) %>% 
    rowSums() %>% 
    tibble(
      error = ., 
      k = as.factor(k)
    )
  
  # 計算結果を結合
  error_animation_df <- rbind(error_animation_df, tmp_df)
}

# 作図
error_animation <- ggplot(error_animation_df, aes(x = error)) + 
  geom_histogram(binwidth = 0.1) + 
  gganimate::transition_manual(k) + # フレーム
  labs(subtitle = "K={current_frame}")

# gif画像を作成
gganimate::animate(error_animation, nframes = K, fps = 5)


# 2.2 一様乱数を用いた積分 ----------------------------------------------------------

# 2.2.1 一様乱数を用いた円周率の計算 ----------------------------------------------------

### 図2.4

# 扇形のラインのデータフレームを作成
circle_df <- tibble(
  x = seq(0, 1, 0.01), 
  y = sqrt(1 - x^2)
)

# 扇形のラインを作図
ggplot(circle_df, aes(x = x, y = y)) + 
  geom_line(size = 1) + 
  coord_fixed() + # 縦横比を固定
  labs(title = expression(x^2 + y^2 == 1))


# 繰り返し回数(K)を指定
K <- 1000

# 変数を初期化
n_in <- 0
sample_x <- 0
sample_y <- 0

# Main loop
for(i in 1:K) {
  # 乱数を生成
  x <- runif(n = 1, min = 0, max = 1)
  y <- runif(n = 1, min = 0, max = 1)
  
  # サンプルを記録
  sample_x <- c(sample_x, x)
  sample_y <- c(sample_y, y)
  
  # 扇形内のサンプルをカウント
  if(x^2 + y^2 < 1) { # 乱数が扇の中のとき
    n_in <- n_in + 1
  }
  
  # 途中経過を表示
  print(paste0("iteration: ", i, ", rate: ", round(n_in / i, 2)))
}

# 乱数をデータフレームに格納
sample_df <- tibble(
  x = sample_x, 
  y = sample_y, 
  label = dplyr::if_else(
    x^2 + y^2 < 1, true = "in", false = "out"
  )
)

# 結果を作図
ggplot() + 
  geom_line(data = circle_df, aes(x = x, y = y), size = 1) + # 扇形のライン
  geom_point(data = sample_df, aes(x = x, y = y, color = label)) + # 乱数の散布図
  coord_fixed() + # 縦横比を固定
  labs(title = expression(x^2 + y^2 < 1), 
       subtitle = paste0("K=", K, ", N_in=", n_in))


### 図2.5

# シミュレーション回数を指定
n_simu <- 100

# シミュレーション
res_df <- tibble() # 変数を初期化
for(s in 1:n_simu) {
  # Main loopの処理
  tmp_df <- tibble(
    x = runif(n = K, min = 0, max = 1), 
    y = runif(n = K, min = 0, max = 1), 
    label = dplyr::if_else(
      x^2 + y^2 < 1, true = TRUE, false = FALSE
    ), # 扇形の中か判定
    k = 1:K, # 繰り返し回数
    simulation = as.factor(s) # 何回目のシミュレーションか
  ) %>% 
    dplyr::mutate(n_in = cumsum(label) / k) # in率を計算
  
  # 結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # 途中経過を表示
  print(paste0(
    "simulation: ", s, ", rate: ", round(tmp_df[["label"]][K], 2)
  ))
}

# 結果を作図
ggplot(res_df, aes(x = k, y = n_in, color = simulation)) + 
  geom_line(alpha = 0.5) + # in率(面積)の推移
  geom_hline(aes(yintercept = pi / 4), linetype = "dashed") + # 真の面積
  theme(legend.position = "none") + # 凡例
  labs(subtitle = paste0("simulation: ", n_simu), 
       y = "")


# 2.2.2 一様乱数を用いたt定積分 ----------------------------------------------------

### 図2.6

# 定数を指定
a <- 0
b <- 1

# シミュレーション回数を指定
n_simu <- 100

# シミュレーション
res_df <- tibble() # 初期化
for(s in 1:n_simu) {
  # Main loopの処理
  tmp_df <- tibble(
    x = runif(n = K, min = a, max = b), # 乱数を生成
    y = sqrt(1 - x^2), 
    k = 1:K, # 繰り返し回数
    simulation = as.factor(s) # 何回目のシミュレーションか
  ) %>% 
    dplyr::mutate(E_f_x = cumsum(y) / k) # 関数の平均を計算
  
  # 結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # 途中経過を表示
  print(paste0(
    "simulation: ", s, ", rate: ", round(tmp_df[["E_f_x"]][K], 2)
  ))
}

# 結果を作図
ggplot(res_df, aes(x = k, y = E_f_x, color = simulation)) + 
  geom_line(alpha = 0.5) + # in率(面積)の推移
  geom_hline(aes(yintercept = pi / 4), linetype = "dashed") + # 真の面積
  theme(legend.position = "none") + # 凡例
  labs(title = expression(sum(sqrt(1 - x[(k)]^2), k==1, K) / k), 
       subtitle = paste0("simulation: ", n_simu), 
       y = "E[f(x)]")


# 2.2.3 ガウス積分 -------------------------------------------------------------

# 積分範囲を指定
a = 10000

# 扇形のラインのデータフレームを作成
gaussian_distribution_df <- tibble(
  x = seq(-a, a, 0.01), 
  P_x = exp(-x^2 / 2) / sqrt(2 * pi) # 確率密度
)

# 扇形のラインを作図
ggplot(gaussian_distribution_df, aes(x = x, y = P_x)) + 
  geom_line() + 
  labs(title = expression(x^2 + y^2 == 1))


### 図2.7？

# シミュレーション回数を指定
n_simu <- 100

# シミュレーション
res_df <- tibble() # 初期化
for(s in 1:n_simu) {
  # Main loopの処理
  tmp_df <- tibble(
    x = runif(n = K, min = -a, max = a), 
    P_x = exp(-(x)^2 / 2) / sqrt(2 * pi), 
    k = 1:K, 
    simulation = as.factor(s)
  ) %>% 
    dplyr::mutate(E_P_x = cumsum(P_x) / k)
  
  # 結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # 途中経過を表示
  print(paste0(
    "simulation: ", s, ", rate: ", round(tmp_df[["E_P_x"]][K], 2)
  ))
}

# 結果を作図
ggplot(res_df, aes(x = k, y = E_P_x, color = simulation)) + 
  geom_line(alpha = 0.5) + 
  theme(legend.position = "none") + # 凡例
  labs(y = "E[P(x)]")


