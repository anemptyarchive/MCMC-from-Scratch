
# Chapter 3 マルコフ連鎖モンテカルロ法の一般論 ---------------------------------------------

# 利用パッケージ
library(tidyverse)


# 3.1 マルコフ連鎖 --------------------------------------------------------------

### 表3.1 ランダムウォーク -----

# 繰り返し回数を指定
K <- 1000

# 確率を指定
p <- 0.5

# シミュレーション回数を指定
n_simu <- 100

# シミュレーション
res_df <- tibble() # 初期化
for(s in 1:n_simu) {
  # ランダムウォーク
  tmp_df <- tibble(
    x = sample(x = c(1, -1), size = K, replace = TRUE, prob = c(p, 1 - p)), 
    k = 1:K, # 繰り返し番号
    simulation = as.factor(s) # シミュレーション番号
  ) %>% 
    dplyr::mutate(y = cumsum(x)) # k番目までの合計
  
  # 結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # 途中経過を表示
  print(paste0("simuration: ", s, ", y = ", tmp_df[["y"]][K]))
}

# 補助線のデータフレームを作成
addline_df <- tibble(
  k = 1:K, 
  root_k = 2 * sqrt(k)
)

# 出力の推移をプロット
ggplot() + 
  geom_line(data = res_df, aes(x = k, y = y, color = simulation), 
            alpha = 0.5) + # yの推移
  geom_hline(aes(yintercept = 0), linetype = "dashed") + # 期待値
  geom_line(data = addline_df, aes(x = k, y = root_k), linetype = "dashed") + 
  geom_line(data = addline_df, aes(x = k, y = -root_k), linetype = "dashed") + 
  theme(legend.position = "none") + # 凡例
  labs(title = "Random Walk", 
       subtitle = paste0("p=", p, ", K=", k, ", simulation:", n_simu))


# 3.2 規約性 -----------------------------------------------------------------

### ？

# 定数を指定
c <- 1

# xの分布を確認
f_df <- tibble(
  x = seq(-c, c, 0.01), 
  y = exp(-1 / x^2 - x^2)
)
ggplot(f_df, aes(x = x, y = y)) + 
  geom_line()


# 繰り返し回数を指定
K <- 1000

# シミュレーション回数を指定
n_simu <- 10

# 描画範囲を生成
x <- seq(-c, c, 0.01)

# 確率分布を計算
P_x <- exp(-1 / x^2 - x^2)

# シミュレーション
res_df <- tibble() # 初期化
for(s in 1:n_simu) {
  # サンプルを生成
  tmp_df <- tibble(
    delta_x = sample(x = x, size = n_iter, replace = TRUE, prob = P_x), 
    k = 1:K, 
    simulation = as.factor(s)
  ) %>% 
    dplyr::mutate(x = cumsum(delta_x) / k) # k番目までの合計
  
  # 結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # 途中経過を表示
  print(paste0("simuration: ", s, ", y = ", tmp_df[["x"]][K]))
}

# 作図
ggplot(res_df, aes(x = k, y = x, color = simulation)) + 
  geom_line(alpha = 0.5) + # yの推移
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") + # 期待値
  theme(legend.position = "none") + # 凡例
  labs(title = expression(x^(k+1) == x^(k) + exp(-1 / x^2 - x^2)), 
       subtitle = paste0("simulation: ", n_simu), 
       x = "k")

# 3.3 非周期性 --------------------------------------------------------------

# 繰り返し回数を指定
K <- 1000

# 一様分布の幅を指定
c <- 3

# シミュレーション回数を指定
n_simu <- 100

# シミュレーション
res_df <- tibble() # 初期化
for(s in 1:n_simu) {
  # ランダムウォーク
  tmp_df <- tibble(
    x = runif(n = K, min = -c, max = c), 
    k = 1:K, # 繰り返し番号
    simulation = as.factor(s) # シミュレーション番号
  ) %>% 
    dplyr::mutate(y = cumsum(x)) # k番目までの合計
  
  # 結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # 途中経過を表示
  print(paste0("simuration: ", s, ", y = ", tmp_df[["y"]][K]))
}

# 出力の推移をプロット
ggplot(res_df, aes(x = k, y = y, color = simulation)) + 
  geom_line(alpha = 0.5) + # yの推移
  geom_hline(aes(yintercept = 0), linetype = "dashed") + # 期待値
  theme(legend.position = "none") + # 凡例
  labs(title = "Random Walk", 
       subtitle = paste0("p=", p, ", K=", k, "simulation:", n_simu))

