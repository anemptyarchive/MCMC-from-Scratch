
# Chapter 3 マルコフ連鎖モンテカルロ法の一般論 ---------------------------------------------

# 3章で利用するパッケージ
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
    k = 1:K, # 繰り返し番号
    x = sample(x = c(1, -1), size = K, replace = TRUE, prob = c(p, 1 - p)), # ステップ幅を生成
    y = cumsum(x), # k番目までの合計
    simulation = as.factor(s) # シミュレーション番号
  )
  
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
  geom_line(data = addline_df, aes(x = k, y = root_k), linetype = "dashed") + # 補助線(上)
  geom_line(data = addline_df, aes(x = k, y = -root_k), linetype = "dashed") + # 補助線(下)
  theme(legend.position = "none") + # 凡例
  labs(title = "Random Walk", 
       subtitle = paste0("p=", p, ", K=", K, ", simulation:", n_simu))


# 3.2 規約性 -----------------------------------------------------------------

### ？

# 確率分布を指定
fn_P <- function(x) {
  exp(-1 / x^2 - x^2)
}

# 一様分布の範囲(ステップ幅)を指定
c <- 2

# 真の分布を計算
true_df <- tibble(
  x = seq(-c, c, 0.01), 
  P_x = fn_P(x)
)

# 真の分布をプロット
ggplot(true_df, aes(x = x, y = P_x)) + 
  geom_line()


# 作用を指定
fn_S <- function(x) {
  -1 / x^2 - x^2
}

# 初期値を指定
x <- -2


# 繰り返し回数を指定
K <- 1000

# Main loop
trace_x <- rep(0, K) # 初期化
for(k in 1:K) {
  # xの更新候補を計算
  dash_x <- x + runif(n = 1, min = -c, max = c)
  
  # テスト値を生成
  r <- runif(n = 1, min = 0, max = 1)
  
  # xを更新
  if(r < exp(fn_S(x) - fn_S(dash_x))) { # メトロポリステスト
    x <- dash_x
  }
  
  # 更新値を記録
  trace_x[k] <- x
}

# 更新値のデータフレームを作成
update_df <- tibble(
  k = 1:K, 
  x = trace_x
)

# ヒストグラムを作成
ggplot() + 
  geom_histogram(data = update_df, aes(x = x, y = ..density..), 
                 binwidth = 0.05) + # xのヒストグラム
  geom_line(data = true_df, aes(x = x, y = P_x), 
            linetype = "dashed", color = "red") + # 真の確率密度
  labs(title = "Metropolis Method", 
       subtitle = paste0("K=", K))

# 更新推移をグラフ化
ggplot(update_df, aes(x = k, y = x)) + 
  geom_line() + 
  labs(title = "Metropolis Method")


# 期待値を計算
expected_df <- tibble(
  k = 1:K, 
  x = trace_x
) %>% 
  mutate(E_x = cumsum(x) / k) %>% 
  mutate(E_x2 = cumsum(x^2) / k)

# 期待値の推移をプロット
ggplot(expected_df, aes(x = k)) + 
  geom_line(aes(y = E_x), color = "purple") + # xの期待値
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "purple") + # 真のxの期待値
  geom_line(aes(y = E_x2), color = "blue") + # x^2の期待値
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "blue") + # 真のx^2の期待値
  labs(title = "Metropolis Method", 
       y = "E[f(x)]")


# 3.3 非周期性 --------------------------------------------------------------

### 連続値をとるランダムウォーク

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
    k = 1:K, # 繰り返し番号
    x = runif(n = K, min = -c, max = c), # ステップ幅を生成
    y = cumsum(x), # k番目までの合計
    simulation = as.factor(s) # シミュレーション番号
  )
  
  # 結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # 途中経過を表示
  print(paste0("simuration: ", s, ", y = ", tmp_df[["y"]][K]))
}

# 補助線のデータフレームを作成
addline_df <- tibble(
  k = 1:K, 
  root_k = c * sqrt(k)
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
       subtitle = paste0("c=", c, ", K=", K, ", simulation:", n_simu))

