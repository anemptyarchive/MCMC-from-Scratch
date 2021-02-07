
# Chapter 4 メトロポリス法 -------------------------------------------------------

# 利用パッケージ
library(tidyverse)
library(gganimate)


# 4.1 メトロポリス法 -------------------------------------------------------------

### 指数関数 -----

# 指数関数のグラフをプロット
tibble(
  x = seq(-1, 1, 0.01), 
  exp_x = exp(x)
) %>% 
  ggplot(aes(x = x, y = exp_x)) + 
    geom_line() + 
    labs(title = expression(y == exp(x)), 
         y = "y")


### 対数関数 -----

# 対数関数のグラフをプロット
tibble(
  x = seq(-1, 2, 0.01), 
  log_x = log(x)
) %>% 
  ggplot(aes(x = x, y = log_x)) + 
    geom_line() + 
    labs(title = expression(y == log(x)), 
         y = "y")


# 4.2 期待値計算の具体例 ------------------------------------------------------------

### 図4.1 メトロポリス法 -----

# 繰り返し回数を指定
K <- 100000

# 一様乱数の範囲(ステップ幅)を指定
c <- 4

# 初期値を指定
x <- 0

# 更新推移の記録用の受け皿を作成
trace_x <- rep(0, K)
n_accept <- 0

# Main loop
for(k in 1:K) {
  # S(x)を計算
  S_x <- 0.5 * x^2
  
  # xの変化量を生成
  delta_x <- runif(n = 1, min = -c, max = c)
  
  # xの更新候補を計算
  dash_x <- x + delta_x
  
  # S(x')を計算
  S_dash_x <- 0.5 * dash_x^2
  
  # テスト値を生成
  r <- runif(n = 1, min = 0, max = 1)
  
  # xを更新
  if(r < exp(S_x - S_dash_x)) { # メトロポリステスト
    x <- dash_x
    n_accept <- n_accept + 1
  }
  
  # 更新値を記録
  trace_x[k] <- x
}

# 更新値のデータフレームを作成
update_df <- tibble(
  k = 1:K, 
  x = trace_x
)

# 真の分布を計算
true_df <- tibble(
  x = seq(-c, c, 0.01), # 描画範囲
  P_x = exp(-x^2 / 2) / sqrt(2 * pi) # 確率密度
  #P_x = dnorm(x = x, mean = 0, sd = 1) # 確率密度
)

# ヒストグラムを作成
ggplot() + 
  geom_histogram(data = update_df, aes(x = x, y = ..density..), 
                 binwidth = 0.05) + # xのヒストグラム
  geom_line(data = true_df, aes(x = x, y = P_x), 
            linetype = "dashed", color = "red") + # 真の確率密度
  labs(title = "Gaussian Distribution", 
       subtitle = paste0("K=", K, ", mu=", 0, ", sigma=", 1))


# 更新値の推移をグラフ化
update_df %>% 
  filter(k >= 1, k <= 2000) %>% # 表示範囲を指定
  ggplot(aes(x = k, y = x)) + 
  geom_line() + 
  labs(title = "Metropolis Method")


# 期待値を計算
expected_df <- tibble(
  k = 1:K, 
  x = trace_x, 
  E_x = cumsum(x) / k, 
  E_x2 = cumsum(x^2) / k
)

# 期待値の推移をプロット
ggplot(expected_df, aes(x = k)) + 
  geom_line(aes(y = E_x), color = "purple") + # xの期待値
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "purple") + # 真のxの期待値
  geom_line(aes(y = E_x2), color = "blue") + # x^2の期待値
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "blue") + # 真のx^2の期待値
  labs(title = "Metropolis Method", 
       y = "E[f(x)]")
  

### 図4.1 メトロポリス法:gif -----

# 繰り返し回数を指定
K <- 10000

# 一様乱数の範囲(ステップ幅)を指定
c <- 4

# xの初期値を指定
x <- 0

# (重いので)表示する間隔を指定
n_interval <- 25

# フレーム数を計算:(割り切れるように要設定)
n_frame <- K / n_interval
n_frame

# Main loop
trace_x <- rep(0, n_frame)
anime_df <- tibble()
for(k in 1:K) {
  # xの更新候補を計算
  dash_x <- x + runif(n = 1, min = -c, max = c)
  
  # テスト値を生成
  r <- runif(n = 1, min = 0, max = 1)
  
  # メトロポリステストによりxを更新
  if(r < exp((0.5 * x^2) - (0.5 * dash_x^2))) {
    x <- dash_x
  }
  
  # 更新値を記録
  trace_x[k] <- x
  
  # 結果を記録
  if(k %% n_interval == 0) {
    
    # k番目までの更新値を格納
    tmp_df <- tibble(
      k = k, 
      x = trace_x
    )
    
    # 更新結果を結合
    anime_df <- rbind(anime_df, tmp_df)
    
    # 途中経過を表示
    print(paste0(k, " (", round(k / K * 100, 1), "%)"))
  }
}

# アニメーション用のヒストグラムを作図
hist_anime <- ggplot() + 
  geom_histogram(data = anime_df, aes(x = x, y = ..density..), 
                 binwidth = 0.05, fill = "purple", color = "purple") + # ヒストグラム
  geom_line(data = true_df, aes(x = x, y = P_x), 
            linetype = "dashed", color = "red") + # 真の確率密度
  gganimate::transition_manual(k) + # フレーム
  ylim(c(0, 1)) + # y軸の範囲
  labs(title = "Metropolis Method", 
       subtitle = "k={current_frame}", 
       x = "x", y = "P(x))")

# gif画像を作成
gganimate::animate(hist_anime, nframes = n_frame, fps = 20)


# 4.3 自己相関 ----------------------------------------------------------------

# 4.3.1 初期値との相関とシミュレーションの熱化 ----------------------------------------------

### 図4.4 初期値との相関 -----

# 確率分布を指定
fn_P <- function(x) {
  exp(-x^2 / 2) / sqrt(2 * pi)
  #dnorm(x = x, mean = 0, sd = 1)
}

# 一様分布の範囲(ステップ幅)を指定
c <- 4

# 真の分布を計算
true_df <- tibble(
  x = seq(-c, c, 0.01), 
  P_x = fn_P(x)
)

# 真の分布をプロット
ggplot(true_df, aes(x = x, y = P_x)) + 
  geom_line()


# 作用(負の対数尤度)を指定
fn_S <- function(x) {
  0.5 * x^2
}

# 確率分布と負の対数尤度をプロット
tibble(
  x = seq(-c, c, 0.01), 
  P_x = fn_P(x), 
  S_x = fn_S(x)
) %>% 
  ggplot(aes(x)) + 
    geom_line(aes(y = P_x)) + 
    geom_line(aes(y = S_x), linetype = "dashed") + 
    labs(y = "f(x)")


# 繰り返し回数を指定
K <- 10000

# 初期値を指定
x <- 20

# Main loop
trace_x <- rep(0, K)
for(k in 1:K) {
  # xの更新候補を計算
  dash_x <- x + runif(n = 1, min = -c, max = c)
  
  # テスト値を生成
  r <- runif(n = 1, min = 0, max = 1)
  
  # メトロポリステストによりxを更新
  if(r < exp(fn_S(x) - fn_S(dash_x))) {
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
update_df %>% 
  filter(k >= 1, k <= 2000) %>% # 表示範囲を指定
  ggplot(aes(x = k, y = x)) + 
  geom_line() + 
  labs(title = "Metropolis Method")


# 期待値を計算
expected_df <- tibble(
  k = 1:K, 
  x = trace_x, 
  E_x = cumsum(x) / k, 
  E_x2 = cumsum(x^2) / k
)

# 期待値の推移をプロット
ggplot(expected_df, aes(x = k)) + 
  geom_line(aes(y = E_x), color = "purple") + # xの期待値
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "purple") + # 真のxの期待値
  geom_line(aes(y = E_x2), color = "blue") + # x^2の期待値
  geom_hline(aes(yintercept = 1), linetype = "dashed", color = "blue") + # 真のx^2の期待値
  ylim(c(-5, 5)) + # y軸の表示範囲
  labs(title = "Metropolis Method", 
       y = "E[f(x)]")


# 4.3.3 自己相関長とジャックナイフ法 ----------------------------------------------------

### ジャックナイフ誤差？ -----

sample_x <- expected_df[["E_x2"]]
# 1グループの要素数を指定
w <- 50

# グループ数を計算
n <- K / w

# 
f_tilde <- rep(0, n)
for(l in 1:n) {
  # l番目のグループの要素を取得
  tmp_x <- sample_x[((l - 1) * w + 1):(l * w)]
  
  # l番目のグループの平均値を計算
  f_tilde[l] <- sum(tmp_x^2) / w
}

f_bar <- sum(sample_x^2) / K
delta_w <- sqrt(sum((f_tilde - f_bar)^2) / n / (n - 1))

jk_w_df <- tibble(
  l = 1:n, 
  x = f_tilde
)

ggplot(jk_w_df, aes(x = l, y = x)) + 
  geom_line()


max_w <- 1000
trace_delta_w <- rep(0, max_w)
for(w in 1:max_w) {
  
  # グループ数を計算
  n <- K / w
  
  # 
  f_tilde <- rep(0, n)
  for(l in 1:n) {
    # l番目のグループの要素を取得
    tmp_x <- sample_x[((l - 1) * w + 1):(l * w)]
    
    # l番目のグループの平均値を計算
    f_tilde[l] <- sum(tmp_x^2) / w
  }
  
  f_bar <- sum(sample_x^2) / K
  trace_delta_w[w] <- sqrt(sum((f_tilde - f_bar)^2) / n / (n - 1))
}

jk_df <- tibble(
  w = 1:max_w, 
  delta_w = trace_delta_w
)

ggplot(jk_df, aes(x = w)) + 
  geom_line(aes(y = delta_w)) + 
  geom_line(aes(y = -delta_w))


# 4.4 ガウス分布以外の例 -----------------------------------------------------------

### 図4.9 混合ガウス分布 -----

# 確率分布を指定
fn_P <- function(x) {
  (exp(-0.5 * (x - 3)^2) + exp(-0.5 * (x + 3)^2)) * 0.5 / sqrt(2 * pi)
}

# 一様分布の範囲(ステップ幅)を指定
c <- 7

# 真の分布を計算
true_df <- tibble(
  x = seq(-c, c, 0.01), 
  P_x = fn_P(x)
)

# 真の分布をプロット
ggplot(true_df, aes(x = x, y = P_x)) + 
  geom_line()


# 作用(負の対数尤度)を指定
fn_S <- function(x) {
  -log(exp(-0.5 * (x - 3)^2) + exp(-0.5 * (x + 3)^2))
}

# 確率分布と負の対数尤度をプロット
tibble(
  x = seq(-c, c, 0.01), 
  P_x = fn_P(x), 
  S_x = fn_S(x)
) %>% 
  ggplot(aes(x)) + 
  geom_line(aes(y = P_x)) + 
  geom_line(aes(y = S_x), linetype = "dashed") + 
  labs(y = "f(x)")


# 繰り返し回数を指定
K <- 100000

# 初期値を指定
x <- 4

# Main loop
trace_x <- rep(0, K)
for(k in 1:K) {
  # xの更新候補を計算
  dash_x <- x + runif(n = 1, min = -c, max = c)
  
  # テスト値を生成
  r <- runif(n = 1, min = 0, max = 1)
  
  # メトロポリステストによりxを更新
  if(r < exp(fn_S(x) - fn_S(dash_x))) {
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


### 図4.10 混合分布 -----

# 確率分布を指定
fn_P <- function(x) {
  # データ数を取得
  n <- length(x)
  
  # 受け皿を初期化
  P_x <- rep(0, n)
  
  # 確率密度を計算
  for(i in 1:n) {
    x_i <- x[i]
    if(x_i >= 0) {
      P_x_i <- exp(-0.5 * x_i^2) / sqrt(2 * pi)#dnorm(x = x_i, mean = 0, sd = 1)
    } else if(x_i >= -1) {
      P_x_i <- 2 / pi * sqrt(1 - x_i^2)
    } else if(x_i < -1) {
      P_x_i <- 0
    }
    P_x[i] <- P_x_i
  }
  
  # 出力
  P_x
}

# 一様分布の範囲(ステップ幅)を指定
c <- 3

# 真の分布を計算
true_df <- tibble(
  x = seq(-c, c, 0.01), 
  P_x = fn_P(x)
)

# 真の分布をプロット
ggplot(true_df, aes(x = x, y = P_x)) + 
  geom_line()


# 作用(負の対数尤度)を指定
fn_S <- function(x) {
  # データ数を取得
  n <- length(x)
  
  # 受け皿を初期化
  S_x <- rep(0, n)
  
  # 確率密度を計算
  for(i in 1:n) {
    x_i <- x[i]
    if(x_i >= 0) {
      S_x_i <- 0.5 * x_i^2 + log(sqrt(2 * pi))
    } else if(x_i >= -1) {
      S_x_i <- -log(2 / pi * sqrt(1 - x_i^2))
    } else if(x_i < -1) {
      S_x_i <- -log(0)
    }
    S_x[i] <- S_x_i
  }
  
  # 出力
  S_x
}

# 確率分布と負の対数尤度をプロット
tibble(
  x = seq(-c, c, 0.01), 
  P_x = fn_P(x), 
  S_x = fn_S(x)
) %>% 
  ggplot(aes(x)) + 
    geom_line(aes(y = P_x)) + 
    geom_line(aes(y = S_x), linetype = "dashed") + 
    labs(y = "f(x)")


# 繰り返し回数を指定
K <- 100000

# 初期値を指定
x <- 0

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


# メトロポリス法:gif -----

# 繰り返し回数を指定
K <- 10000

# 一様分布の範囲(ステップ幅)を指定
c <- 3

# 初期値を指定
x <- 1

# アニメーション表示間隔を指定
n_interval <- 40

# フレーム数を計算:(割り切れるように要設定)
n_frame <- K / n_interval
n_frame

# Main loop
trace_x <- rep(0, K)
anime_df <- tibble()
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
  
  if(k %% n_interval == 0) {
    # k番目までの更新値を格納
    tmp_df <- tibble(
      k = k, 
      x = trace_x
    )
    
    # 結合
    anime_df <- rbind(anime_df, tmp_df)
    
    # 途中経過を表示
    print(paste0(k, " (", round(k / K * 100, 1), "%)"))
  }
}

# アニメーション用のヒストグラムを作図
anime_graph <- ggplot() + 
  geom_histogram(data = anime_df, aes(x = x, y = ..density..), 
                 binwidth = 0.05, fill = "purple", color = "purple") + # xのヒストグラム
  geom_line(data = true_df, aes(x = x, y = P_x), 
            linetype = "dashed", color = "red") + # 真の確率密度
  ylim(c(0, 1)) + # y軸の表示範囲
  gganimate::transition_manual(k) + # フレーム
  labs(title = "Metropolis Method", 
       subtitle = "k={current_frame}")

# gif画像を作成
gganimate::animate(anime_graph, nframes = n_frame, fps = 20)


# 4.5 複雑な数値積分への応用 ---------------------------------------------------------

### 規格化因子の積分 -----


# 作用(負の対数尤度)を指定
fn_S <- function(x) {
  # データ数を取得
  n <- length(x)
  
  # 受け皿を初期化
  S_x <- rep(0, n)
  
  # 確率密度を計算
  for(i in 1:n) {
    x_i <- x[i]
    if(x_i >= 0) {
      S_x_i <- 0.5 * x_i^2 + log(sqrt(2 * pi))
    } else if(x_i >= -1) {
      S_x_i <- -log(2 / pi * sqrt(1 - x_i^2))
    } else if(x_i < -1) {
      S_x_i <- -log(0)
    }
    S_x[i] <- S_x_i
  }
  
  # 出力
  S_x
}

