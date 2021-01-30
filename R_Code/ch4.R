
# Chapter 4 メトロポリス法 -------------------------------------------------------

# 利用パッケージ
library(tidyverse)
library(gganimate)


# 4.1 メトロポリス法 -------------------------------------------------------------

### 図4.1

# 一様乱数のとり得る範囲を指定
c <- 4

# 繰り返し回数を指定
K <- 100000

# 初期値を指定
x <- 0

# 更新推移の記録用の受け皿を作成
trace_x <- rep(0, K)
n_accept <- 0

# maine loop
for(k in 1:K) {
  # S(x)を計算
  S_x <- x^2 / 2
  
  # xの変化量を生成
  delta_x <- runif(n = 1, min = -c, max = c)
  
  # xの更新候補を計算
  dash_x <- x + delta_x
  
  # S(x')を計算
  S_dash_x <- dash_x^2 / 2
  
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

# 更新推移をグラフ化
ggplot(update_df, aes(x = k, y = x)) + 
    geom_line()


# 真の分布を計算
ture_df <- tibble(
  x = seq(-c, c, 0.01), # 描画範囲
  density = exp(-x^2 / 2) / sqrt(2 * pi) # 確率密度
)

# ヒストグラムを作成
ggplot() + 
  geom_histogram(data = update_df, aes(x = x, y = ..density..), 
                 binwidth = 0.1) + # xのヒストグラム
  geom_line(data = ture_df, aes(x = x, y = density), 
            linetype = "dashed", color = "red") # 真の確率密度


### 図4.1:アニメver

# シミュレーション数を指定
K <- 10000

# maine loop
res_df <- tibble()
for(k in 1:K) {
  # 前回のxを取得
  if(k > 1) {
    old_x <- tmp_df[["update_x"]]
  } else if(k == 1) { # 初回の場合
    old_x <- 0 # 初期値を指定
  }
  
  # xを更新
  tmp_df <- tibble(
    k = k, # 繰り返し番号
    x = old_x, 
    S_x = x^2 / 2, # S(x)を計算
    delta_x = runif(n = 1, min = -c, max = c), # xの変化量を生成
    S_dash_x = (x + delta_x)^2 / 2, # S(x')を計算
    r = runif(n = 1, min = 0, max = 1), # テスト値を生成
    update_x = dplyr::if_else(
      r < exp(S_x - S_dash_x), true = x + delta_x, false = x
    ) # メトロポリステスト
  )
  
  # 更新結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # n回ごとに途中経過を表示
  if(k_idx %% 20 == 0) {
    print(paste0(k_idx, " (", round(k_idx / K * 100, 1), "%)"))
  }
}

# gif画像作成用にデータフレームを変形
n_interval <- 20 # 表示する間隔を指定(重いので)
anime_df <- tibble()
for(k_idx in 1:(K %/% n_interval) * n_interval) {
  # k回目までの結果のデータフレームを作成
  tmp_df <- res_df %>% 
    dplyr::filter(k <= k_idx) %>% # k回目までの行を抽出
    dplyr::select(x = update_x) %>% # 更新値列を抽出して列名を変更
    dplyr::mutate(k = k_idx) # フレーム番号列を追加
  
  # 結合
  anime_df <- rbind(anime_df, tmp_df)
  
  # 途中経過を表示
  print(paste0(k_idx, " (", round(k_idx / K * 100, 1), "%)"))
}

# ヒストグラムを作成
hist_anime <- ggplot() + 
  geom_histogram(data = anime_df, aes(x = x, y = ..density..), binwidth = 0.1) + # ヒストグラム
  geom_line(data = ture_df, aes(x = x, y = density), 
            linetype = "dashed", color = "red") + # 真の確率密度
  gganimate::transition_manual(k) + # フレーム
  ylim(c(0, 1)) + # y軸の範囲
  theme(legend.position = "None") + # 凡例
  labs(subtitle = "k = {current_frame}", 
       x = "x", y = "density")

# gif画像を作成
gganimate::animate(hist_anime, nframes = (K %/% n_interval), fps = 10, width = 400, height = 400)


### 更新値の推移:複数ver

# シミュレーション数を指定
n_simu <- 100

# 繰り返し回数を指定
K <- 1000

# main loop
res_df <- tibble() # 受け皿を初期化
for(k in 1:K) {
  # 前回のxを取得
  if(k > 1) {
    old_x <- tmp_df[["update_x"]]
  } else if(k == 1) { # 初回の場合
    old_x <- rep(0, n_simu) # 初期値を指定
  }
  
  # xを更新
  tmp_df <- tibble(
    x = old_x, # 初期値を指定
    S_x = x^2 / 2, # S(x)を計算
    delta_x = runif(n = n_simu, min = -c, max = c), # xの変化量を生成
    S_dash_x = (x + delta_x)^2 / 2, # S(x')を計算
    r = runif(n = n_simu, min = 0, max = 1), # テスト値を生成
    k = k, # 繰り返し番号
    simulation = as.factor(1:n_simu) # シミュレーション番号
  ) %>% 
    dplyr::mutate(
      update_idx = dplyr::if_else(
        r < exp(S_x - S_dash_x), true = 1, false = 0
      ) # メトロポリステスト
    ) %>% 
    dplyr::mutate(update_x = x + update_idx * delta_x) # xを更新
  
  # 結果を結合
  res_df <- rbind(res_df, tmp_df)
  
  # n回ごとに途中経過を表示
  if(k %% 20 == 0) {
    print(paste0("k = ", k, " (", round(k / K * 100, 1), "%)"))
  }
}

# シミュレーションごと推移のグラフを作成
ggplot(res_df, aes(x = k, y = update_x, color = simulation)) + 
  geom_line(alpha = 0.5) + # 推移
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") + # 期待値
  theme(legend.position = "None") + # 凡例
  labs(subtitle = paste0("simulation: ", n_simu), 
       y = "x")

