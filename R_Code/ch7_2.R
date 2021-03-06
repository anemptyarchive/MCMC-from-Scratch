
# Chapter 7 マルコフ連鎖モンテカルロ法の応用例 ---------------------------------------------

# 7.2 イジング模型 ---------------------------------------------------------------

# 7.2.1 メトロポリス法 -----------------------------------------------------------

# 7.2.1項で利用するパッケージ
library(tidyverse)
library(gganimate)


### 関数定義 -----

# エネルギーの計算関数を作成:式(7.45)
fn_E <- function(spin_mat, J = 1, h = 0) {
  # 1辺のスピン数を取得
  N <- nrow(spin_mat)
  
  # 変数を初期化
  sum_s_i <- 0  # 全てのスピンの和
  sum_s_ij <- 0 # 隣接する全ての組み合わせのスピンの積の和
  
  # 列インデックスを順番に生成
  for(i_x in 1:N) {
    # 右隣のインデックスを取得
    i_x_plus <- i_x + 1
    
    # 行インデックスを順番に生成
    for(i_y in 1:N) {
      # 下隣のインデックスを取得
      i_y_plus <- i_y + 1
      
      # i番目のスピンを加算
      sum_s_i <- sum_s_i + spin_mat[i_y, i_x]
      
      # i番目と右隣のスピンとの積を加算
      if(i_x < N) { # N列目のとき右隣はない
        sum_s_ij <- sum_s_ij + spin_mat[i_y, i_x] * spin_mat[i_y, i_x_plus]
      }
      
      # i番目と下隣のスピンとの積を加算
      if(i_y < N) { # N行目のとき下隣はない
        sum_s_ij <- sum_s_ij + spin_mat[i_y, i_x] * spin_mat[i_y_plus, i_x]
      }
    }
  }
  
  # エネルギーを計算:式(7.45)
  energy <- - J * sum_s_ij - h * sum_s_i
  return(energy)
}

# エネルギーの計算関数を作成:式(7.49)
fn_simpleE <- function(spin_mat, J = 1, h = 0) {
  # 1辺のスピン数を取得
  N <- nrow(spin_mat)
  
  # 隣接する全ての組み合わせのスピンの積の和を初期化
  sum_s_ij <- 0
  
  # i番目と右隣のスピンとの積を加算
  if(i_x < N) { # N列目のとき右隣はない
    sum_s_ij <- sum_s_ij + spin_mat[i_y, i_x] * spin_mat[i_y, i_x + 1]
  }
  
  # i番目と下隣のスピンとの積を加算
  if(i_y < N) { # N行目のとき下隣はない
    sum_s_ij <- sum_s_ij + spin_mat[i_y, i_x] * spin_mat[i_y + 1, i_x]
  }

  # i番目と左隣のスピンとの積を加算
  if(i_x > 1) { # N列目のとき左隣はない
    sum_s_ij <- sum_s_ij + spin_mat[i_y, i_x] * spin_mat[i_y, i_x - 1]
  }
  
  # i番目と上隣のスピンとの積を加算
  if(i_y > 1) { # N行目のとき上隣はない
    sum_s_ij <- sum_s_ij + spin_mat[i_y, i_x] * spin_mat[i_y - 1, i_x]
  }
  
  # エネルギーを計算:式(7.45)
  energy <- - J * sum_s_ij - h * spin_mat[i_y, i_x]
  return(energy)
}


### 初期値の設定 -----

# 1辺のスピン数を指定
N <- 50

# スピンのマトリクスを初期化
spin_mat <- matrix(sample(x = c(-1, 1), size = N^2, replace = TRUE), nrow = N, ncol = N)

# 作図用にスピンのデータフレームを作成
spin_df <- tibble(
  i_y = rep(1:N, times = N), 
  i_x = rep(1:N, each = N), 
  spin = as.vector(spin_mat), 
  label = as.factor(paste0("iter:", 0, ", rate:", sum(spin_mat) / N^2))
)

# 初期値を作図
ggplot(spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  labs(subtitle = paste0("rate:", sum(spin_mat) / N^2))


### メトロポリス法によるシミュレーション -----

# パラメータを指定
J <- 1
h <- -0.5
temperature <- 1.5

# 最大試行回数を指定
max_iter <- 50

# 試行回数を初期化
iter <- 0

# メトロポリス法
while(abs(sum(spin_mat)) / N^2 < 0.9) { # 指定したレートに達するまで
  # 試行回数を加算
  iter <- iter + 1

  # 更新するスピン番号を生成
  i_vec <- 1:N^2 # 順番に選択
  #i_vec <- sample(1:N^2, size = N^2, replace = FALSE) # 順番をランダムに選択
  #i_vec <- sample(1:N^2, size = N^2, replace = TRUE) # ランダムに選択
  
  # 1試行における更新回数を初期化
  n_accept <- 0
  
  # 配位を更新
  for(i in i_vec) {
    # スピンのインデックスを計算
    i_y <- ifelse(i %% N == 0, yes = N, no = i %% N) # 行番号
    i_x <- (i - 1) %/% N + 1 # 列番号
    
    # エネルギーを計算:式(7.49)
    energy <- fn_simpleE(spin_mat, J, h)
    
    # i番目のスピンを反転させたマトリクスを作成
    spin_dash_mat <- spin_mat
    spin_dash_mat[i_y, i_x] <- spin_dash_mat[i_y, i_x] * (-1)
    
    # エネルギーを計算:式(7.49)
    energy_dash <- fn_simpleE(spin_dash_mat, J, h)
    
    # 判定用の確率値を計算:式(4.46')
    prob <- exp((energy - energy_dash) / temperature)
    
    # テスト値を生成
    metropolis <- runif(n = 1, min = 0, max = 1)
    
    # メトロポリステストにより配位を更新
    if(prob > metropolis) {
      # スピンを更新
      spin_mat <- spin_dash_mat
      
      # 更新回数を加算
      n_accept <- n_accept + 1
    }
    
  }
  
  # スピンのデータフレームを作成
  tmp_spin_df <- tibble(
    i_y = rep(1:N, times = N), 
    i_x = rep(1:N, each = N), 
    spin = as.vector(spin_mat), 
    label = as.factor(paste0(
      "J=", J, ", h=", h, ", T=", temperature, 
      ", iter:", iter, ", rate:", sum(spin_mat) / N^2
    ))
  )
  
  # 結果を結合
  spin_df <- rbind(spin_df, tmp_spin_df)
  
  # 途中経過を表示
  print(paste0("iter:", iter, ", accept:", n_accept, ", rate:", sum(spin_mat) / N^2))
  
  # 最大回数に達したら終了
  if(iter == max_iter) break
}

# 最終結果を作図
ggplot(tmp_spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  labs(title = "Metropolis Method", 
       subtitle = paste0("J=", J, ", h=", h, ", T=", temperature, 
                         ", iter:", iter, ", rate:", sum(spin_mat) / N^2))


### アニメーションの作成 -----

# アニメーション用のグラフを作成
graph <- ggplot(spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  gganimate::transition_manual(label) + # フレーム
  labs(title = "Metropolis Method", 
       subtitle = "{current_frame}")

# gif画像を出力
gganimate::animate(graph, nframes = iter + 1, fps = 10)



# 7.2.2 ギブスサンプリング ---------------------------------------------------------

# 7.2.2項で利用するパッケージ
library(tidyverse)
library(gganimate)


### 関数定義 -----

# エネルギーの計算関数を作成:式(7.45)
fn_E <- function(spin_mat, J = 1, h = 0) {
  # 1辺のスピン数を取得
  N <- nrow(spin_mat)
  
  # 全てのスピンの和
  sum_s_i <- sum(spin_mat)
  
  # 右隣のスピンとの積の和
  sum_s_ij_x <- sum(spin_mat[, -N] * spin_mat[, -1])
  
  # 下隣のスピンとの積の和
  sum_s_ij_y <- sum(spin_mat[-1, ] * spin_mat[-N, ])
  
  # エネルギーを計算:式(7.45)
  energy <- - J * (sum_s_ij_x + sum_s_ij_y) - h * sum_s_i
  return(energy)
}


### 初期値の設定 -----

# 1辺のスピン数を指定
N <- 50

# スピンのマトリクスを初期化
spin_mat <- matrix(sample(x = c(-1, 1), size = N^2, replace = TRUE), nrow = N, ncol = N)

# 作図用にスピンのデータフレームを作成
spin_df <- tibble(
  i_y = rep(1:N, times = N), 
  i_x = rep(1:N, each = N), 
  spin = as.vector(spin_mat), 
  label = as.factor(paste0("iter:", 0, ", rate:", sum(spin_mat) / N^2))
)

# 初期値を作図
ggplot(spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  labs(subtitle = paste0("rate:", sum(spin_mat) / N^2))


### ギブスサンプリングによるシミュレーション -----

# パラメータを指定
J <- 1
h <- 0
temperature <- 2

# 最大試行回数を指定
max_iter <- 50

# 試行回数を初期化
iter <- 0

# ギブスサンプリング
while(abs(sum(spin_mat)) / N^2 < 0.9) { # 指定したレートに達するまで
  # 試行回数を加算
  iter <- iter + 1
  
  # 更新するスピン番号を生成
  i_vec <- 1:N^2 # 順番に選択
  #i_vec <- sample(1:N^2, size = N^2, replace = FALSE) # 順番をランダムに選択
  #i_vec <- sample(1:N^2, size = N^2, replace = TRUE) # ランダムに選択
  
  # 配位を更新
  for(i in i_vec) {
    # スピンのインデックスを計算
    i_y <- ifelse(i %% N == 0, yes = N, no = i %% N) # 行番号
    i_x <- (i - 1) %/% N + 1 # 列番号
    
    # i番目のスピンを反転させたマトリクスを作成
    spin_plus_mat <- spin_mat
    spin_plus_mat[i_y, i_x] <- 1
    spin_minus_mat <- spin_mat
    spin_minus_mat[i_y, i_x] <- -1
    
    # エネルギーを計算:式(7.45)と式(7.46)の分子の指数部分
    energy_plus  <- - fn_E(spin_plus_mat, J, h) / temperature
    energy_minus <- - fn_E(spin_minus_mat, J, h) / temperature
    
    # エネルギーの指数を計算:式(7.46)の分子
    if(max(energy_plus, energy_minus) > 700) {
      # 最大値を取得
      max_energy <- max(energy_plus, energy_minus)
      
      # オーバーフロー対策
      exp_energy_plus  <- exp(energy_plus - max_energy)
      exp_energy_minus <- exp(energy_minus - max_energy)
    } else if(min(energy_plus, energy_minus) < - 700) {
      # 最小値を取得
      min_energy <- min(energy_plus, energy_minus)
      
      # アンダーフロー対策
      exp_energy_plus  <- exp(energy_plus - min_energy)
      exp_energy_minus <- exp(energy_minus - min_energy)
    } else {
      exp_energy_plus  <- exp(energy_plus)
      exp_energy_minus <- exp(energy_minus)
    }
    
    # 確率を計算:式(7.46)
    p_plus <- exp_energy_plus / (exp_energy_plus + exp_energy_minus)
    p_minus <- exp_energy_minus / (exp_energy_plus + exp_energy_minus)
    
    # 確率に従いスピンを更新
    spin_mat[i_y, i_x] <- sample(c(-1, 1), size = 1, prob = c(p_minus, p_plus))
  }
  
  # スピンのデータフレームを作成
  tmp_spin_df <- tibble(
    i_y = rep(1:N, times = N), 
    i_x = rep(1:N, each = N), 
    spin = as.vector(spin_mat), 
    label = as.factor(paste0(
      "J=", J, ", h=", h, ", T=", temperature, 
      ", iter:", iter, ", rate:", sum(spin_mat) / N^2
    ))
  )
  
  # 結果を結合
  spin_df <- rbind(spin_df, tmp_spin_df)
  
  # 途中経過を表示
  print(paste0("iter:", iter, ", rate:", sum(spin_mat) / N^2))
  
  # 最大回数に達したら終了
  if(iter == max_iter) break
}

# 最終結果を作図
ggplot(tmp_spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  labs(title = "Gibbs Sampling", 
       subtitle = paste0("J=", J, ", h=", h, ", T=", temperature, 
                         ", iter:", iter, ", rate:", sum(spin_mat) / N^2))


### アニメーションの作成 -----

# アニメーション用のグラフを作成
graph <- ggplot(spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  gganimate::transition_manual(label) + # フレーム
  labs(title = "Gibbs Sampling", 
       subtitle = "{current_frame}")

# gif画像を出力
gganimate::animate(graph, nframes = iter + 1, fps = 10)


# 7.2.5 Wolffのアルゴリズム ------------------------------------------------------------

# 7.2.5項で利用するパッケージ
library(tidyverse)
library(gganimate)


### 初期値の設定 -----

# 1辺のスピン数を指定
N <- 50

# スピンのマトリクスを初期化
spin_mat <- matrix(sample(x = c(-1, 1), size = N^2, replace = TRUE), nrow = N, ncol = N)

# 作図用にスピンのデータフレームを作成
spin_df <- tibble(
  i_y = rep(1:N, times = N), 
  i_x = rep(1:N, each = N), 
  spin = as.vector(spin_mat), 
  label = as.factor(paste0(
    "iter:", 0, ", rate:", sum(spin_mat) / N^2
  ))
)

# 初期値を作図
ggplot(spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  labs(subtitle = paste0("rate:", sum(spin_mat) / N^2))


### Wolffのアルゴリズムによるシミュレーション -----

# パラメータを指定
J <- 1
h <- 0
temperature <- 2

# クラスタの追加判定用の確率値を計算
prob <- 1 - exp(-2 * J / temperature)

# 最大試行回数を指定
max_iter <- 50

# 試行回数を初期化
iter <- 0

# Wolffアルゴリズム
while(abs(sum(spin_mat)) / N^2 < 0.9) { # 指定したレートに達するまで
  # 試行回数を加算
  iter <- iter + 1

  # クラスタのマトリクスを初期化
  cluster_mat <- matrix(1, nrow = N, ncol = N)
  
  # クラスタの起点となるスピンのインデックスを生成
  i_y <- sample(x = 1:N, size = 1)
  i_x <- sample(x = 1:N, size = 1)
  
  # 起点のクラスタを設定
  cluster_mat[i_y, i_x] <- 0
  
  # 起点のクラスタのスピンを取得
  spin_cluster <- spin_mat[i_y, i_x]
  
  # クラスタのインデックスを記録
  cluster_idx_mat <- matrix(c(i_y, i_x), nrow = 1, ncol = 2)
  
  # クラスタ数を初期化
  n_cluster <- 1
  
  # 処理済みのクラスタ数を初期化
  k <- 0
  
  # 配位を更新
  while(k < n_cluster) {
    # 処理済みのクラスタ数を加算
    k <- k + 1
    
    # クラスタのインデックスを取得
    i_y <- cluster_idx_mat[k, 1]
    i_x <- cluster_idx_mat[k, 2]
    
    # 隣接するインデックスを計算
    i_y_plus <- i_y + 1
    i_x_plus <- i_x + 1
    i_y_minus <- i_y - 1
    i_x_minus <- i_x - 1
    
    # 上隣を処理
    if(i_y_minus >= 1) { # 枠内である
      if(spin_mat[i_y_minus, i_x] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[i_y_minus, i_x] == 1) { # 未処理である
          # 判定用の確率値を生成
          r <- runif(n = 1, min = 0, max = 1)
          if(prob > r) { # 確率で判定
            # クラスタを追加
            cluster_mat[i_y_minus, i_x] <- 0
            
            # 新たなクラスタのインデックスを記録
            cluster_idx_mat <- rbind(cluster_idx_mat, c(i_y_minus, i_x))
            
            # クラスタ数を加算
            n_cluster <- n_cluster + 1
          }
        }
      }
    }
    
    # 下隣を処理
    if(i_y_plus <= N) { # 枠内である
      if(spin_mat[i_y_plus, i_x] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[i_y_plus, i_x] == 1) { # 未処理である
          # 判定用の確率値を生成
          r <- runif(n = 1, min = 0, max = 1)
          if(prob > r) { # 確率で判定
            # クラスタを追加
            cluster_mat[i_y_plus, i_x] <- 0
            
            # 新たなクラスタのインデックスを記録
            cluster_idx_mat <- rbind(cluster_idx_mat, c(i_y_plus, i_x))
            
            # クラスタ数を加算
            n_cluster <- n_cluster + 1
          }
        }
      }
    }
    
    # 左隣を処理
    if(i_x_minus >= 1) { # 枠内である
      if(spin_mat[i_y, i_x_minus] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[i_y, i_x_minus] == 1) { # 未処理である
          # 判定用の確率値を生成
          r <- runif(n = 1, min = 0, max = 1)
          if(prob > r) { # 確率で判定
            # クラスタを追加
            cluster_mat[i_y, i_x_minus] <- 0
            
            # 新たなクラスタのインデックスを記録
            cluster_idx_mat <- rbind(cluster_idx_mat, c(i_y, i_x_minus))
            
            # クラスタ数を加算
            n_cluster <- n_cluster + 1
          }
        }
      }
    }
    
    # 右隣を処理
    if(i_x_plus <= N) { # 枠内である
      if(spin_mat[i_y, i_x_plus] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[i_y, i_x_plus] == 1) { # 未処理である
          # 判定用の確率値を生成
          r <- runif(n = 1, min = 0, max = 1)
          if(prob > r) { # 確率で判定
            # クラスタを追加
            cluster_mat[i_y, i_x_plus] <- 0
            
            # 新たなクラスタのインデックスを記録
            cluster_idx_mat <- rbind(cluster_idx_mat, c(i_y, i_x_plus))
            
            # クラスタ数を加算
            n_cluster <- n_cluster + 1
          }
        }
      }
    }
    
  }
  
  # クラスタのスピンを反転
  spin_mat[cluster_mat == 0] <- spin_cluster * (-1)
  
  # スピンのデータフレームを作成
  tmp_spin_df <- tibble(
    i_y = rep(1:N, times = N), 
    i_x = rep(1:N, each = N), 
    spin = as.vector(spin_mat), 
    label = as.factor(paste0(
      "J=", J, ", h=", h, ", T=", temperature, 
      ", iter:", iter, ", rate:", sum(spin_mat) / N^2
    ))
  )
  
  # 結果を結合
  spin_df <- rbind(spin_df, tmp_spin_df)
  
  # 結果を表示
  print(paste0("iter:", iter, ", cluster:", k, ", rate:", sum(spin_mat) / N^2))
  
  # 最大回数に達したら終了
  if(iter == max_iter) break
}

# 最終結果を作図
ggplot(tmp_spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  labs(title = "Wolff Algorithm", 
       subtitle = paste0("J=", J, ", h=", h, ", T=", temperature, 
                         ", iter:", iter, ", rate:", sum(spin_mat) / N^2))


### アニメーションの作成 -----

# アニメーション用のグラフを作成
graph <- ggplot(spin_df, aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  gganimate::transition_manual(label) + 
  coord_fixed(ratio = 1) + # アスペクト比
  labs(title = "Wolff Algorithm", 
       subtitle = "{current_frame}")

# gif画像を出力
gganimate::animate(graph, nframes = iter + 1, fps = 10)

