
# Chapter 7 マルコフ連鎖モンテカルロ法の応用例 ---------------------------------------------

# 7.2 イジング模型 ---------------------------------------------------------------

# 7.2.1 メトロポリス法 -----------------------------------------------------------

# 7.2.1項で利用するパッケージ
library(tidyverse)
library(ggaimate)


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


### 初期値の設定 -----

# 1辺のスピン数を指定
N <- 25

# スピンのマトリクスを初期化
spin_mat <- matrix(sample(x = c(-1, 1), size = N^2, replace = TRUE), nrow = N, ncol = N)

# 作図用にスピンのデータフレームを作成
spin_df <- tibble(
  i_x = rep(1:N, times = N), 
  i_y = rep(1:N, each = N), 
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


### メトロポリス法による更新 -----

# パラメータを指定
J <- 1
h <- 0
temperature <- 1.5

# メトロポリス法
for(iter in 1:100) {
  # 更新するスピン番号を生成
  n_vec <- 1:N^2 # 順番に選択
  #n_vec <- sample(1:N^2, size = N^2, replace = FALSE) # 順番をランダムに選択
  #n_vec <- sample(1:N^2, size = N^2, replace = TRUE) # ランダムに選択
  
  # 1試行における更新回数を初期化
  n_accept <- 0
  
  # 配位を更新
  for(n in n_vec) {
    # スピンのインデックスを計算
    i_y <- ifelse(n %% N == 0, yes = N, no = n %% N) # 行番号
    i_x <- (n - 1) %/% N + 1 # 列番号
    
    # エネルギーを計算:式(7.45)
    energy <- fn_E(spin_mat, J, h)
    
    # i番目のスピンを反転させたマトリクスを作成
    spin_dash_mat <- spin_mat
    spin_dash_mat[i_y, i_x] <- spin_dash_mat[i_y, i_x] * (-1)
    
    # エネルギーを計算:式(7.45)
    energy_dash <- fn_E(spin_dash_mat, J, h)
    
    # 式(4.46')
    p <- exp((energy - energy_dash) / temperature)
    
    # テスト値を生成
    metropolis <- runif(n = 1, min = 0, max = 1)
    
    # メトロポリステストにより配位を更新
    if(p > metropolis) {
      # スピンを更新
      spin_mat <- spin_dash_mat
      
      # 更新回数を加算
      n_accept <- n_accept + 1
    }
    
  }
  
  # スピンのデータフレームを作成
  tmp_spin_df <- tibble(
    i_x = rep(1:N, times = N), 
    i_y = rep(1:N, each = N), 
    spin = as.vector(spin_mat), 
    label = as.factor(paste0("iter:", iter, ", rate:", sum(spin_mat) / N^2))
  )
  
  # 結果を結合
  spin_df <- rbind(spin_df, tmp_spin_df)
  
  # 途中経過を表示
  print(paste0("iter:", iter, ", accept:", n_accept))
}

# 最終結果を作図
spin_df %>% 
  filter(label == as.character(tmp_spin_df[["label"]][1])) %>% # 最後の試行を抽出
  ggplot(aes(x = i_x, y = i_y, fill = spin)) + 
    geom_tile() + # ヒートマップ
    scale_fill_gradient(low = "red", high = "yellow", 
                        breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  labs(title = "Metropolis Method", 
         subtitle = paste0("iter:", iter, ", rate:", sum(spin_mat) / N^2))


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


# 7.2.5 Wolffのアルゴリズム ------------------------------------------------------------

# 7.2.5項で利用するパッケージ
library(tidyverse)
library(gganimate)


### 初期値の設定 -----

# 1辺のスピン数を指定
N <- 250

# スピンのマトリクスを初期化
spin_mat <- matrix(sample(x = c(-1, 1), size = N^2, replace = TRUE), nrow = N, ncol = N)

# 作図用にスピンのデータフレームを作成
spin_df <- tibble(
  i_x = rep(1:N, times = N), 
  i_y = rep(1:N, each = N), 
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


### Wolffのアルゴリズムによる更新 -----

# パラメータを指定
J <- 1
h <- 0
temperature <- 1

# クラスタの追加判定用の確率値を計算
prob <- 1 - exp(-2 * J / temperature)


# 試行回数を初期化
iter <- 0

# Wolffアルゴリズム
while(abs(sum(spin_mat)) / N^2 < 0.95) { # 指定したレートに達するまで
  # 試行回数を加算
  iter <- iter + 1

  # クラスタのマトリクスを初期化
  cluster_mat <- matrix(1, nrow = N, ncol = N)
  
  # 起点のクラスタのインデックスを生成
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
    
    # 右隣を処理
    if(i_y_plus <= N) { # 枠内である
      if(spin_mat[i_y_plus, i_x] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[i_y_plus, i_x] == 1) { # 未処理である
          # 判定用の確率値を生成
          p <- sample(seq(0, 1, 0.001), size = 1)
          if(p < prob) { # 確率で判定
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
    
    # 上隣を処理
    if(i_x_plus <= N) { # 枠内である
      if(spin_mat[i_y, i_x_plus] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[i_y, i_x_plus] == 1) { # 未処理である
          # 判定用の確率値を生成
          p <- sample(seq(0, 1, 0.001), size = 1)
          if(p < prob) { # 確率で判定
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
    
    # 左隣を処理
    if(i_y_minus >= 1) { # 枠内である
      if(spin_mat[i_y_minus, i_x] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[i_y_minus, i_x] == 1) { # 未処理である
          # 判定用の確率値を生成
          p <- sample(seq(0, 1, 0.001), size = 1)
          if(p < prob) { # 確率で判定
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
    if(i_x_minus >= 1) { # 枠内である
      if(spin_mat[i_y, i_x_minus] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[i_y, i_x_minus] == 1) { # 未処理である
          p <- sample(seq(0, 1, 0.001), size = 1)
          if(p < prob) { # 確率で判定
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
    
  }
  
  # クラスタのスピンを反転
  spin_mat[cluster_mat == 0] <- spin_cluster * (-1)
  
  # スピンのデータフレームを作成
  tmp_spin_df <- tibble(
    i_x = rep(1:N, times = N), 
    i_y = rep(1:N, each = N), 
    spin = as.vector(spin_mat), 
    label = as.factor(paste0("iter:", iter, ", rate:", sum(spin_mat) / N^2))
  )
  
  # 結果を結合
  spin_df <- rbind(spin_df, tmp_spin_df)
  
  # 結果を表示
  print(paste0(
    "iter:", iter, ", cluster:", k, ", rate:", sum(spin_mat) / N^2
  ))
}

# 最終結果を作図
spin_df %>% 
  filter(label == as.character(tmp_spin_df[["label"]][1])) %>% # 最後の試行を抽出
  ggplot(aes(x = i_x, y = i_y, fill = spin)) + 
  geom_tile() + # ヒートマップ
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + # タイルの色
  coord_fixed(ratio = 1) + # アスペクト比
  labs(title = "Wolff Algorithm", 
       subtitle = paste0("iter:", iter, ", rate:", sum(spin_mat) / N^2))


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

