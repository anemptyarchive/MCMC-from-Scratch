
# ch7.2 イジング模型 ---------------------------------------------------------------

# 7.2節で利用するパッケージ
library(tidyverse)


# 1辺のスピン数を指定
N <- 10

# スピンのマトリクスを初期化
spin_mat <- matrix(sample(x = c(-1, 1), size = N^2, replace = TRUE), nrow = N, ncol = N)

# スピンのデータフレームを作成
spin_df <- tibble(
  x = rep(1:N, times = N), 
  y = rep(1:N, each = N), 
  spin = as.vector(spin_mat)
)

# 作図
ggplot(spin_df, aes(x = x, y = y, fill = spin)) + 
  geom_tile() + 
  scale_fill_gradient(low = "red", high = "yellow", 
                       breaks = c(-1, 1), guide = guide_legend())


# パラメータを指定
J <- 1
h <- 0
t <- 1


# クラスタのマトリクスを初期化
cluster_mat <- matrix(1, nrow = N, ncol = N)

# 初回のクラスタのインデックスを生成
x_idx <- sample(x = 1:N, size = 1)
y_idx <- sample(x = 1:N, size = 1)

# 初回のクラスタを設定
cluster_mat[x_idx, y_idx] <- 0

# 初回のクラスタのスピンを取得
spin_cluster <- spin_mat[x_idx, y_idx]

# クラスタのインデックスを記録
cluster_idx_mat <- matrix(c(x_idx, y_idx), nrow = 1, ncol = 2)

# クラスタ数をカウント
n_cluster <- 1

# 追加判定用の確率値を計算
prob <- 1 - exp(-2 * J / t)

# 処理済みのクラスタ数
k <- 1
while(k <= n_cluster) {
  # クラスタのインデックスを取得
  x_idx <- cluster_idx_mat[k, 1]
  y_idx <- cluster_idx_mat[k, 2]
  
  # 隣接するインデックスを計算
  x_idx_plus <- x_idx + 1
  y_idx_plus <- y_idx + 1
  x_idx_minus <- x_idx - 1
  y_idx_minus <- y_idx - 1
  
  # 右隣を処理
  if(x_idx_plus <= N) { # 枠内である
    if(spin_mat[x_idx_plus, y_idx] == spin_cluster) { # 現クラスタと同じスピンである
      if(cluster_mat[x_idx_plus, y_idx] == 1) { # 未処理である
        # 判定用の確率値を生成
        p <- sample(seq(0, 1, 0.001), size = 1)
        if(p < prob) { # 確率で判定
          # 新たなクラスタのインデックスを記録
          cluster_idx_mat <- rbind(cluster_idx_mat, c(x_idx_plus, y_idx))
          
          # クラスタ数を加算
          n_cluster <- n_cluster + 1
          
          # クラスタを追加
          cluster_mat[x_idx_plus, y_idx] <- 0
        }
      }
    }
  }
  
  # 上隣を処理
  if(y_idx_plus <= N) { # 枠内である
    if(spin_mat[x_idx, y_idx_plus] == spin_cluster) { # 現クラスタと同じスピンである
      if(cluster_mat[x_idx, y_idx_plus] == 1) { # 未処理である
        # 判定用の確率値を生成
        p <- sample(seq(0, 1, 0.001), size = 1)
        if(p < prob) { # 確率で判定
          # 新たなクラスタのインデックスを記録
          cluster_idx_mat <- rbind(cluster_idx_mat, c(x_idx, y_idx_plus))
          
          # クラスタ数を加算
          n_cluster <- n_cluster + 1
          
          # クラスタを追加
          cluster_mat[x_idx, y_idx_plus] <- 0
        }
      }
    }
  }
  
  # 左隣を処理
  if(x_idx_minus >= 1) { # 枠内である
    if(spin_mat[x_idx_minus, y_idx] == spin_cluster) { # 現クラスタと同じスピンである
      if(cluster_mat[x_idx_minus, y_idx] == 1) { # 未処理である
        # 判定用の確率値を生成
        p <- sample(seq(0, 1, 0.001), size = 1)
        if(p < prob) { # 確率で判定
          # 新たなクラスタのインデックスを記録
          cluster_idx_mat <- rbind(cluster_idx_mat, c(x_idx_minus, y_idx))
          
          # クラスタ数を加算
          n_cluster <- n_cluster + 1
          
          # クラスタを追加
          cluster_mat[x_idx_minus, y_idx] <- 0
        }
      }
    }
  }
  
  # 下隣を処理
  if(y_idx_minus >= 1) { # 枠内である
    if(spin_mat[x_idx, y_idx_minus] == spin_cluster) { # 現クラスタと同じスピンである
      if(cluster_mat[x_idx, y_idx_minus] == 1) { # 未処理である
        p <- sample(seq(0, 1, 0.001), size = 1)
        if(p < prob) { # 確率で判定
          # 新たなクラスタのインデックスを記録
          cluster_idx_mat <- rbind(cluster_idx_mat, c(x_idx, y_idx_minus))
          
          # クラスタ数を加算
          n_cluster <- n_cluster + 1
          
          # クラスタを追加
          cluster_mat[x_idx, y_idx_minus] <- 0
        }
      }
    }
  }
  
  # 処理済みのクラスタ数を加算
  k <- k + 1
  print(paste0("k=", k, " (", k, " / ", n_cluster, ")"))
}

# スピンのマトリクスを更新:(スピンを反転)
spin_mat[cluster_mat == 0] <- spin_cluster * (-1)

# スピンのデータフレームを更新
spin_df <- tibble(
  x = rep(1:N, times = N), 
  y = rep(1:N, each = N), 
  spin = as.vector(spin_mat)
)

# 作図
ggplot(spin_df, aes(x = x, y = y, fill = spin)) + 
  geom_tile() + 
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend())


# ・アニメーション ----------------------------------------------------------------

# 利用するパッケージ
library(tidyverse)
library(gganimate)


# 1辺のスピン数を指定
N <- 100

# スピンのマトリクスを初期化
spin_mat <- matrix(sample(x = c(-1, 1), size = N^2, replace = TRUE), nrow = N, ncol = N)

# スピンのデータフレームを作成
spin_df <- tibble(
  x = rep(1:N, times = N), 
  y = rep(1:N, each = N), 
  spin = as.vector(spin_mat), 
  label = as.factor(paste0("i=", 0))
)

# 作図
ggplot(spin_df, aes(x = x, y = y, fill = spin)) + 
  geom_tile() + 
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + 
  labs(subtitle = paste0("i=", 0))


# パラメータを指定
J <- 1
h <- 0
t <- 1

# 追加判定用の確率値を計算
prob <- 1 - exp(-2 * J / t)


# 試行回数を初期化
i <- 0

# イジング模型
while(abs(sum(spin_mat)) / N^2 < 0.95) {
  # 試行回数を加算
  i <- i + 1

  # クラスタのマトリクスを初期化
  cluster_mat <- matrix(1, nrow = N, ncol = N)
  
  # 起点のクラスタのインデックスを生成
  x_idx <- sample(x = 1:N, size = 1)
  y_idx <- sample(x = 1:N, size = 1)
  
  # 起点のクラスタを設定
  cluster_mat[x_idx, y_idx] <- 0
  
  # 起点のクラスタのスピンを取得
  spin_cluster <- spin_mat[x_idx, y_idx]
  
  # クラスタのインデックスを記録
  cluster_idx_mat <- matrix(c(x_idx, y_idx), nrow = 1, ncol = 2)
  
  # クラスタ数をカウントを初期化
  n_cluster <- 1
  
  # 処理済みのクラスタ数を初期化
  k <- 0
  
  # 試行
  while(k < n_cluster) {
    # 処理済みのクラスタ数を加算
    k <- k + 1
    
    # クラスタのインデックスを取得
    x_idx <- cluster_idx_mat[k, 1]
    y_idx <- cluster_idx_mat[k, 2]
    
    # 隣接するインデックスを計算
    x_idx_plus <- x_idx + 1
    y_idx_plus <- y_idx + 1
    x_idx_minus <- x_idx - 1
    y_idx_minus <- y_idx - 1
    
    # 右隣を処理
    if(x_idx_plus <= N) { # 枠内である
      if(spin_mat[x_idx_plus, y_idx] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[x_idx_plus, y_idx] == 1) { # 未処理である
          # 判定用の確率値を生成
          p <- sample(seq(0, 1, 0.001), size = 1)
          if(p < prob) { # 確率で判定
            # 新たなクラスタのインデックスを記録
            cluster_idx_mat <- rbind(cluster_idx_mat, c(x_idx_plus, y_idx))
            
            # クラスタ数を加算
            n_cluster <- n_cluster + 1
            
            # クラスタを追加
            cluster_mat[x_idx_plus, y_idx] <- 0
          }
        }
      }
    }
    
    # 上隣を処理
    if(y_idx_plus <= N) { # 枠内である
      if(spin_mat[x_idx, y_idx_plus] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[x_idx, y_idx_plus] == 1) { # 未処理である
          # 判定用の確率値を生成
          p <- sample(seq(0, 1, 0.001), size = 1)
          if(p < prob) { # 確率で判定
            # 新たなクラスタのインデックスを記録
            cluster_idx_mat <- rbind(cluster_idx_mat, c(x_idx, y_idx_plus))
            
            # クラスタ数を加算
            n_cluster <- n_cluster + 1
            
            # クラスタを追加
            cluster_mat[x_idx, y_idx_plus] <- 0
          }
        }
      }
    }
    
    # 左隣を処理
    if(x_idx_minus >= 1) { # 枠内である
      if(spin_mat[x_idx_minus, y_idx] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[x_idx_minus, y_idx] == 1) { # 未処理である
          # 判定用の確率値を生成
          p <- sample(seq(0, 1, 0.001), size = 1)
          if(p < prob) { # 確率で判定
            # 新たなクラスタのインデックスを記録
            cluster_idx_mat <- rbind(cluster_idx_mat, c(x_idx_minus, y_idx))
            
            # クラスタ数を加算
            n_cluster <- n_cluster + 1
            
            # クラスタを追加
            cluster_mat[x_idx_minus, y_idx] <- 0
          }
        }
      }
    }
    
    # 下隣を処理
    if(y_idx_minus >= 1) { # 枠内である
      if(spin_mat[x_idx, y_idx_minus] == spin_cluster) { # 現クラスタと同じスピンである
        if(cluster_mat[x_idx, y_idx_minus] == 1) { # 未処理である
          p <- sample(seq(0, 1, 0.001), size = 1)
          if(p < prob) { # 確率で判定
              # 新たなクラスタのインデックスを記録
            cluster_idx_mat <- rbind(cluster_idx_mat, c(x_idx, y_idx_minus))
          
            # クラスタ数を加算
            n_cluster <- n_cluster + 1
            
            # クラスタを追加
            cluster_mat[x_idx, y_idx_minus] <- 0
          }
        }
      }
    }
    
  }
  
  # スピンのマトリクスを更新:(スピンを反転)
  spin_mat[cluster_mat == 0] <- spin_cluster * (-1)
  
  # スピンのデータフレームを作成
  tmp_spin_df <- tibble(
    x = rep(1:N, times = N), 
    y = rep(1:N, each = N), 
    spin = as.vector(spin_mat), 
    label = as.factor(paste0("i=", i))
  )
  
  # 結果を結合
  spin_df <- rbind(spin_df, tmp_spin_df)
  
  # 結果を表示
  print(paste0(
    "i=", i, ", k=", k, ", rate:", round(abs(sum(spin_mat)) / N^2, 3)
  ))
}


# 作図
graph <- ggplot(spin_df, aes(x = x, y = y, fill = spin)) + 
  geom_tile() + 
  scale_fill_gradient(low = "red", high = "yellow", 
                      breaks = c(-1, 1), guide = guide_legend()) + 
  gganimate::transition_manual(label) + 
  labs(subtitle = "{current_frame}")

# gif画像を出力
gganimate::animate(graph, nframes = i + 1, fps = 10)

