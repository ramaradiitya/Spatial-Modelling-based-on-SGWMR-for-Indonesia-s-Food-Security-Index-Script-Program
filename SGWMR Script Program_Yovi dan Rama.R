library(lmtest) #Untuk Model Linier test
library(spdep) #Untuk Spatial Dependency
library(car) #Untuk Multivariat Multiple Regression  
library(psych) #Untuk uji Bartlett Test
library(MVN) #Untuk uji multivariate normality
library(MASS)   # Untuk fungsi ginv() (Generalized Inverse)
library(readxl) # Untuk membaca file Excel

# ----------------------------------------------------
# IMPOR DATASET
# ----------------------------------------------------
df <- read_excel("C:/Users/LENOVO/Downloads/Dataset IKP Besar.xlsx", sheet = "Dataset")

# Mapping variabel
df$Y1 <- df$`Indeks Ketersediaan`
df$Y2 <- df$`Indeks Keterjangkauan`
df$Y3 <- df$`Indeks Pemanfaatan`
df$X1 <- df$`Produksi Beras Per Kapita (Ton/Ribu Jiwa)`
df$X2 <- df$`PDRB atas Harga Berlaku`
df$X3 <- df$`Persentase Air Minum Layak`
df$X4 <- df$`PoU`

n <- nrow(df)
m <- 3 # Jumlah variabel respons (Y1, Y2, Y3)

# ----------------------------------------------------
# PRA-PEMROSESAN (STANDARDISASI & PEMBENTUKAN MATRIKS)
# ----------------------------------------------------
# Standardisasi variabel X.
X1_s <- scale(df$X1)
X2_s <- scale(df$X2)
X3_s <- scale(df$X3)
X4_s <- scale(df$X4)
# Menyusun Matriks Desain (X) dengan tambahan kolom Intersep
X <- cbind(1, X1_s, X2_s, X3_s, X4_s)
colnames(X) <- c("Intercept", "X1", "X2", "X3", "X4")

# Menyusun Matriks Respons (Y)
Y <- as.matrix(df[, c("Y1", "Y2", "Y3")])

# Menyiapkan data frame untuk fungsi lm()
df_model <- data.frame(Y1=Y[,1], Y2=Y[,2], Y3=Y[,3], 
                       X1=X1_s, X2=X2_s, X3=X3_s, X4=X4_s)

# Estimasi Model Global MMR menggunakan OLS
fit_MMR <- lm(cbind(Y1, Y2, Y3) ~ X1 + X2 + X3 + X4, data = df_model)
summary(fit_MMR)

cat("\n=== UJI ASUMSI KLASIK MODEL GLOBAL (MMR) ===\n")
residual_MMR <- residuals(fit_MMR)

# --------------------------------
# 1. UJI MULTIKOLINIERITAS (VIF)
# --------------------------------
cat("\n1. UJI MULTIKOLINIERITAS (VIF)\n")
fit_vif <- lm(Y1 ~ X1 + X2 + X3 + X4, data = df_model)
fit_vif1 <- lm(Y2 ~ X1 + X2 + X3 + X4, data = df_model)
fit_vif2 <- lm(Y3 ~ X1 + X2 + X3 + X4, data = df_model)
vif(fit_vif)
vif(fit_vif1)
vif(fit_vif2)

# ----------------------------------
# UJI BARTLETT'S TEST OF SPHERICITY
# ----------------------------------
cat("\n UJI BARTLETT'S TEST OF SPHERICITY (Korelasi antar Respons)\n")
# Uji ini diterapkan pada matriks korelasi dari residual MMR
# H0: Matriks korelasi identitas (Tidak ada korelasi antar respons Y)
# H1: Matriks korelasi BUKAN identitas (Ada korelasi antar respons Y)

cor_matrix_resid <- cor(residual_MMR)
bartlett_test <- cortest.bartlett(cor_matrix_resid, n = nrow(df_model))
cor_matrix_Y <- cor(df_model[, c("Y1","Y2","Y3")])
bartlett_Y <- cortest.bartlett(cor_matrix_Y, n = nrow(df_model))
print(bartlett_Y)
cat("Chi-Square Stat :", round(bartlett_test$chisq, 4), "\n")
cat("P-value         :", round(bartlett_test$p.value, 5), "\n")
cat("Interpretasi    : Jika p-value < 0.05, maka H0 ditolak. Artinya terdapat korelasi antar\n",
    "                 respons (Y1, Y2, Y3), sehingga pemodelan MULTIVARIAT sangat tepat digunakan.\n")

# ------------------------------
# UJI NORMALITAS MULTIVARIAT
# ------------------------------
cat("\n UJI NORMALITAS MULTIVARIAT RESIDUAL\n")
# H0: Residual berdistribusi normal multivariat
# H1: Residual tidak berdistribusi normal multivariat
Henze_Zirkler <- mvn(data = as.data.frame(residual_MMR))
print(Henze_Zirkler$multivariate_normality)

cat("\n=== UJI DIAGNOSTIK RESIDUAL MODEL GLOBAL (MMR) ===\n")
# 1. Membentuk Matriks Bobot Spasial untuk Moran's I
# Menggunakan metode k-Nearest Neighbors (misal k=4 tetangga terdekat)
k_neighbors <- 3
coords <- as.matrix(df[, c("Longitude", "Latitude")])
nb <- knn2nb(knearneigh(coords, k = k_neighbors))
listw <- nb2listw(nb, style = "W") # Row-standardized weights
nb_check <- knn2nb(knearneigh(coords, k = 2))
for(i in 1:38) {
  cat(df$Provinsi[i], "->", df$Provinsi[nb_check[[i]]], "\n")
}

# Menyiapkan data frame untuk rekapitulasi hasil uji
diagnostik_results <- data.frame(
  Respons = character(),
  BP_Stat = numeric(),
  BP_PValue = numeric(),
  Moran_I_Stat = numeric(),
  Moran_PValue = numeric(),
  stringsAsFactors = FALSE
)

# 2. Eksekusi Uji untuk Masing-masing Variabel Respons (Y1, Y2, Y3)
respons_names <- c("Y1", "Y2", "Y3")

for (i in 1:m) {
  # Membentuk model regresi univariat (karena uji spdep/lmtest butuh objek 1D)
  formula_uni <- as.formula(paste(respons_names[i], "~ X1 + X2 + X3 + X4"))
  fit_uni <- lm(formula_uni, data = df_model)
  
  # A. Uji Breusch-Pagan (Heteroskedastisitas)
  # H0: Varians residual konstan (Homoskedastis)
  # H1: Terdapat heteroskedastisitas
  bp_test <- bptest(fit_uni)
  
  # B. Uji Moran's I pada Residual (Autokorelasi Spasial)
  # H0: Tidak ada autokorelasi spasial pada residual
  # H1: Terdapat autokorelasi spasial
  moran_test <- lm.morantest(fit_uni, listw)
  
  # Menyimpan hasil ke dalam data frame
  diagnostik_results <- rbind(diagnostik_results, data.frame(
    Respons = respons_names[i],
    BP_Stat = round(bp_test$statistic, 4),
    BP_PValue = round(bp_test$p.value, 5),
    Moran_I_Stat = round(moran_test$estimate[1], 4),
    Moran_PValue = round(moran_test$p.value, 5)
  ))
}

# Menampilkan hasil kompilasi uji diagnostik
print(diagnostik_results)

# Memberikan interpretasi otomatis
cat("\nInterpretasi Diagnostik:\n")
for(i in 1:nrow(diagnostik_results)) {
  cat("Persamaan", diagnostik_results$Respons[i], ":\n")
  if(diagnostik_results$BP_PValue[i] < 0.05) {
    cat(" - Terdapat Heteroskedastisitas (Breusch-Pagan p < 0.05).\n")
  } else {
    cat(" - Varians Konstan (Breusch-Pagan p >= 0.05).\n")
  }
  
  if(diagnostik_results$Moran_PValue[i] < 0.05) {
    cat(" - Terdapat Autokorelasi Spasial (Moran's I p < 0.05).\n")
  } else {
    cat(" - Tidak ada Autokorelasi Spasial (Moran's I p >= 0.05).\n")
  }
}

cat("\nKesimpulan: Jika terdapat p-value < 0.05 pada salah satu uji di atas, penggunaan model GWMR sangat direkomendasikan.\n\n")

# GOODNESS OF FIT MMR MODEL
# Menghitung Residual Sum of Squares (RSS) untuk MMR
E_MMR <- residuals(fit_MMR)
RSS_MMR <- sum(E_MMR^2)

# Menghitung AICc untuk MMR
Sigma_MMR <- (t(E_MMR) %*% E_MMR) / n
det_Sigma_MMR <- det(Sigma_MMR)
p_params <- 5 # Jumlah parameter per respons (Intercept + 3 Prediktor)

AICc_MMR <- (n * log(det_Sigma_MMR)) + (n * m * log(2 * pi)) + 
  (n * (n + p_params) * m) / (n - p_params - m - 1)
AICc_MMR


# ===============================================
# ALGORITMA GWMR UNTUK 3 RESPONS DAN 4 PREDIKTOR
# Kombinasi: PDRB, Produksi Beras, Air Minum Layak, PoU
# ===============================================
library(MASS)
library(readxl)
df <- read_excel("C:/Users/LENOVO/Downloads/Dataset IKP Besar.xlsx", sheet = "Dataset")

# Mapping variabel
df$Y1 <- df$`Indeks Ketersediaan`
df$Y2 <- df$`Indeks Keterjangkauan`
df$Y3 <- df$`Indeks Pemanfaatan`
df$X1 <- df$`Produksi Beras Per Kapita (Ton/Ribu Jiwa)`
df$X2 <- df$`PDRB atas Harga Berlaku`
df$X3 <- df$`Persentase Air Minum Layak`
df$X4 <- df$`PoU`

n <- nrow(df)
m <- 3

#
X1_s <- scale(df$X1)   
X2_s <- scale(df$X2)   
X3_s <- scale(df$X3)   
X4_s <- scale(df$X4)   

X <- cbind(1, X1_s, X2_s, X3_s, X4_s)
colnames(X) <- c("Intercept", "X1", "X2", "X3", "X4")

Y <- as.matrix(df[, c("Y1", "Y2", "Y3")])
coords <- as.matrix(df[, c("Longitude", "Latitude")]) 
dist_matrix <- as.matrix(dist(coords))

adaptive_bisquare <- function(distance_vec, h) {
  weight_vec <- ifelse(distance_vec <= h, (1 - (distance_vec / h)^2)^2, 0)
  return(diag(weight_vec))
}
# ------------------------------------------------------------------------------
# 3. OPTIMASI BANDWIDTH (MEMINIMALKAN AICc)
# ------------------------------------------------------------------------------
calculate_aicc <- function(h) {
  tr_S <- 0
  E_mat <- matrix(0, nrow = n, ncol = m)
  
  for (i in 1:n) {
    W_i <- adaptive_bisquare(dist_matrix[i, ], h)
    
    if (sum(diag(W_i)) <= 4) return(Inf) 
    
    X_Wi_X_inv <- tryCatch(ginv(t(X) %*% W_i %*% X), error = function(e) return(NULL))
    if (is.null(X_Wi_X_inv)) return(Inf)
    
    S_ii <- (X[i, , drop=FALSE] %*% X_Wi_X_inv %*% t(X[i, , drop=FALSE])) * W_i[i, i]
    tr_S <- tr_S + as.numeric(S_ii)
    
    B_hat_i <- X_Wi_X_inv %*% t(X) %*% W_i %*% Y
    Y_hat_i <- X[i, , drop=FALSE] %*% B_hat_i
    E_mat[i, ] <- Y[i, ] - Y_hat_i
  }
  
  Sigma_hat <- (t(E_mat) %*% E_mat) / n
  det_Sigma <- det(Sigma_hat)
  
  if (det_Sigma <= 0) return(Inf)
  
  AICc <- (n * log(det_Sigma)) + (n * m * log(2 * pi)) + 
    (n * (n + tr_S) * m) / (n - tr_S - m - 1)
  return(list(AICc = AICc, tr_S = tr_S))
}

cat("Mencari bandwidth optimal...\n")
h_candidates <- seq(max(dist_matrix) * 0.1, max(dist_matrix), length.out = 30)
best_h <- NULL
best_AICc <- Inf
global_eta <- NULL

for (h in h_candidates) {
  res <- calculate_aicc(h)
  if (is.list(res) && res$AICc < best_AICc) {
    best_AICc <- res$AICc
    best_h <- h
    global_eta <- res$tr_S
  }
}

cat("Bandwidth Optimal (h):", round(best_h, 4), "\n")
cat("AICc Terkecil:", round(best_AICc, 4), "\n")
cat("Effective Parameters (eta):", round(global_eta, 4), "\n\n")

# ------------------------------------------------------------------------------
# 4. FUNGSI UJI HIPOTESIS MGLH LOKAL
# ------------------------------------------------------------------------------
gwmr_mglh_local <- function(X, Y, Wi, L, M, C, eta) {
  X_Wi_X <- t(X) %*% Wi %*% X
  inv_X_Wi_X <- ginv(X_Wi_X) 
  B_hat <- inv_X_Wi_X %*% t(X) %*% Wi %*% Y
  
  L_B_M_minus_C <- (L %*% B_hat %*% M) - C
  H_mat <- t(L_B_M_minus_C) %*% ginv(L %*% inv_X_Wi_X %*% t(L)) %*% L_B_M_minus_C
  
  E_mat <- t(M) %*% (t(Y) %*% Wi %*% Y - t(B_hat) %*% X_Wi_X %*% B_hat) %*% M
  
  V <- sum(diag(H_mat %*% ginv(H_mat + E_mat)))
  
  a <- nrow(L); b <- ncol(M); q <- min(a, b)
  theta <- (abs(a - b) - 1) / 2
  lambda <- (eta - b - 1) / 2
  
  num <- 2 * lambda + q + 1
  den <- 2 * theta + q + 1
  
  if(abs(q - V) < 1e-10) {
    p_val <- NA
  } else {
    F_star <- (num / den) * (V / (q - V))
    df1 <- a * b; df2 <- q * num
    p_val <- pf(F_star, df1, df2, lower.tail = FALSE)
  }
  
  return(list(B_hat = B_hat, p_value = p_val))
}

# ------------------------------------------------------------------------------
# 5. EKSEKUSI PENGUJIAN MGLH (PARSIAL UNTUK SELURUH KOMBINASI X DAN Y)
# ------------------------------------------------------------------------------
cat("\n=== TAHAP 5: PENGUJIAN HIPOTESIS PARSIAL MGLH ===\n")

# Menyiapkan matriks kosong untuk menyimpan 12 koefisien dan 12 p-value tiap lokasi
n_combinations <- 4 * 3
coef_matrix <- matrix(0, nrow = n, ncol = n_combinations)
pval_matrix <- matrix(0, nrow = n, ncol = n_combinations)

col_names <- c()
for(p in 1:4) {
  for(r in 1:3) {
    col_names <- c(col_names, paste0("X", p, "_Y", r))
  }
}
colnames(coef_matrix) <- paste0("Coef_", col_names)
colnames(pval_matrix) <- paste0("Pval_", col_names)

cat("Memproses pengujian untuk", n, "lokasi...\n")

for (i in 1:n) {
  Wi <- adaptive_bisquare(dist_matrix[i, ], best_h)
  X_Wi_X <- t(X) %*% Wi %*% X
  inv_X_Wi_X <- ginv(X_Wi_X) 
  B_hat <- inv_X_Wi_X %*% t(X) %*% Wi %*% Y
  Y_Wi_Y <- t(Y) %*% Wi %*% Y
  B_X_Wi_X_B <- t(B_hat) %*% X_Wi_X %*% B_hat
  
  col_idx <- 1
  
  # 2. Iterasi untuk 4 Prediktor dan 3 Respons
  for (p in 1:4) {
    for (r in 1:3) {
      # H0: beta_pr = 0
      L_mat <- matrix(0, nrow = 1, ncol = 5)
      L_mat[1, p + 1] <- 1  # Index p+1 karena kolom 1 adalah Intersep
      
      M_mat <- matrix(0, nrow = 3, ncol = 1)
      M_mat[r, 1] <- 1
      
      C_mat <- matrix(0, nrow = 1, ncol = 1)
      
      L_B_M_minus_C <- (L_mat %*% B_hat %*% M_mat) - C_mat
      H_mat <- t(L_B_M_minus_C) %*% ginv(L_mat %*% inv_X_Wi_X %*% t(L_mat)) %*% L_B_M_minus_C
      E_mat <- t(M_mat) %*% (Y_Wi_Y - B_X_Wi_X_B) %*% M_mat
      
      V <- sum(diag(H_mat %*% ginv(H_mat + E_mat)))
      
      a <- 1; b <- 1; q <- 1; theta <- 0
      lambda <- (global_eta - b - 1) / 2
      num <- 2 * lambda + q + 1
      den <- 2 * theta + q + 1
      
      if (abs(q - V) < 1e-10) {
        p_val <- NA
      } else {
        F_star <- (num / den) * (V / (q - V))
        p_val <- pf(F_star, df1 = a * b, df2 = q * num, lower.tail = FALSE)
      }
      
      coef_matrix[i, col_idx] <- B_hat[p + 1, r]
      pval_matrix[i, col_idx] <- p_val
      
      col_idx <- col_idx + 1
    }
  }
}

# ------------------------------------------------------------------------------
# 6. OUTPUT HASIL AKHIR
# ------------------------------------------------------------------------------
output_final <- data.frame(Provinsi = df$Provinsi)
output_final <- cbind(output_final, round(coef_matrix, 4), round(pval_matrix, 5))

print(output_final)

# write.csv(output_final, "Hasil_Estimasi_MGLH_GWMR.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. PERBANDINGAN DENGAN MODEL GLOBAL (MMR)
# ------------------------------------------------------------------------------
df_model <- data.frame(Y1=Y[,1], Y2=Y[,2], Y3=Y[,3], 
                       X1=X1_s, X2=X2_s, X3=X3_s, X4=X4_s)

fit_MMR <- lm(cbind(Y1, Y2, Y3) ~ X1 + X2 + X3 + X4, data = df_model)

E_MMR <- residuals(fit_MMR)
RSS_MMR <- sum(E_MMR^2)

Sigma_MMR <- (t(E_MMR) %*% E_MMR) / n
det_Sigma_MMR <- det(Sigma_MMR)
p_params <- 5 # Jumlah parameter per respons (Intercept + 4 Prediktor)

AICc_MMR <- (n * log(det_Sigma_MMR)) + (n * m * log(2 * pi)) + 
  (n * (n + p_params) * m) / (n - p_params - m - 1)

# Menghitung RSS untuk GWMR (Berdasarkan bandwidth optimal)
E_GWMR <- matrix(0, n, m)
B_GWMR_array <- array(0, dim=c(n, p_params, m))

for(i in 1:n) {
  W_i <- adaptive_bisquare(dist_matrix[i, ], best_h)
  B_hat_i <- ginv(t(X) %*% W_i %*% X) %*% t(X) %*% W_i %*% Y
  
  B_GWMR_array[i, , ] <- B_hat_i
  
  Y_hat_i <- X[i, , drop=FALSE] %*% B_hat_i
  E_GWMR[i, ] <- Y[i, ] - Y_hat_i
}

RSS_GWMR <- sum(E_GWMR^2)

cat("RSS MMR Global   :", round(RSS_MMR, 4), "\n")
cat("RSS GWMR Lokal   :", round(RSS_GWMR, 4), "\n")
cat("AICc MMR Global  :", round(AICc_MMR, 4), "\n")
cat("AICc GWMR Lokal  :", round(best_AICc, 4), "\n")

# -----------------------------------------------------------------
# UJI HETEROGENITAS SPASIAL (MONTE CARLO TEST UNTUK SEMUA VARIABEL)
# -----------------------------------------------------------------
var_orig <- matrix(0, nrow = 4, ncol = 3)
rownames(var_orig) <- c("X1", "X2", "X3", "X4")
colnames(var_orig) <- c("Y1", "Y2", "Y3")

for(p in 1:4) {
  for(r in 1:3) {
    var_orig[p, r] <- var(B_GWMR_array[, p + 1, r]) 
  }
}

B_sim <- 999
var_perm_array <- array(0, dim = c(B_sim, 4, 3))

set.seed(2026)
cat("Menjalankan", B_sim, "permutasi spasial untuk", 4*3, "kombinasi variabel...\n")

for(b in 1:B_sim) {
  idx_perm <- sample(1:n)
  coords_perm <- coords[idx_perm, ]
  dist_perm <- as.matrix(dist(coords_perm))
  
  coef_perm_tmp <- array(0, dim = c(n, 4, 3)) 
  
  for(i in 1:n) {
    W_i_perm <- adaptive_bisquare(dist_perm[i, ], best_h)
    XWX_inv <- tryCatch(ginv(t(X) %*% W_i_perm %*% X), error=function(e) NULL)
    
    if(!is.null(XWX_inv)) {
      B_hat_perm <- XWX_inv %*% t(X) %*% W_i_perm %*% Y
      # Ekstrak baris prediktor (2 s.d 5) untuk semua kolom respons (1 s.d 3)
      coef_perm_tmp[i, , ] <- B_hat_perm[2:5, ] 
    }
  }
  
  for(p in 1:4) {
    for(r in 1:3) {
      var_perm_array[b, p, r] <- var(coef_perm_tmp[, p, r], na.rm = TRUE)
    }
  }
  
  if(b %% 10 == 0) cat("Iterasi", b, "/", B_sim, "selesai...\n")
}

result_mc_df <- data.frame(
  Predictor = character(),
  Response = character(),
  Original_Var = numeric(),
  P_Value = numeric(),
  Status = character(),
  stringsAsFactors = FALSE
)

for(p in 1:4) {
  for(r in 1:3) {
    p_val <- (1 + sum(var_perm_array[, p, r] >= var_orig[p, r])) / (1 + B_sim)
    
    status <- ifelse(p_val < 0.05, 
                     "Signifikan (Heterogen Spasial)", 
                     "Tidak Signifikan (Stasioner/Global)")
    
    result_mc_df <- rbind(result_mc_df, data.frame(
      Predictor = paste0("X", p),
      Response = paste0("Y", r),
      Original_Var = round(var_orig[p, r], 4),
      P_Value = round(p_val, 4),
      Status = status
    ))
  }
}
cat("\n=== HASIL UJI HETEROGENITAS SPASIAL ===\n")
print(result_mc_df)


cat("\n=== EVALUASI KEBAIKAN MODEL (R2, MSE, MAPE) ===\n")
local_R2 <- numeric(n)
local_MSE <- numeric(n)
local_MAPE <- numeric(n)

Y_hat_global <- matrix(0, nrow = n, ncol = m)
colnames(Y_hat_global) <- c("Y1_hat", "Y2_hat", "Y3_hat")

Y_bar <- colMeans(Y)

cat("Menghitung metrik evaluasi untuk", n, "lokasi...\n")

for (i in 1:n) {
  Wi <- adaptive_bisquare(dist_matrix[i, ], best_h)
  X_Wi_X <- t(X) %*% Wi %*% X
  inv_X_Wi_X <- ginv(X_Wi_X) 
  B_hat_i <- inv_X_Wi_X %*% t(X) %*% Wi %*% Y
  
  Y_hat_all_s <- X %*% B_hat_i
  
  Y_hat_i <- X[i, , drop = FALSE] %*% B_hat_i
  Y_hat_global[i, ] <- Y_hat_i
  
  num_R2 <- 0
  den_R2 <- 0
  for(s in 1:n) {
    if(Wi[s, s] > 0) {
      num_R2 <- num_R2 + Wi[s, s] * sum((Y[s, ] - Y_hat_all_s[s, ])^2)
      den_R2 <- den_R2 + Wi[s, s] * sum((Y[s, ] - Y_bar)^2)
    }
  }
  local_R2[i] <- 1 - (num_R2 / den_R2)
  
  local_MSE[i] <- mean((Y[i, ] - Y_hat_i)^2)
  
  local_MAPE[i] <- mean(abs((Y[i, ] - Y_hat_i) / Y[i, ])) * 100
}

global_MSE <- mean((Y - Y_hat_global)^2)
global_MAPE <- mean(abs((Y - Y_hat_global) / Y)) * 100

SST_global <- sum(sweep(Y, 2, Y_bar)^2)
SSE_global <- sum((Y - Y_hat_global)^2)
global_R2 <- 1 - (SSE_global / SST_global)

cat("\n--- RINGKASAN METRIK GLOBAL ---\n")
cat("Global R-Squared :", round(global_R2, 4), "\n")
cat("Global MSE       :", round(global_MSE, 4), "\n")
cat("Global MAPE      :", round(global_MAPE, 2), "%\n")

output_final$Local_R2 <- round(local_R2, 4)
output_final$Local_MSE <- round(local_MSE, 4)
output_final$Local_MAPE <- round(local_MAPE, 2)

print(output_final[, c("Provinsi", "Local_R2", "Local_MSE", "Local_MAPE")])

cat("\n=== PROGRAM SEMIPARAMETRIC GWMR (MIXED GWMR) ===\n")
# ------------------------------------------------------------------------------
cat("\n=== PROGRAM SEMIPARAMETRIC GWMR (MIXED GWMR) ===\n")
# ------------------------------------------------------------------------------
X1_s <- scale(df$X1)
X2_s <- scale(df$X2)
X3_s <- scale(df$X3)
X4_s <- scale(df$X4)
X <- cbind(1, X1_s, X2_s, X3_s, X4_s)
colnames(X) <- c("Intercept", "X1", "X2", "X3", "X4")

# TAMBAHAN 1: Rebuild df_model dari nilai X yang baru
df_model <- data.frame(Y1=Y[,1], Y2=Y[,2], Y3=Y[,3],
                       X1=X1_s, X2=X2_s, X3=X3_s, X4=X4_s)

# TAMBAHAN 2: Cari ulang bandwidth optimal untuk kombinasi X4,6,7,8 ini
# (persis bagian "OPTIMASI BANDWIDTH" di script GWMR — jalankan calculate_aicc(), grid search h)
h_candidates <- seq(max(dist_matrix) * 0.1, max(dist_matrix), length.out = 30)
best_h <- NULL
best_AICc <- Inf
global_eta <- NULL

for (h in h_candidates) {
  res <- calculate_aicc(h)
  if (is.list(res) && res$AICc < best_AICc) {
    best_AICc <- res$AICc
    best_h <- h
    global_eta <- res$tr_S
  }
}
cat("Bandwidth Optimal (h) untuk kombinasi X1,2,3,4:", round(best_h, 4), "\n")
cat("AICc Terkecil:", round(best_AICc, 4), "\n\n")

# TAMBAHAN 3: Hitung ulang RSS_GWMR & B_GWMR_array (dipakai nanti di uji F GWMR vs SGWMR)
E_GWMR <- matrix(0, n, m)
p_params <- 5
B_GWMR_array <- array(0, dim=c(n, p_params, m))

for(i in 1:n) {
  W_i <- adaptive_bisquare(dist_matrix[i, ], best_h)
  B_hat_i <- ginv(t(X) %*% W_i %*% X) %*% t(X) %*% W_i %*% Y
  B_GWMR_array[i, , ] <- B_hat_i
  Y_hat_i <- X[i, , drop=FALSE] %*% B_hat_i
  E_GWMR[i, ] <- Y[i, ] - Y_hat_i
}
RSS_GWMR <- sum(E_GWMR^2)
# TAMBAHAN 4: Hitung ulang MMR dan Global R2 GWMR untuk kombinasi X4,6,7,8
# ------------------------------------------------------------------------------
# 4a. Model MMR (global)
fit_MMR <- lm(cbind(Y1, Y2, Y3) ~ X1 + X2 + X3 + X4, data = df_model)
E_MMR <- residuals(fit_MMR)
RSS_MMR <- sum(E_MMR^2)

Sigma_MMR <- (t(E_MMR) %*% E_MMR) / n
det_Sigma_MMR <- det(Sigma_MMR)

AICc_MMR <- (n * log(det_Sigma_MMR)) + (n * m * log(2 * pi)) + 
  (n * (n + p_params) * m) / (n - p_params - m - 1)

cat("RSS MMR Global   :", round(RSS_MMR, 4), "\n")
cat("AICc MMR Global  :", round(AICc_MMR, 4), "\n\n")

# 4b. Global R2 untuk GWMR
Y_hat_global_gwmr <- matrix(0, nrow = n, ncol = m)
Y_bar <- colMeans(Y)

for (i in 1:n) {
  W_i <- adaptive_bisquare(dist_matrix[i, ], best_h)
  B_hat_i <- ginv(t(X) %*% W_i %*% X) %*% t(X) %*% W_i %*% Y
  Y_hat_global_gwmr[i, ] <- X[i, , drop = FALSE] %*% B_hat_i
}

SST_global_gwmr <- sum(sweep(Y, 2, Y_bar)^2)
SSE_global_gwmr <- sum((Y - Y_hat_global_gwmr)^2)
global_R2 <- 1 - (SSE_global_gwmr / SST_global_gwmr)

cat("Global R2 GWMR   :", round(global_R2, 4), "\n\n")


# 1. SETTING VARIABEL GLOBAL (Z) DAN LOKAL (X)
# ------------------------------------------------------------------------------
# Skenario (berdasarkan hasil Uji Heterogenitas Spasial Monte Carlo):

# Matriks Prediktor Global (Z)
Z_global <- cbind(df_model$X1, df_model$X2, df_model$X4)
colnames(Z_global) <- c("Z_X1", "Z_X2", "Z_X4")

# Matriks Prediktor Lokal (X_L)
X_local <- cbind(1, df_model$X3)
colnames(X_local) <- c("Intercept_Loc", "X3_Loc")

n_loc_vars <- ncol(X_local)
n_glob_vars <- ncol(Z_global)

# ----------------------------------
# LANGKAH 1: ELIMINASI EFEK SPASIAL 
# ----------------------------------
E_Y <- matrix(0, nrow = n, ncol = m)
E_Z <- matrix(0, nrow = n, ncol = n_glob_vars)

for (i in 1:n) {
  W_i <- adaptive_bisquare(dist_matrix[i, ], best_h)
  
  inv_XWX_L <- tryCatch(ginv(t(X_local) %*% W_i %*% X_local), error = function(e) NULL)
  
  if (!is.null(inv_XWX_L)) {
    H_i <- X_local[i, , drop = FALSE] %*% inv_XWX_L %*% t(X_local) %*% W_i
    
    Y_hat_loc_i <- H_i %*% Y
    Z_hat_loc_i <- H_i %*% Z_global
    
    E_Y[i, ] <- Y[i, ] - Y_hat_loc_i
    E_Z[i, ] <- Z_global[i, ] - Z_hat_loc_i
  }
}

# ------------------------------------------------------------------------------
# 3. LANGKAH 2: ESTIMASI PARAMETER GLOBAL (GAMMA)
# ------------------------------------------------------------------------------
cat("Langkah 2: Mengestimasi parameter global...\n")

Gamma_hat <- ginv(t(E_Z) %*% E_Z) %*% t(E_Z) %*% E_Y
rownames(Gamma_hat) <- colnames(Z_global)
colnames(Gamma_hat) <- c("Y1", "Y2", "Y3")

cat("\n--- Koefisien Global (Sama untuk semua wilayah) ---\n")
print(round(Gamma_hat, 4))

# ------------------------------------------------------------------------------
# 4. LANGKAH 3: PEMBENTUKAN RESPONS PARSIAL (Y*)
# ------------------------------------------------------------------------------
cat("\nLangkah 3: Mengekstraksi efek spasial murni (Y*)...\n")

Y_star <- Y - (Z_global %*% Gamma_hat)

# ------------------------------------------------------------------------------
# 5. LANGKAH 4: ESTIMASI PARAMETER LOKAL AKHIR
# ------------------------------------------------------------------------------
cat("Langkah 4: Mengestimasi parameter lokal yang bervariasi secara geografis...\n")

B_local_array <- array(0, dim = c(n, n_loc_vars, m))
dimnames(B_local_array) <- list(df$Provinsi, colnames(X_local), colnames(Y))

for (i in 1:n) {
  W_i <- adaptive_bisquare(dist_matrix[i, ], best_h)
  inv_XWX_L <- tryCatch(ginv(t(X_local) %*% W_i %*% X_local), error = function(e) NULL)
  
  if (!is.null(inv_XWX_L)) {
    B_local_hat_i <- inv_XWX_L %*% t(X_local) %*% W_i %*% Y_star
    B_local_array[i, , ] <- B_local_hat_i
  }
}

# ------------------------------------------------------------------------------
# 6. LANGKAH 5: KALKULASI PREDIKSI SGWMR & EVALUASI MODEL (R2, MSE, MAPE)
# ------------------------------------------------------------------------------
cat("\nLangkah 5: Mengevaluasi Goodness-of-Fit (Lokal & Global)...\n")

Y_hat_sgwmr <- matrix(0, nrow = n, ncol = m)
colnames(Y_hat_sgwmr) <- c("Y1_hat", "Y2_hat", "Y3_hat")

local_R2_sgwmr <- numeric(n)
local_MSE_sgwmr <- numeric(n)
local_MAPE_sgwmr <- numeric(n)

Y_bar <- colMeans(Y)

for (i in 1:n) {
  pred_global_i <- Z_global[i, , drop = FALSE] %*% Gamma_hat
  pred_local_i  <- X_local[i, , drop = FALSE] %*% B_local_array[i, , ]
  
  Y_hat_sgwmr[i, ] <- pred_global_i + pred_local_i
  
  local_MSE_sgwmr[i] <- mean((Y[i, ] - Y_hat_sgwmr[i, ])^2)
  local_MAPE_sgwmr[i] <- mean(abs((Y[i, ] - Y_hat_sgwmr[i, ]) / Y[i, ])) * 100
  
  W_i <- adaptive_bisquare(dist_matrix[i, ], best_h)
  num_R2 <- 0
  den_R2 <- 0
  
  for(s in 1:n) {
    if(W_i[s, s] > 0) {
      Y_hat_s_from_i <- (Z_global[s, , drop=FALSE] %*% Gamma_hat) +
        (X_local[s, , drop=FALSE] %*% B_local_array[i, , ])
      
      num_R2 <- num_R2 + W_i[s, s] * sum((Y[s, ] - Y_hat_s_from_i)^2)
      den_R2 <- den_R2 + W_i[s, s] * sum((Y[s, ] - Y_bar)^2)
    }
  }
  local_R2_sgwmr[i] <- 1 - (num_R2 / den_R2)
}

global_MSE_sgwmr <- mean((Y - Y_hat_sgwmr)^2)
global_MAPE_sgwmr <- mean(abs((Y - Y_hat_sgwmr) / Y)) * 100

SST_global <- sum(sweep(Y, 2, Y_bar)^2)
SSE_global <- sum((Y - Y_hat_sgwmr)^2)
global_R2_sgwmr <- 1 - (SSE_global / SST_global)

cat("\n--- RINGKASAN METRIK GLOBAL SGWMR ---\n")
cat("Global R-Squared :", round(global_R2_sgwmr, 4), "\n")
cat("Global MSE       :", round(global_MSE_sgwmr, 4), "\n")
cat("Global MAPE      :", round(global_MAPE_sgwmr, 2), "%\n")

# ------------------------------------------------------------------------------
# 7. KOMPILASI HASIL AKHIR SGWMR
# ------------------------------------------------------------------------------
sgwmr_results <- data.frame(
  Provinsi = df$Provinsi,
  
  Global_X4_Y1 = rep(Gamma_hat[1, 1], n),
  Global_X4_Y2 = rep(Gamma_hat[1, 2], n),
  Global_X4_Y3 = rep(Gamma_hat[1, 3], n),
  Global_X6_Y1 = rep(Gamma_hat[2, 1], n),
  Global_X6_Y2 = rep(Gamma_hat[2, 2], n),
  Global_X6_Y3 = rep(Gamma_hat[2, 3], n),
  Global_X8_Y1 = rep(Gamma_hat[3, 1], n),
  Global_X8_Y2 = rep(Gamma_hat[3, 2], n),
  Global_X8_Y3 = rep(Gamma_hat[3, 3], n),
  
  Local_X7_Y1  = B_local_array[, "X3_Loc", "Y1"],
  Local_X7_Y2  = B_local_array[, "X3_Loc", "Y2"],
  Local_X7_Y3  = B_local_array[, "X3_Loc", "Y3"],
  
  Local_R2     = round(local_R2_sgwmr, 4),
  Local_MSE    = round(local_MSE_sgwmr, 4),
  Local_MAPE   = round(local_MAPE_sgwmr, 2)
)

cat("\n=== CUPLIKAN HASIL ESTIMASI DAN EVALUASI SGWMR ===\n")
print(sgwmr_results[, c("Provinsi", "Local_X3_Y1", "Local_X3_Y2", "Local_X3_Y3", "Local_R2", "Local_MSE", "Local_MAPE")])

# ================================================================
# UJI KESESUAIAN MODEL — LENGKAP (MMR vs GWMR vs SGWMR)
# ================================================================

# ----------------------------------------------------------------
# BAGIAN A: UJI KESESUAIAN MMR vs GWMR (F-test)
# ----------------------------------------------------------------
# H0: Model global (MMR) sudah cukup
# H1: GWMR signifikan lebih baik daripada MMR

df_MMR  <- (n * m) - (p_params * m)
df_GWMR <- (n * m) - (global_eta * m)

F_uji_mg <- ((RSS_MMR - RSS_GWMR) / (df_MMR - df_GWMR)) / (RSS_GWMR / df_GWMR)
p_val_mg <- pf(F_uji_mg, df1 = (df_MMR - df_GWMR), df2 = df_GWMR, lower.tail = FALSE)

cat("\n=== UJI KESESUAIAN MODEL: MMR vs GWMR ===\n")
cat("F-statistik :", round(F_uji_mg, 4), "\n")
cat("df1, df2    :", round(df_MMR - df_GWMR, 2), ",", round(df_GWMR, 2), "\n")
cat("P-value     :", round(p_val_mg, 5), "\n")
cat("Interpretasi: Jika p-value < 0.05, GWMR signifikan lebih baik dari MMR.\n\n")

# ----------------------------------------------------------------
# BAGIAN B: AICc UNTUK SGWMR
# ----------------------------------------------------------------
# 1. Trace bagian lokal (dari X_local saja)
tr_S_local <- 0
for (i in 1:n) {
  W_i <- adaptive_bisquare(dist_matrix[i, ], best_h)
  inv_XWX_L <- tryCatch(ginv(t(X_local) %*% W_i %*% X_local), error = function(e) NULL)
  if (!is.null(inv_XWX_L)) {
    H_i <- X_local[i, , drop = FALSE] %*% inv_XWX_L %*% t(X_local) %*% W_i
    tr_S_local <- tr_S_local + H_i[1, i]
  }
}

# 2. Parameter efektif SGWMR = trace bagian lokal + jumlah parameter global (OLS)
eta_SGWMR <- tr_S_local + n_glob_vars

# 3. RSS SGWMR
RSS_SGWMR <- SSE_global   # dari Langkah 6 SGWMR (SSE_global = sum((Y - Y_hat_sgwmr)^2))

# 4. AICc SGWMR (pakai residual multivariat penuh, bukan 1 respons)
E_SGWMR_full <- Y - Y_hat_sgwmr
Sigma_SGWMR <- (t(E_SGWMR_full) %*% E_SGWMR_full) / n
det_Sigma_SGWMR <- det(Sigma_SGWMR)

AICc_SGWMR <- (n * log(det_Sigma_SGWMR)) + (n * m * log(2 * pi)) + 
  (n * (n + eta_SGWMR) * m) / (n - eta_SGWMR - m - 1)

cat("Effective params SGWMR (eta):", round(eta_SGWMR, 4), "\n")
cat("RSS SGWMR                   :", round(RSS_SGWMR, 4), "\n")
cat("AICc SGWMR                  :", round(AICc_SGWMR, 4), "\n\n")

# ----------------------------------------------------------------
# BAGIAN C: UJI KESESUAIAN GWMR vs SGWMR (F-test)
# ----------------------------------------------------------------
# H0: SGWMR (lebih sederhana) sudah cukup
# H1: GWMR (penuh lokal) signifikan lebih baik daripada SGWMR

df_SGWMR <- (n * m) - (eta_SGWMR * m)

F_uji_gs <- ((RSS_SGWMR - RSS_GWMR) / (df_SGWMR - df_GWMR)) / (RSS_GWMR / df_GWMR)
p_val_gs <- pf(F_uji_gs, df1 = (df_SGWMR - df_GWMR), df2 = df_GWMR, lower.tail = FALSE)

cat("=== UJI KESESUAIAN MODEL: GWMR vs SGWMR ===\n")
cat("F-statistik :", round(F_uji_gs, 4), "\n")
cat("df1, df2    :", round(df_SGWMR - df_GWMR, 2), ",", round(df_GWMR, 2), "\n")
cat("P-value     :", round(p_val_gs, 5), "\n")

if (p_val_gs < 0.05) {
  cat("Interpretasi: P-value <", 0.05, "-> H0 ditolak.\n")
  cat("              GWMR signifikan lebih baik daripada SGWMR.\n")
  cat("              Variasi spasial pada prediktor yang dipaksa global (X1, X3, X4)\n")
  cat("              ternyata tetap penting -> model GWMR PENUH yang dipilih,\n")
  cat("              bukan SGWMR.\n\n")
} else {
  cat("Interpretasi: P-value >=", 0.05, "-> H0 gagal ditolak.\n")
  cat("              SGWMR sudah cukup dan lebih parsimoni dibanding GWMR penuh.\n")
  cat("              Model SGWMR yang dipilih.\n\n")
}

# ----------------------------------------------------------------
# BAGIAN D: TABEL PERBANDINGAN RINGKAS — MMR vs GWMR vs SGWMR
# ----------------------------------------------------------------
tabel_perbandingan_final <- data.frame(
  Model = c("MMR", "GWMR", "SGWMR"),
  RSS = c(round(RSS_MMR, 4), round(RSS_GWMR, 4), round(RSS_SGWMR, 4)),
  AICc = c(round(AICc_MMR, 4), round(best_AICc, 4), round(AICc_SGWMR, 4)),
  Eff_Params = c(p_params, round(global_eta, 4), round(eta_SGWMR, 4)),
  R2_Global = c("-", round(global_R2, 4), round(global_R2_sgwmr, 4))
)

cat("=== TABEL PERBANDINGAN MODEL (MMR vs GWMR vs SGWMR) ===\n")
print(tabel_perbandingan_final)