library(tidyverse)


# dossier contenant les données
path_data <- "C:/Users/saout/OneDrive/Bureau/M2IA/P1_Non_Supervise/article/MicrobeDMMv1.0/Data"
# dossier résultats pgm C
path_res <- "C:/Users/saout/OneDrive/Bureau/M2IA/P1_Non_Supervise/article/MicrobeDMMv1.0"

# Fonction Q (découplée selon les K composantes) : à maximiser dans l'étape M

# première version codé 
Q_k <- function(alpha_k, ti_k, df ,nu ,eta) {
  
  res = 0
  beta_2 = sum( lgamma(alpha_k ))  - lgamma( sum(alpha_k ) )
  
  for(i in seq(1,N)){
    xi = as.numeric(df[i, ])
    beta_1 = sum( lgamma(xi + alpha_k) ) - lgamma(sum(xi + alpha_k))
    
    res = res + ti_k[i] *  ( beta_1 - beta_2 ) 
  }
  
  res = res - nu * sum(alpha_k) + eta * sum(log(alpha_k))
  
  return(res)
}


Q_k <- function(alpha_k, ti_k, X ,nu ,eta) {
  
  res = 0
  beta_2 = sum( lgamma(alpha_k ))  - lgamma( sum(alpha_k ) )
  X_alpha = sweep(X,2,alpha_k,FUN="+")
  beta_1 = rowSums(lgamma(X_alpha))-lgamma(rowSums(X_alpha))
  res_0 = sum(ti_k * (beta_1-beta_2) ) 
  res = res_0 - nu * sum(alpha_k) + eta * sum(log(alpha_k))
  
  return(res)
}


Q_k <- function(alpha_k, ti_k, X, n_i, nu, eta) {
  # alpha_k : vecteur de longueur S
  # ti_k    : vecteur de longueur N (responsabilités pour le cluster k)

  Nk <- sum(ti_k)
  A  <- sum(alpha_k)

  # terme log B(alpha_k)
  beta_2 <- sum(lgamma(alpha_k)) - lgamma(A)

  # lgamma(x_i + alpha_k) pour toute la matrice (N x S)
  # sweep ajoute alpha_k à chaque colonne
  G <- lgamma( sweep(X, 2, alpha_k, "+") )  # N x S
  # pour chaque i : sum_s lgamma(x_is + alpha_s)
  g1 <- rowSums(G)

  # pour chaque i : lgamma(sum_s x_is + sum_s alpha_s)
  g2 <- lgamma(n_i + A)

  beta1 <- g1 - g2

  res <- sum(ti_k * (beta1 - beta_2)) - nu * sum(alpha_k) + eta * sum(log(alpha_k))
  return(res)
}

# Gradient de la fonction Q (par rapport à un lambda_k)
g_Q_k <- function(alpha_k, ti_k, df, nu, eta) {
  S <- length(alpha_k)
  N <- nrow(df)
  res <- numeric(S)
  
  for (j in seq_len(S)) {
    alpha_j <- alpha_k[j]
    beta_2 <- digamma(alpha_j)       - digamma(sum(alpha_k))
    for (i in seq_len(N)) {
      xi  <- as.numeric(df[i, ])
      xij <- xi[j]
      beta_1 <- digamma(alpha_j + xij) - digamma(sum(xi + alpha_k))
      
      res[j] <- res[j] + ti_k[i] * beta_1 
    }
    Nk <- sum(ti_k)
    res[j] <- res[j] - Nk * beta_2 - nu + eta / alpha_j
  }
  res <- alpha_k * res # dérivée / log(alpha_k)
  return(res)
}

g_Q_k <- function(alpha_k, ti_k, X ,nu ,eta){
  
  S = length(alpha_k)
  res = rep(c(0),S)
  beta_2 = digamma(alpha_k) - digamma(sum(alpha_k))       # vecteur taille S
  X_alpha = sweep(X,2,alpha_k,FUN="+")
  beta_1 = digamma(X_alpha) - digamma(rowSums(X_alpha))   # matrice N*S
  tmp <- beta_1 - matrix(beta_2,
                         nrow = nrow(beta_1),
                         ncol = length(beta_2),
                         byrow = TRUE)
  
  res <- colSums( sweep(tmp, 1, ti_k, "*") ) - nu + eta/alpha_k
  res = alpha_k * res
  return(res)
}

g_Q_k<- function(alpha_k, ti_k, X, n_i, nu, eta) {
  Nk <- sum(ti_k)
  A  <- sum(alpha_k)

  # B = X + alpha_k colonne par colonne
  B <- sweep(X, 2, alpha_k, "+")   # N x S
  dig1 <- digamma(B)               # N x S

  # terme 1 : pour chaque s, sum_i t_i * digamma(x_is + alpha_s)
  # dig1 * ti_k : recyclage de ti_k sur les colonnes
  term1 <- colSums(dig1 * ti_k)

  # terme 2 : sum_i t_i * digamma(n_i + A), scalaire
  dig_nA <- digamma(n_i + A)
  term2_scalar <- sum(ti_k * dig_nA)

  # termes en alpha_k et A
  dig_alpha <- digamma(alpha_k)    # longueur S
  dig_A     <- digamma(A)          # scalaire

  # dQ / d alpha_k (vecteur de longueur S)
  dQ_dalpha <- term1 - term2_scalar - Nk * dig_alpha + Nk * dig_A - nu + eta / alpha_k

  # dQ / d lambda_k = alpha_k * dQ/dalpha_k (chain rule)
  grad_lambda <- alpha_k * dQ_dalpha
  return(grad_lambda)
}

#_______________________________________________________________________________

# Hyper Parametres

K = 4
nu=0.1 
eta = 0.1
set.seed(42)  

#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
# 1 . Lecture données

data <- read.csv(file.path(path_data, "Genera.csv"),header = TRUE,sep = ",")

df_mat <- data %>% column_to_rownames(var = colnames(data)[1])

# 2) Transposer
df_t <- t(df_mat)

# 3) Transformation en data.frame
df <- as.data.frame(df_t)
X  <- as.matrix(df)        # N x S
n_i <- rowSums(X)   

rm(data,df_mat,df_t)

N <- nrow(df)
S <- ncol(df)

#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
# Algo EM

# Initialisation K means
X_log <- log1p(X)
km <- kmeans(X_log, centers = K, nstart = 20)

init_clusters <- km$cluster
pi_k <- table(init_clusters) / N
print(pi_k)

alphas_list <- lapply(1:K, function(k){
  colMeans(df[init_clusters == k, , drop = FALSE]) + 0.1
})
print(pi_k)

# Boucle EM
n_limit = 1000
compt=0
conv = FALSE
Q_old <- -Inf
cst_pi = 11 # controle valeur Pi_k (pour empécher qu'un cluster disparaisse)
while(compt<=n_limit & conv == FALSE){
  
  pi_k_old = pi_k
  alphas_list_old = alphas_list
  
  cat("Iter", compt, "- début\n")
  
  #_____________________________________________________________________________
  # Etape E
  
  cat("  E-step...\n")
  
  # Calcul tik
  tik <- matrix(nrow = N, ncol = K) 
  log_tik_i <- numeric(K)
  for(i in seq(1,N)){
    for(j in seq(1,K)){
      alpha_k = alphas_list[[j]]
      xi = as.numeric(df[i, ])
      
      num = sum( lgamma(xi + alpha_k) ) - lgamma(sum(xi + alpha_k))
      den = sum( lgamma(alpha_k ))  - lgamma( sum(alpha_k ))
      
      
      # tik[i,j] = exp(log(pi_k[j]) + num-den) # cela peut exploser, du coup le log
      
      log_tik_i[j] <- log(pi_k[j]) + num - den
    }
    m <- max(log_tik_i)    
    w <- exp(log_tik_i - m)   # empeche les valeurs d'exploser
    tik[i, ] <- w / sum(w)    # exp(a-m)/somme_k(exp(a_k-m)) = exp(a)/somme_k(exp(a_k))
  }
  
  cat("  E-step OK\n")
  
  Nk <- colSums(tik)

  dead <- which(Nk < 0.05 * N)   # composantes mortes
  for(k in dead){
    # réinitialisation du cluster
    alphas_list[[k]] <- colMeans(df[sample(1:N, 10), ]) + 0.1
    pi_k[k] <- 1/K
  }
  pi_k <- pi_k / sum(pi_k)
 

  #_____________________________________________________________________________
  # Etape C
  est = apply(tik,1,which.max)
  
  
  #_____________________________________________________________________________
  # Etape M
  
  cat("  M-step (alphas)...\n")
  
  # Quasi Newton sur les alpha (lambda)
  for(j in seq(1,K)){
      alpha_k = alphas_list[[j]] 
      ti_k = tik[,j]
      alpha_k_new = optim(
            par = log(alpha_k),         # valeur initiale de la variable d’optimisation lambda
            
            fn = function(lambda){            # on optimise par rapport à lambda
              a = exp(lambda)                 # on passe exp(lambda)=alpha à la fonction, elle va maximiser par rapport à lambda -> les alpha_k seront  >0
              -Q_k(a, ti_k, X ,n_i, nu, eta)  # - Q  car optim() est un minimiseur
            } ,
            
            gr  = function(lambda) {          # on optimise par rapport lambda
              a <- exp(lambda)
              -g_Q_k(a, ti_k, X ,n_i, nu, eta)
            },
            
            method = "BFGS",
            control = list(maxit = 400, reltol = 1e-4)
      )
      alphas_list[[j]] = exp(alpha_k_new$par)
      
  }
  
  # mise à jour pi_k
  Nk <- colSums(tik)
  pi_k <- Nk + cst_pi
  pi_k <- pi_k / sum(pi_k)
  
  
  #_____________________________________________________________________________
  # Controle Convergence

  # Q somme sur K des Q_k
  Q_vals <- numeric(K)
  for (j in seq_len(K)) {
    Q_vals[j] <- Q_k(
      alpha_k = alphas_list[[j]],
      ti_k    = tik[, j],
      X       = X,
      n_i     = n_i,
      nu      = nu,
      eta     = eta
    )
  }
  Q_new <- sum(Q_vals)
  
  if (!is.infinite(Q_old)) {
    if (abs(Q_new - Q_old) < 1e-2) {
      conv <- TRUE
      print("les Q CV")
    }
  }
  Q_old <- Q_new
  
  if (max(abs(pi_k_old - pi_k)) < 1e-3) {
    conv <- TRUE
    print("les pi CV")
  }
  
  cat("iter", compt,
      "| Q:", round(Q_new, 2),
      "| pi_k:", paste(round(pi_k, 3), collapse = " / "), "\n")

  compt = compt+1
  
}

print(pi_k)


#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
# Variance cluster (somme des alpha_j_k pour chacun des K cluster)

theta_K <- sapply(alphas_list, sum)
theta_K <- round(theta_K)
print(theta_K)

# cluster 4 & 2 : variance inversement proportionnelle à 19 et 27 (+ forte variance)
# cluster 1 et 3 : variance inversement proportionnelle 48 et 54  (+ faible variance)


#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
# Matrice de confusion des individus par cluster (Entre les résultats du code C et le code ci dessus)

path_z <- file.path(path_res, "result.z")

# 1) Lire result.z (Holmes programme C)
L <- readLines(path_z)
L <- L[-1]                         # enlever la 1ère ligne "278,4"
parts <- strsplit(L, ",")
ids_C <- sapply(parts, function(x) x[1])                    # IDs échantillons côté C
probs_C <- t(sapply(parts, function(x) as.numeric(x[-1])))  # N x K
labels_C <- max.col(probs_C, ties.method = "first")       # clusters C (1..K)
tab_RC <- table(est, labels_C)              # tableau de fréaquence
tab_RC

#           programme C       programme R

# cluster       2                 1
#               1                 2
#               4                 3
#               3                 4


#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
# Compraraison des CLuster
Theta_R <- t(sapply(alphas_list, function(a) a / sum(a)))
colnames(Theta_R) <- colnames(df)

top5_R <- do.call(rbind, lapply(1:nrow(Theta_R), function(k) {
  idx <- order(Theta_R[k, ], decreasing = TRUE)[1:5]
  data.frame(
    cluster = paste0("R", k),
    genus   = colnames(df)[idx],
    theta   = Theta_R[k, idx],
    row.names = NULL
  )
}))
top5_R


mix <- read.csv(  file.path(path_res, "result.mixture"),header = FALSE,stringsAsFactors = FALSE)
meta<- read.table(file.path(path_data, "TwinStudy.t"), header = FALSE)
K <- as.integer(mix[1,1])  
S <- as.integer(mix[1,2])

alpha_C <- as.matrix(mix[3:(2+S), 2:(K+1)])
Theta_C <- t(apply(alpha_C, 2, function(a) a / sum(a)))
colnames(Theta_C) <- mix[3:(2+S),1]

top5_C <- do.call(rbind, lapply(1:nrow(Theta_C), function(k) {
  idx <- order(Theta_C[k, ], decreasing = TRUE)[1:5]
  data.frame(
    cluster = paste0("C", k),
    genus   = colnames(df)[idx],
    theta   = Theta_C[k, idx],
    row.names = NULL
  )
}))
top5_C

# mapping des clusters R -> C (d’après la matrice de confusion)
map_R_to_C <- c(`1` = 2,  
                `2` = 1,  
                `3` = 4,  
                `4` = 3)  

# listes par cluster
top5_R_list <- split(top5_R, top5_R$cluster)
top5_C_list <- split(top5_C, top5_C$cluster)

# Fonction de comparaison
compare_clusters <- function(k) {
  c_target <- map_R_to_C[as.character(k)]
  list(
    R = top5_R_list[[paste0("R", k)]],
    C = top5_C_list[[paste0("C", c_target)]]
  )
}

# Comparaison cluster
compare_top5 <- lapply(1:4, compare_clusters)

compare_top5

#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
#_______________________________________________________________________________
# tableau de fréauence IMC / cluster R

# on récupère les IMC des individus

table(meta$V1)
cluster_R <- max.col(tik) 

tab <- prop.table(table(meta$V1, cluster_R), 1)
rownames(tab) <- c("Lean", "Obese", "Overweight")
round(100 * tab, 1)

# On a retrouvé les clusters








