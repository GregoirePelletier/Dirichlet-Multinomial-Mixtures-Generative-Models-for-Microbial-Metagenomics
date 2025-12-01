# Script de reproduction des résultats
# Auteurs : Saout, Pelletier

library(DirichletMultinomial)
library(lattice)
library(parallel)

# --- 1. CHARGEMENT ET PREPARATION ---
fl <- system.file(package="DirichletMultinomial", "extdata", "Twins.csv")
counts <- t(as.matrix(read.csv(fl, row.names=1)))

cat("Dimensions initiales : ", dim(counts), "\n")

# --- 2. ENTRAINEMENT DES MODÈLES ---
k_range <- 1:7
set.seed(123)
fits <- lapply(k_range, function(k) {
  cat("Fitting K =", k, "...\n")
  dmn(counts, k = k, verbose = FALSE)
})

# --- 3. SELECTION DE MODELE ---
laplace_scores <- sapply(fits, laplace)
bic_scores <- sapply(fits, BIC)
aic_scores <- sapply(fits, AIC)

best_idx <- which.min(laplace_scores)
best_k <- k_range[best_idx]
best_model <- fits[[best_idx]]

# Visualisation Sélection
cat("\n--- RÉSULTATS DE SÉLECTION ---\n")
cat("Meilleur K selon Laplace :", best_k, "\n")

df_scores <- data.frame(Laplace = laplace_scores, AIC = aic_scores, BIC = bic_scores)
df_plot_scaled <- scale(df_scores)

matplot(k_range, df_plot_scaled, type="b", pch=19, lwd=2,
        ylab="Score normalisé (Plus bas est mieux)", 
        xlab="Nombre de clusters (K)",
        main="Sélection de modèle : Laplace vs AIC/BIC")
legend("topleft", legend=colnames(df_scores), col=1:3, pch=19, lty=1, lwd=2)
grid()

# --- 4. ANALYSE DU PRINCIPE D'ANNA KARENINE ---

params <- mixturewt(best_model)
thetas <- params[, "theta"]
weights <- params[, "pi"]

# Graphique : Stabilité des clusters
names(thetas) <- paste0("C", 1:best_k, "\n(w=", round(weights, 2), ")")

cols <- ifelse(thetas > mean(thetas), "skyblue", "orange")

barplot(thetas, col=cols, 
        main="Principe d'Anna Karénine (Stabilité)",
        ylab="Précision Theta (Homogénéité)", las=1)
abline(h=mean(thetas), col="red", lty=2)
legend("topright", legend=c("Stable (>Moy)", "Variable (<Moy)"), 
       fill=c("skyblue", "orange"), cex=0.8)

# --- 5. CLASSIFICATION SUPERVISÉE (Lean vs Obese) ---

cat("\n--- Classification Supervisée (Lean vs Obese) ---\n")

f1 <- system.file(package="DirichletMultinomial", "extdata", "TwinStudy.t")
pheno0 <- scan(f1, quiet=TRUE)
pheno <- factor(c("Lean", "Obese", "Overwt")[pheno0 + 1])
names(pheno) <- rownames(counts)

keep <- pheno %in% c("Lean", "Obese")
counts_sub <- counts[keep, ]
pheno_sub <- droplevels(pheno[keep])

# Entraînement du Classifieur "Groupé"
targets_k <- c(Lean = 1, Obese = 3)
bestgrp <- dmngroup(counts_sub, group = pheno_sub, k = targets_k, verbose = FALSE)

# Prédictions (Classes et Probabilités)
predictions_class <- predict(bestgrp, counts_sub, assign = TRUE)
predictions_prob  <- predict(bestgrp, counts_sub, assign = FALSE) # Pour le ROC

cat("\n--- Matrice de Confusion ---\n")
conf_mat <- table(Actual = pheno_sub, Predicted = predictions_class)
print(conf_mat)

# Comparaison des Modèles (Laplace)

laplace_global <- laplace(best_model) 
laplace_grouped <- sum(sapply(bestgrp, laplace))

cat("\n--- Comparaison Laplace ---\n")
cat("Score Global (Single) :", round(laplace_global, 0), "\n")
cat("Score Groupé (Lean+Obese) :", round(laplace_grouped, 0), "\n")
cat("Gain :", round(laplace_global - laplace_grouped, 0), "\n")

# --- 6. COURBES ROC ---
#install.packages("pROC")
library(pROC)

# 1. Extraction des probabilités pour la classe "Obese"
# predict(..., assign=FALSE) retourne une matrice avec les probas pour chaque classe
prob_obese <- predictions_prob[, "Obese"]

# levels : on précise l'ordre (Négatif, Positif) -> Lean est le contrôle, Obese le cas
# direction : "auto" laisse R deviner, mais "<" assure que si Proba augmente, on va vers Obese
roc_obj <- roc(response = pheno_sub, 
               predictor = prob_obese, 
               levels = c("Lean", "Obese"))

# 3. Affichage
par(mfrow=c(1,1))
plot(roc_obj, 
     col="blue", lwd=3, legacy.axes=TRUE, # legacy.axes met 1-Specificité en X (standard)
     main="Fig 5: Performance (Training Set)",
     xlab="Taux Faux Positifs (1 - Specificité)",
     ylab="Taux Vrais Positifs (Sensibilité)",
     print.auc = TRUE, # Affiche l'AUC directement sur le graphe
     print.auc.y = 0.4, 
     print.auc.col = "blue")

grid()
abline(a=0, b=1, lty=2, col="gray")

cat("AUC calculé : ", auc(roc_obj), "\n")
