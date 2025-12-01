# Dirichlet-Multinomial Mixtures for Microbial Metagenomics

Ce projet propose une reproduction et une analyse des résultats de l'article de **Holmes et al. (2012)** concernant l'utilisation de modèles de mélanges Dirichlet-Multinomial (DMM) pour l'analyse de données de métagénomique microbienne.

Il contient une implémentation "from scratch" de l'algorithme EM et une implémentation utilisant le package R dédié qui reproduit les figures clés de l'article.

## 👥 Auteurs
* **Rémi Saout** : Implémentation de l'algorithme EM "from scratch" et reproduction des résultats en langage C.
* **Grégoire Pelletier** : Reproduction des résultats via le package R `DirichletMultinomial` et analyse comparative.

## 📂 Structure du projet

* **`EM.R`** : Implémentation manuelle de l'algorithme Espérance-Maximisation (EM) pour les mélanges Dirichlet-Multinomial. Inclut le calcul de la fonction Q, des gradients, et la comparaison avec les résultats du code C.
* **`Reproduction_resultat.R`** : Script utilisant le package `DirichletMultinomial` pour :
  * La sélection de modèle (Critères Laplace vs AIC/BIC).
  * L'analyse de stabilité des clusters (Principe d'Anna Karénine).
  * La classification supervisée (Lean vs Obese) et courbes ROC.
* **`Résultat C/`** : Contient les sorties brutes du programme C (fichiers `.fit`, `.mixture`, `.z`) utilisées pour la validation croisée.
* **`figures/`** : Graphiques générés (Sélection de modèle, Homogénéité des clusters, ROC).

## 🛠️ Prérequis

Le projet nécessite **R** et les librairies suivantes :

```r
install.packages(c("tidyverse", "DirichletMultinomial", "lattice", "parallel", "pROC"))
