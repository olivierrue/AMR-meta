#!/usr/bin/env Rscript

library(Matrix)
library(stringr)
library(glmnet)

args = commandArgs(trailingOnly=TRUE)

progdir = args[1]   # base directory
in_file = args[2]   # input feature file
outdir = args[3]    # output directory

# === Lire les features avec noms de reads ===
all_lines <- readLines(in_file)
all_lines <- all_lines[!grepl("^END$", all_lines)]  # Supprime la ligne END éventuelle

read_ids <- character()
rows <- list()

for (line in all_lines) {
  parts <- str_split(line, ",")[[1]]
  read_ids <- c(read_ids, parts[1])
  if (length(parts) > 1) {
    indices <- as.integer(parts[-1])
    rows[[length(read_ids)]] <- indices
  }
}

# === Convertir en matrice creuse binaire ===
tot_samples <- length(read_ids)
non_empty <- which(sapply(rows, length) > 0)
ij <- do.call(rbind, lapply(non_empty, function(i) {
  cbind(rep(i, length(rows[[i]])), rows[[i]])
}))

kmer_data <- sparseMatrix(i = ij[,1], j = ij[,2], x = TRUE)

# Compléter les colonnes s'il en manque
if (ncol(kmer_data) < 138260) {
  add_zeros <- Matrix(FALSE, nrow = nrow(kmer_data), ncol = 138260 - ncol(kmer_data), sparse = TRUE)
  kmer_data <- cbind(kmer_data, add_zeros)
}

# Ajouter lignes vides si certains reads n'ont aucun kmer
if (nrow(kmer_data) < tot_samples) {
  add_zeros <- Matrix(FALSE, ncol = 138260, nrow = tot_samples - nrow(kmer_data), sparse = TRUE)
  kmer_data <- rbind(kmer_data, add_zeros)
}

# === Calcul des métacaractéristiques ===
metaf_matrix <- readRDS(paste0(progdir, '/data/metaf.matrix.rds'))
metaf_data <- kmer_data %*% metaf_matrix

# === Définir les classes à prédire ===
classes <- c("Aminoglycosides", "betalactams", "Drug_and_biocide_resistance",
             "Fluoroquinolones", "Glycopeptides", "Lipopeptides", "MLS",
             "Multi-biocide_resistance", "Multi-drug_resistance", "Multi-metal_resistance",
             "Phenicol", "Sulfonamides", "Tetracyclines")

# === Prédictions ===
predictions_kmer <- predictions_metaf <- NULL
for (i in 1:length(classes)) {
  cv.fit.lasso <- readRDS(paste0(progdir, '/models/', classes[i], '.lasso.kmer.model.rds'))
  cv.fit.ridge <- readRDS(paste0(progdir, '/models/', classes[i], '.ridge.metaf.model.rds'))

  predictions_kmer <- cbind(predictions_kmer, predict(cv.fit.lasso, newx = kmer_data, type="response"))
  predictions_metaf <- cbind(predictions_metaf, predict(cv.fit.ridge, newx = metaf_data, type="response"))
}
colnames(predictions_kmer) <- classes
colnames(predictions_metaf) <- classes

# === Écriture des fichiers de sortie ===
in_file_name <- str_split(basename(in_file), "_")[[1]][3]  # extraire l'identifiant
results_kmer <- cbind(ReadID = read_ids, format(predictions_kmer, digits=4))
results_metaf <- cbind(ReadID = read_ids, format(predictions_metaf, digits=4))

write.table(results_kmer, file = paste0(outdir, '/kmer_predictions_', in_file_name),
            sep = ',', quote = FALSE, row.names = FALSE)

write.table(results_metaf, file = paste0(outdir, '/metaf_predictions_', in_file_name),
            sep = ',', quote = FALSE, row.names = FALSE)
