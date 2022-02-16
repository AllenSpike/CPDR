#' Compute sRGES
#'
#' @param dz_signature disease signature. Make sure input data frame has a gene Symbol column, otherwise an error is produced. It must be an UPPERCASE gene symbol.
#' @param max_gene_size maximum number of disease genes used for drug prediction. By default 50 for each side (up/down).
#' @param LINCS_data LINCS data. By default set to NULL, which means use LINCS_978.
#' @param permutations number of permutations, by default 10000.
#'
#' @return A data.frame containing scores and p.values for every instance. data.frame contains drug id in pert_iname collumn, n contains the number of instances for this drug, mean, median and sd of sRGES RGES sores.
#' @export

get_reverse_score = function(dz_signature = NULL, max_gene_size = 500,
                             LINCS_data = NULL, permutations = 10000) {

  getsRGES <- function(RGES, cor, pert_dose, pert_time, diff,
                       max_cor) {
    sRGES <- RGES
    pert_time <- ifelse(as.numeric(pert_time) < 24, "short", "long")
    pert_dose <- ifelse(pert_dose < 10, "low", "high")
    if (pert_time == "short" & pert_dose == "low") {
      sRGES <- sRGES + diff[4]
    }
    if (pert_dose == "low" & pert_time == "long") {
      sRGES <- sRGES + diff[2]
    }
    if (pert_dose == "high" & pert_time == "short") {
      sRGES <- sRGES + diff[1]
    }
    return(sRGES*cor/max_cor)
  }

  cmap_score_ultimate = function(sig_up, sig_down, drug_signature) {
    num_genes <- length(drug_signature)
    ks_up <- 0
    ks_down <- 0
    connectivity_score <- 0

    drug_signature = rank(drug_signature)
    up_tags_rank = drug_signature[as.vector(sig_up)]
    down_tags_rank = drug_signature[as.vector(sig_down)]
    up_tags_position = sort(up_tags_rank)
    down_tags_position = sort(down_tags_rank)
    num_tags_up <- length(up_tags_position)
    num_tags_down <- length(down_tags_position)

    if (num_tags_up > 1) {
      a_up <- 0
      b_up <- 0
      a_up <- max(sapply(seq_len(num_tags_up), function(j) {
        j/num_tags_up - up_tags_position[j]/num_genes
      }))
      b_up <- max(sapply(seq_len(num_tags_up), function(j) {
        up_tags_position[j]/num_genes - (j - 1)/num_tags_up
      }))
      if (a_up > b_up) {
        ks_up <- a_up
      }
      else {
        ks_up <- -b_up
      }
    }
    else {
      ks_up <- 0
    }
    if (num_tags_down > 1) {
      a_down <- 0
      b_down <- 0
      a_down <- max(sapply(seq_len(num_tags_down), function(j) {
        j/num_tags_down - down_tags_position[j]/num_genes
      }))
      b_down <- max(sapply(seq_len(num_tags_down), function(j) {
        down_tags_position[j]/num_genes - (j - 1)/num_tags_down
      }))
      if (a_down > b_down) {
        ks_down <- a_down
      }
      else {
        ks_down <- -b_down
      }
    }
    else {
      ks_down <- 0
    }
    if (ks_up == 0 & ks_down != 0) {
      connectivity_score <- -ks_down
    }
    else if (ks_up != 0 & ks_down == 0) {
      connectivity_score <- ks_up
    }
    else if (sum(sign(c(ks_down, ks_up))) == 0) {
      connectivity_score <- ks_up - ks_down
    }
    else {
      connectivity_score <- ks_up - ks_down
    }
    return(connectivity_score)
  }


  if (missing(dz_signature)) {
    stop("Disease signature input not found")
  }
  if (is.null(dz_signature$Symbol) | is.null(dz_signature$log2FoldChange)) {
    stop("Either Symbol or log2FoldChange collumn in Disease signature is missing")
  }

  if (is.null(perturbation)){
    lincs_signatures = octad.db::lincs_signatures
    lincs_sig_info = octad.db::lincs_sig_info
  } else {
    lincs_signatures = perturbation$matrix
    row.names(lincs_signatures) = perturbation$row_id[["pr_gene_symbol"]]
    lincs_sig_info = perturbation$col_id
  }

  lincs_sig_info <- subset(lincs_sig_info, id %in% colnames(lincs_signatures))
  lincs_sig_info <- lincs_sig_info[!duplicated(lincs_sig_info$id),]

  sig.ids <- lincs_sig_info$id
  gene.list <- toupper(rownames(lincs_signatures))
  dz_signature <- subset(dz_signature, Symbol %in% gene.list)
  dz_genes_up <- subset(dz_signature, log2FoldChange > 0)
  dz_genes_up = dz_genes_up[order(dz_genes_up$log2FoldChange,
                                  decreasing = TRUE), ]
  dz_genes_down <- subset(dz_signature, log2FoldChange < 0)
  dz_genes_down = dz_genes_down[order(dz_genes_down$log2FoldChange,
                                      decreasing = TRUE), ]
  if (nrow(dz_genes_up) > max_gene_size) {
    dz_genes_up <- dz_genes_up %>% head(max_gene_size)
  }
  if (nrow(dz_genes_down) > max_gene_size) {
    dz_genes_down <- dz_genes_down %>% head(max_gene_size)
  }

  dz_cmap_scores <- NULL
  cmap_exp_sig <- Rfast::colRanks(-1 * lincs_signatures,
                                  method = "max")
  names.list <- list(rownames(lincs_signatures), colnames(lincs_signatures))
  dimnames(cmap_exp_sig) <- names.list
  cat(paste("Started sRGES computation. Average computation time ~1-3mins."),
      "\n")
  start_time = Sys.time()
  pb <- txtProgressBar(min = 1, max = permutations, style = 3)
  i = 0
  time_vector = 0
  cmap_exp_signature = data.frame(ids = gene.list, rank = cmap_exp_sig[,
                                                                       as.character(sig.ids)])
  cmap_exp_sig = as.data.frame(cmap_exp_sig)
  cmap_exp_sig$ids = NULL

  dz_cmap_scores = apply(cmap_exp_sig[as.character(sig.ids)],
                         2, FUN = function(x) cmap_score_ultimate(dz_genes_up$Symbol,
                                                                  dz_genes_down$Symbol, drug_signature = x))

  random_sig_ids <- sample(colnames(lincs_signatures),
                           permutations, replace = TRUE)
  random_cmap_scores <- NULL
  cmap_exp_signature = as.data.frame(Rfast::colRanks(-1 *
                                                       lincs_signatures[, as.character(random_sig_ids)],
                                                     method = "max"))
  random_cmap_scores = apply(cmap_exp_signature, 2, FUN = function(x) cmap_score_ultimate(sample(1:length(dz_genes_up$Symbol),
                                                                                                 replace = TRUE), sample(1:length(dz_genes_down$Symbol),
                                                                                                                         replace = TRUE), drug_signature = x))
  p <- sapply(dz_cmap_scores, function(score) {
    sum(random_cmap_scores < score)/length(random_cmap_scores)
  })
  padj <- p.adjust(p, "fdr")
  results <- data.frame(id = sig.ids, cmap_score = dz_cmap_scores,
                        p, padj)
  results <- merge(results, lincs_sig_info, by = "id")
  results <- results[order(results$cmap_score), ]

  lincs_drug_prediction <- results

  lincs_drug_prediction_subset <- subset(lincs_drug_prediction,
                                         pert_dose > 0 & pert_time %in% c(6, 24))
  lincs_drug_prediction_pairs <- merge(lincs_drug_prediction_subset,
                                       lincs_drug_prediction_subset, by = c("pert_iname",
                                                                            "cell_id"))
  lincs_drug_prediction_pairs <- subset(lincs_drug_prediction_pairs,
                                        id.x != id.y & pert_time.x == 24 & pert_dose.x ==
                                          10)


  lincs_drug_prediction_pairs$cmap_diff <- lincs_drug_prediction_pairs$cmap_score.x -
    lincs_drug_prediction_pairs$cmap_score.y
  lincs_drug_prediction_pairs$dose <- round(log(lincs_drug_prediction_pairs$pert_dose.y,
                                                2), 1)
  lincs_drug_prediction_pairs$dose_bin <- ifelse(lincs_drug_prediction_pairs$pert_dose.y <
                                                   10, "low", "high")

  diff <- tapply(lincs_drug_prediction_pairs$cmap_diff,
                 paste(lincs_drug_prediction_pairs$dose_bin, lincs_drug_prediction_pairs$pert_time.y),
                 mean)
  pred <- lincs_drug_prediction
  pred$medcor <- 1
  pred$RGES <- sapply(1:nrow(pred), function(id) {
    getsRGES(pred$cmap_score[id], pred$medcor[id], pred$pert_dose[id],
             pred$pert_time[id], diff, max(pred$medcor))
  })
  cmpd_freq <- table(pred$pert_iname)
  pred <- subset(pred, pert_iname %in% names(cmpd_freq[cmpd_freq >
                                                         0]))

  pred_merged <- pred %>% group_by(pert_iname) %>% dplyr::summarise(mean = mean(RGES), min = min(RGES), max = max(RGES),
                                                                    n = length(RGES), median = median(RGES), sd = sd(RGES))
  pred_merged$sRGES <- pred_merged$mean
  pred_merged <- pred_merged[order(pred_merged$sRGES),
  ]
  gc()
  return(pred_merged)
}
