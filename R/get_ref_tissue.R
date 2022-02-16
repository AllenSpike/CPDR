#' @import octad.db
#' @importFrom stats IQR quantile
#' @importFrom dplyr select arrange mutate desc
#' @importFrom magrittr %>%
get_ref_tissue = function (case_id = NULL, adjacent = FALSE, source = "octad", cor_cutoff = "0",
                           n_varGenes = 500, method = "varGenes", expSet = NULL, control_size = length(case_id),
                           biopsy = NULL) {
  if (missing(case_id)) {
    stop("Case ids vector input not found")
  }
  if (source == "octad") {
    expSet = octad.db::EncoderDF
    case_id = case_id[case_id %in% colnames(expSet)]

    if (!is.null(biopsy)){
      phenoDF_f = filter(phenoDF, data.source=="GTEX")
      phenoDF_f = phenoDF_f[grep(biopsy,phenoDF_f$biopsy.site),]
      if (adjacent == TRUE) {
        adjacent_ids = as.vector(subset(phenoDF_f,
                                        sample.type == "adjacent")$sample.id)
        normal_id = as.vector(subset(phenoDF_f,
                                     sample.type == "normal")$sample.id)
        normal_id = c(adjacent_ids, normal_id)
      }
      else {
        normal_id = as.vector(subset(phenoDF_f,
                                     sample.type == "normal")$sample.id)
      }
    } else{
      if (adjacent == TRUE) {
        adjacent_ids = as.vector(subset(octad.db::phenoDF,
                                        sample.type == "adjacent")$sample.id)
        normal_id = as.vector(subset(octad.db::phenoDF,
                                     sample.type == "normal")$sample.id)
        normal_id = c(adjacent_ids, normal_id)
      }
      else {
        normal_id = as.vector(subset(octad.db::phenoDF,
                                     sample.type == "normal")$sample.id)
      }
    }

  }
  else if (source != "octad" & missing(expSet)) {
    stop("expSet is not supported")
  }
  else if (source != "octad") {
    normal_id = colnames(expSet)[!colnames(expSet) %in%
                                   case_id]
  }
  normal_id = normal_id[normal_id %in% colnames(expSet)]
  if (method == "random") {
    GTEXid <- sample(normal_id, size = control_size)
    return(GTEXid)
  }
  else if (method == "varGenes") {
    expSet_normal <- expSet[, as.vector(normal_id)]
    expSet_case <- expSet[, as.vector(case_id)]
    iqr_gene <- apply(expSet_normal, 1, stats::IQR)
    varying_genes <- order(iqr_gene, decreasing = TRUE)[1:min(n_varGenes,
                                                              length(iqr_gene))]
    normal_dz_cor <- cor(expSet_normal[varying_genes, ],
                         expSet_case[varying_genes, ], method = "spearman")
    normal_dz_cor_each <- apply(normal_dz_cor, 1, median)
    normal_dz_cor_eachDF = data.frame(cor = sort(normal_dz_cor_each,
                                                 decreasing = TRUE)) %>% mutate(sample.id = row.names(.)) %>%
      select(sample.id, cor)
    cutoff = stats::quantile(normal_dz_cor_eachDF$cor, probs = seq(0,
                                                                   1, 0.05), na.rm = TRUE)[paste0(cor_cutoff, "%")]
    GTEXid_temp <- subset(normal_dz_cor_eachDF, cor >= cutoff)
    GTEXid_temp = GTEXid_temp[order(GTEXid_temp$cor, decreasing = TRUE),
    ]
    GTEXid = GTEXid_temp$sample.id
    GTEXid <- GTEXid[seq_len(min(control_size, length(GTEXid)))]

    return(GTEXid)
  }
}
