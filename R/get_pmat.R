#' @importFrom magrittr %>%
#' @importFrom foreach foreach %do%
#'
get_pmat = function(pmat_counts,clinical_profile_set,topk=1500){

  if (is.vector(clinical_profile_set)) {
    pmat <- pmat_counts[intersect(row.names(pmat_counts),
                                  names(clinical_profile_set)),]
    pmat.median = apply(pmat,
                        1, median)

    pmat.expressed.gene <- names(pmat.median)[pmat.median >
                                                1]
    tmp <- pmat[pmat.expressed.gene,
    ]
    tmp.rank <- apply(tmp, 2, rank)
    rank.mean <- apply(tmp.rank, 1, mean)
    rank.sd <- apply(tmp.rank, 1, sd)
    marker.gene <- names(sort(rank.sd, decreasing = TRUE))[1:topk]

    marker.gene <- intersect(names(clinical_profile_set),
                             (marker.gene))
    marker.gene <- intersect(rownames(pmat_counts),
                             (marker.gene))
    correlation.matrix <- cor(clinical_profile_set[marker.gene],
                              pmat_counts[marker.gene, ], method = "spearman")
  } else {
    pmat <- pmat_counts[intersect(row.names(pmat_counts),
                                  row.names(clinical_profile_set)),]
    pmat.median = apply(pmat,
                        1, median)

    pmat.expressed.gene <- names(pmat.median)[pmat.median >
                                                1]
    tmp <- pmat[pmat.expressed.gene,
    ]
    tmp.rank <- apply(tmp, 2, rank)
    rank.mean <- apply(tmp.rank, 1, mean)
    rank.sd <- apply(tmp.rank, 1, sd)
    marker.gene <- names(sort(rank.sd, decreasing = TRUE))[1:topk]

    marker.gene <- intersect(rownames(clinical_profile_set),
                             (marker.gene))
    marker.gene <- intersect(rownames(pmat_counts),
                             (marker.gene))
    correlation.matrix <- cor(clinical_profile_set[marker.gene,
    ], pmat_counts[marker.gene, ], method = "spearman")
  }

  correlation.matrix[is.na(correlation.matrix)] = 0
  median.cor <- apply(correlation.matrix, 2,
                      median) %>% sort(decreasing = TRUE)

  best.cell.line <- names(median.cor)[1]
  p.value.vec <- foreach::foreach(cell.line = setdiff(names(median.cor),
                                                      best.cell.line), .combine = "c") %do% {
                                                        v <- correlation.matrix[, cell.line]
                                                        p.value <- wilcox.test(correlation.matrix[, best.cell.line],
                                                                               v, alternative = "greater", paired = TRUE)$p.value
                                                      }

  names(p.value.vec) <- setdiff(names(median.cor),
                                best.cell.line)
  fdr.vec <- p.adjust(p.value.vec, method = "fdr")
  return(list(median.cor = median.cor, best.sample = best.cell.line, compare.p.value.vec = p.value.vec,
              compare.fdr.vec = fdr.vec, correlation.matrix = correlation.matrix))
}
