#' Differential expression analysis
#'
#' @param data The purified profiles of subgoups and control groups for each clinical query sample
#' @param DE_method edgeR, DESeq2 or limma DE analysis.
#' @param normalize_samples if TRUE, RUVSeq normalization is applied to either EdgeR or DESeq. No normalization needed for limma+voom.
#' @param threshold_log2foldchange threshold for log2foldchange. 2.5 by default.
#' @param threshold_pval threshold for pvalue. 0.05 by default.
#' @param threshold_adjpval threshold for adjusted pvalue. 0.001 by default.
#' @param k either k=1 (by default), k=2 or k=3, number of factors used in model matrix construction in RUVSeq normalization if normalize_samples=TRUE.
#' @param n_topGenes number of empiricaly differentially expressed genes estimated for RUVSeq normalization. Default is 5000.
#' @param parallel_cores number of cores to be used for parallel computing in DESeq2.
#'
#' @return A list with list of differentially expressed genes.
#' @export
#' @importFrom octad diffExp
#' @import dplyr

get_diff = function(data = NULL,DE_method = "edgeR",normalize_samples = TRUE,
                    threshold_log2foldchange = 2.5, threshold_pval = 0.05, threshold_adjpval = 0.001,
                    k = 1, n_topGenes = 500, parallel_cores = 2) {

  de = DE_method
  ns = normalize_samples
  nt = n_topGenes
  pc = parallel_cores

  diff_analysis = function(i,d,th1,th2,th3,de,ns,nt,pc,k){
    case_id <- colnames(d[[i]]$case)
    control_id <- colnames(d[[i]]$control)
    gene <- intersect(row.names(d[[i]]$case_puify),row.names(d[[i]]$control))
    d = cbind(log2(d[[i]]$case_puify[gene,] + 1),log2(d[[i]]$control[gene,] + 1))
    res_group <- octad::diffExp(case_id,control_id,source='side',expSet=d,
                                DE_method = de, k = k, n_topGenes = nt, normalize_samples = ns,
                                output=F, annotate=F, parallel_cores = pc)
    merged_gene_info = octad.db::merged_gene_info
    merged_gene_info$ensembl <- as.vector(merged_gene_info$ensembl)
    merged_gene_info$V1 = NULL
    merged_gene_info$Symbol = merged_gene_info$gene
    merged_gene_info$gene = NULL
    res_group <- left_join(res_group, merged_gene_info, by = c(identifier = "Symbol_autho"))
    res_group <- arrange(res_group, log2FoldChange)

    res_group <- subset(res_group,abs(log2FoldChange) > th1&pvalue < th2&padj < th3)
    return(res_group)
  }

  res <- lapply(seq(length(data)),function(i,d,th1,th2,th3,de,ns,nt,pc,k){diff_analysis(i,d,th1,th2,th3,de,ns,nt,pc,k)},
                d=data,th1=threshold_log2foldchange,th2=threshold_pval,th3=threshold_adjpval,
                de=de,ns=ns,nt=nt,pc=pc,k=k)
  names(res) <- names(data)
  return(res)
}

