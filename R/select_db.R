#'  Preprocess transcriptome profiling and filter TCGA samples
#'
#' @param Assay The TCGA data set of selected cancer type.
#' @param cmat The RNA-seq count expression profiles of clinical query samples.
#' @param MSI_status Filter TCGA samples with the MSI status (MSI or MSS).
#' @param driver_gene  A character with the driver gene symbol.
#' @param MUT_status Filter TCGA samples with the MUT status of driver gene (MUT or NULL).
#' @param CNA_status Filter TCGA samples with the CNA status of driver gene (AMP, GAIN, HETLOSS, HOMDEL).
#' @param geneList A data.frame of 2 column with gene and status. Only used in multiple driver genes.
#' @param minSampleSize Minimal size of the outpputted samples (10 by default).
#' @param removeBatchEffect Remove batch effects from expression data using limma (TRUE by default).
#' @param OrgDb The human annotation db for ID convert.
#'
#' @return A large list with processed profiles of clinical query samples (cmat), and processed and filtered profiles of selected TCGA data set (mat).
#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom clusterProfiler bitr
#' @importFrom limma removeBatchEffect
#' @import dplyr
#' @import octad.db

select_db = function(Assay, cmat, MSI_status = NULL, driver_gene = NULL,
                     MUT_status = NULL, CNA_status = NULL, geneList = NULL, minSampleSize = 10,
                     removeBatchEffect = TRUE, OrgDb){

  sample_s = NULL
  sampleID = octad.db::phenoDF$sample.id
  clinical <- as.data.frame(Assay@colData@listData)
  mat <- SummarizedExperiment::assay(Assay,"mrna_seq_v2_rsem")
  mut <- SummarizedExperiment::assay(Assay,"mutations")
  cna <- SummarizedExperiment::assay(Assay,"cna")

  row.names(mut) <- apply(mut, 1, function(x){x[which(!is.na(x))[1]]})
  cna <- ifelse(cna == 2, 'AMP',
                ifelse(cna == 1, 'GAIN',
                       ifelse(cna == -1, 'HETLOSS',
                              ifelse(cna == -2, 'HOMDEL', NA))))

  if (sum(is.na(suppressWarnings(as.numeric(row.names(cmat))))) == 0 | sum(grepl('^ENSG',row.names(cmat))) > 0) {
    stop('Row names of cmat should be SYMBOL.')
  }
  if (!is.null(MSI_status)){
    if (MSI_status == 'MSI') {
      sample_s <- clinical$PATIENT_ID[clinical$MSI_SCORE_MANTIS > 0.6]
    }
    if (MSI_status == 'MSS') {
      sample_s <- clinical$PATIENT_ID[clinical$MSI_SCORE_MANTIS < 0.4]
    }
    if (length(sample_s) < minSampleSize) {
      if (length(sample_s) == 0) {
        stop(paste0('No patients with ',MSI_status))
      }
      warning(paste0('Insufficient patients with ',MSI_status))
    }
  }

  if (!is.null(driver_gene)){
    if (!is.null(MUT_status)){
      sample_s <- colnames(mut)[which(!is.na(mut[driver_gene,]))]
      if (length(sample_s) < minSampleSize) {
        if (is.na(sample_s)) {
          stop(paste0('No patients with ',driver_gene,' ',MUT_status))
        }
        warning(paste0('Insufficient patients with ',driver_gene,' ',MUT_status))
      }
    }
    if (!is.null(CNA_status)){
      sample_s <- ifelse(CNA_status == 'AMP', colnames(cna)[which(cna[driver_gene,] == 'AMP')],
                         ifelse(CNA_status == 'GAIN', colnames(cna)[which(cna[driver_gene,] == 'GAIN')],
                                ifelse(CNA_status == 'HETLOSS', colnames(cna)[which(cna[driver_gene,] == 'HETLOSS')],
                                       ifelse(CNA_status == 'HOMDEL', colnames(cna)[which(cna[driver_gene,] == 'HOMDEL')]))))
      if (length(sample_s) < minSampleSize) {
        if (is.na(sample_s)) {
          stop(paste0('No patients with ',driver_gene,' ',CNA_status))
        }
        warning(paste0('Insufficient patients with ',driver_gene,' ',CNA_status))
      }
    }
  }

  if (!is.null(geneList)) {
    if (!colnames(geneList) %in% c('gene','status')) {
      stop('GeneList must be a dataframe with `gene` and `status` columns. ')
    }

    mut_list <- geneList[,geneList$status == 'MUT']
    cna_list <- list(
      geneList[,geneList$status == 'AMP'],
      geneList[,geneList$status == 'GAIN'],
      geneList[,geneList$status == 'HETLOSS'],
      geneList[,geneList$status == 'HOMDEL']
    )

    select_sample = function(ls){
      return(ifelse(dim(ls)[2] == 1,
                    colnames(mut)[which(!is.na(mut[ls$gene,]))],
                    colnames(mut)[which(apply(mut[ls$gene,], 2, function(x){sum(is.na(x))}) == 0)]))
    }

    sample_s <- Reduce(intersect,
                       list(select_sample(mut_list),
                            unlist(lapply(cna_list, function(x){select_sample(x)}))))

    if (length(sample_s) < minSampleSize) {
      warning(paste0('Insufficient patients with genetic alterations in `geneList`.'))
    }
  }


  gene <- suppressMessages(suppressWarnings(clusterProfiler::bitr(row.names(mat), fromType='ENTREZID', toType='SYMBOL',OrgDb)))
  gene <- gene[!duplicated(gene$ENTREZID),]
  mat <- mat[row.names(mat) %in% gene$ENTREZID,]
  mat <- mat[gene$ENTREZID,]
  row.names(mat) <- gene$SYMBOL

  if (is.null(sample_s)) {sample_s <- colnames(mat)} else {sample_s <- intersect(paste0(sample_s,'-01'),colnames(mat))}
  sample_s <- intersect(sampleID,sample_s)
  mat <- mat[,sample_s]

  commongene <- intersect(row.names(mat)[complete.cases(mat)],row.names(cmat)[complete.cases(cmat)])
  cmat <- cmat[commongene,]
  mat <- mat[commongene,]

  pmat <- list(
    mat = mat,
    cmat = cmat)

  if (isTRUE(removeBatchEffect)) {
    message('... Removing batch effect ...')
    mat <- log2(mat + 1)
    cmat <- log2(cmat + 1)
    X <- t(cbind(mat, cmat))
    X <- t(limma::removeBatchEffect(t(X), batch = c(rep("TCGA", dim(mat)[2]),
                                                    rep("clinical", dim(cmat)[2]))))
    pmat <- list(
      mat = t(X[1:dim(mat)[2],]),
      cmat = t(X[(dim(mat)[2]+1):dim(X)[1],])
    )
  }
  return(pmat)
}
