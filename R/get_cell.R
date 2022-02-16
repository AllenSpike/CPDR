#' Compute individual-related cell line
#'
#' @param cmat The RNA-seq count expression profiles of clinical query samples.
#' @param db.path A character(1) string with the folder path where the background data was saved.
#' @param removeBatchEffect Remove batch effects from expression data using limma (TRUE by default).
#' @param orgDB The human annotation db for ID convert.
#'
#' @return A cell line information list.
#' @export
#' @importFrom PharmacoGx molecularProfiles cellInfo
#' @importFrom foreach foreach
#' @importFrom CoreGx sensitivityProfiles
#' @importFrom clusterProfiler bitr
#' @importFrom limma removeBatchEffect
#' @import dplyr

get_cell = function(cmat, db.path = tempdir(), removeBatchEffect = TRUE, orgDB) {

  pick.out.cell.line <- function(expr.of.samples, expr.of.cell.lines,
                                 marker.gene) {
    if (is.vector(expr.of.samples)) {
      marker.gene <- intersect(names(expr.of.samples),
                               (marker.gene))
      marker.gene <- intersect(rownames(expr.of.cell.lines),
                               (marker.gene))
      correlation.matrix <- cor(expr.of.samples[marker.gene],
                                expr.of.cell.lines[marker.gene, ], method = "spearman")
    } else {
      marker.gene <- intersect(rownames(expr.of.samples),
                               (marker.gene))
      marker.gene <- intersect(rownames(expr.of.cell.lines),
                               (marker.gene))
      correlation.matrix <- cor(expr.of.samples[marker.gene,
      ], expr.of.cell.lines[marker.gene, ], method = "spearman")
    }
    correlation.matrix[is.na(correlation.matrix)] = 0
    cell.line.median.cor <- apply(correlation.matrix, 2,
                                  median) %>% sort(decreasing = TRUE)
    best.cell.line <- names(cell.line.median.cor)[1]
    p.value.vec <- foreach::foreach(cell.line = setdiff(names(cell.line.median.cor),
                                                        best.cell.line), .combine = "c") %do% {
                                                          v <- correlation.matrix[, cell.line]
                                                          p.value <- wilcox.test(correlation.matrix[, best.cell.line],
                                                                                 v, alternative = "greater", paired = TRUE)$p.value
                                                        }
    names(p.value.vec) <- setdiff(names(cell.line.median.cor),
                                  best.cell.line)
    fdr.vec <- p.adjust(p.value.vec, method = "fdr")
    list(cell.line.median.cor = cell.line.median.cor, best.cell.line = best.cell.line,
         compare.fdr.vec = fdr.vec, correlation.matrix = correlation.matrix)
  }

  choosecell_indv = function(case_id, CCLE.log2.read.count.matrix, case2, cmat, orgDB){

    ###
    case_counts = as.matrix(cmat[,case_id])
    colnames(case_counts) = case_id
    row.names(case_counts) = row.names(cmat)

    ###
    print(paste("computing correlation between cell lines and selected sample",
                sep = " "))
    colnames(CCLE.log2.read.count.matrix) = case2$cellid
    CCLE.median = apply(CCLE.log2.read.count.matrix,
                        1, median)

    CCLE.expressed.gene <- names(CCLE.median)[CCLE.median >
                                                1]
    tmp <- CCLE.log2.read.count.matrix[CCLE.expressed.gene,
    ]
    tmp.rank <- apply(tmp, 2, rank)
    rank.mean <- apply(tmp.rank, 1, mean)
    rank.sd <- apply(tmp.rank, 1, sd)
    CCLE.rna.seq.marker.gene.1000 <- names(sort(rank.sd, decreasing = TRUE))[1:1000]


    TCGA.vs.CCLE.polyA.expression.correlation.result <- pick.out.cell.line(case_counts,
                                                                           CCLE.log2.read.count.matrix, CCLE.rna.seq.marker.gene.1000)
    correlation.dataframe <- TCGA.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor %>%
      as.data.frame()
    colnames(correlation.dataframe) <- "cor"
    correlation.dataframe$cellid = row.names(correlation.dataframe)
    return(correlation.dataframe)
  }

  get_cellline_info = function(cell_line_computed,prism){
    df_sp <- CoreGx::sensitivityProfiles(prism)
    df_sp$key <- row.names(df_sp)
    df_in <- unique(CoreGx::sensitivityInfo(prism)[,c(7,9,10)])
    df_in$key <- row.names(df_in)
    df_al <- left_join(df_sp, df_in, by='key')
    df_al <- filter(df_al, PRISM.drugid %in% res_prism$PRISM_drug_name) %>% left_join(res_prism[,c(2,3)],by=c('PRISM.drugid'='PRISM_drug_name'))

    if(class(cell_line_computed)=="list"){
      cell_line_computed <- lapply(cell_line_computed, function(x){
        return(    x %>% mutate(
          drugnum = unlist(lapply(cellid, function(y){sum(df_al$cellid == y,na.rm = T)}))
        ))
      })
    } else {
      cell_line_computed <- cell_line_computed %>% mutate(
        drugnum = unlist(lapply(cellid, function(y){sum(df_al$cellid == y,na.rm = T)}))
      )
    }

    tmp = list(df_al, cell_line_computed)
    return(tmp)
  }


  ##
  if (is.null(cmat)) {
    stop("Patient profile input not found.")
  }
  if (is.null(db.path)) {
    stop("db.path input not found.")
  }

  message("loading background data......")
  prism <- suppressMessages(readRDS(file.path(db.path,'PRISM.rds')))
  ccle <- suppressMessages(readRDS(file.path(db.path,'CCLE.rds')))
  CCLE.log2.read.count.matrix = PharmacoGx::molecularProfiles(ccle, 'rna')
  case = intersect(colnames(CCLE.log2.read.count.matrix),toupper(PharmacoGx::cellInfo(ccle)$Expression.arrays))
  case2 = dplyr::filter(PharmacoGx::cellInfo(ccle), toupper(Expression.arrays) %in% case)

  gene = suppressMessages(clusterProfiler::bitr(row.names(cmat),fromType = 'SYMBOL',toType = 'ENSEMBL',OrgDb = orgDB))
  gene = gene[!duplicated(gene$SYMBOL),]
  cmat = cmat[gene$SYMBOL,]
  row.names(cmat) = gene$ENSEMBL

  CCLE.log2.read.count.matrix = CCLE.log2.read.count.matrix[intersect(row.names(CCLE.log2.read.count.matrix),
                                                                      row.names(cmat)),
                                                            toupper(case2$Expression.arrays)]
  cmat = cmat[intersect(row.names(CCLE.log2.read.count.matrix),
                        row.names(cmat)),]

  if (isTRUE(removeBatchEffect)) {
    message('... Removing batch effect ...')
    X <- t(cbind(cmat, CCLE.log2.read.count.matrix))
    X <- t(limma::removeBatchEffect(t(X), batch = c(rep("cmat", dim(cmat)[2]),
                                                    rep("ccle", dim(CCLE.log2.read.count.matrix)[2]))))

    cmat = t(X[1:dim(cmat)[2],])
    CCLE.log2.read.count.matrix = t(X[(dim(cmat)[2]+1):dim(X)[1],])
  }

  res=lapply(colnames(cmat), function(x){
    return(choosecell_indv(case_id=x,
                           CCLE.log2.read.count.matrix=CCLE.log2.read.count.matrix,
                           case2=case2,
                           cmat=cmat,
                           orgDB=orgDB))
  })

  message("collecting cell line info......")
  tmp <- get_cellline_info(res,prism)
  sensp <- tmp[[1]]
  cell_line_computed <- tmp[[2]]
  names(cell_line_computed) <- colnames(cmat)

  tmp <- list(
    correlation = cell_line_computed,
    sensp = sensp
  )
  gc()
  return(tmp)
}




