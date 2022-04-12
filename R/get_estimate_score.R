#' Evaluatting non-cancer cell infiltration of each subtype
#'
#' @param mat Preprocessed Profiles of TCGA dataset.
#' @param subtype Subtyping result.
#' @param gs A data.frame of 3 columns with Metagene, Cell_type and Annot.
#' @param doplot A logical value. If TRUE, draw the heatmap of non-cancer cell infiltration.
#'
#' @return A data.frame of GSVA enrichment scores.
#' @importFrom GSVA gsva
#' @importFrom pheatmap pheatmap
#' @importFrom BiocGenerics Reduce
#' @import dplyr
#' @export

get_estimateScore = function(mat, subtype, gs = NULL, doplot = TRUE) {

  if (!is.null(gs)) {
    if (all(colnames(gs) != c("Metagene","Cell type","Annot"))) {
      stop('gs should be a dataframe with Metagene, Cell type and Annot columns')
    }
    gene_set = gs
  }

  if (sum(row.names(mat) %in% gene_set$Metagene) == 0) {
    stop('Row names of mat should be SYMBOL.')
  }
  if (!identical(colnames(subtype$distanceMatrix), colnames(mat))) {
    stop('The patient names in mat and subtype do not match.')
  }


  list <- split(as.matrix(gene_set)[,1], gene_set[,2])
  gsva_matrix <- GSVA::gsva(as.matrix(mat), list, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

  gsva_matrix1<- t(scale(t(gsva_matrix)))
  gsva_matrix1[gsva_matrix1< -2] <- -2
  gsva_matrix1[gsva_matrix1>2] <- 2

  annot = table(gene_set$Annot)
  n = length(annot)
  ls1 = lapply(seq(n), function(i){gsub('^ ','',rownames(gsva_matrix1))%in%unique(gene_set$`Cell type`[which(gene_set$Annot == names(annot[i]))])})
  ls2 = lapply(ls1, function(i){gsva_matrix1[i,]})
  gsva_matrix1 <- BiocGenerics::Reduce(rbind, ls2)

  normalization<-function(x){
    return((x-min(x))/(max(x)-min(x)))}
  nor_gsva_matrix1 <- normalization(gsva_matrix1)


  annotation_col <- data.frame(group=subtype$group,
                               patient=colnames(subtype$distanceMatrix))
  annotation_col <- arrange(annotation_col,group)
  nor_gsva_matrix1 <- nor_gsva_matrix1[,annotation_col$patient]
  annotation_col = data.frame(group=as.character(annotation_col$group))
  rownames(annotation_col)<-colnames(nor_gsva_matrix1)

  if (doplot) {
    bk = unique(c(seq(0,1, length=100)))
    gb = group_by(gene_set,Annot)
    gr = as.data.frame(dplyr::summarise(gb,x=length(unique(`Cell type`))))
    row.names(gr) = gr$Annot
    gr = gr[names(annot),]
    grid.newpage()
    pheatmap::pheatmap(nor_gsva_matrix1,
                       show_colnames = F,
                       cluster_rows = F,
                       cluster_cols = F,
                       annotation_col = annotation_col,
                       breaks=bk,
                       fontsize=6,
                       gaps_row = cumsum(gr$x)[-length(gr$x)])
  }

  return(nor_gsva_matrix1)
}
