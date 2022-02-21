#' Identify subgroup and control group for each clinical query sample
#'
#' @param cmat The RNA-seq count expression profiles of clinical query samples.
#' @param mat Profiles of selected TCGA data set.
#' @param subtype The NMF subtyping result.
#' @param k The subgroup size (10 by default. Optional 3, 5, 10.)
#' @param biopsy Tissue sources of the control group (Use View(GTEX_set) for optional 51 biopsy.sites. By default set to NULL, which will select control group based on distance between TCGA and GTEX samples.)
#' @param adjacent A logical value. By default set to FALSE. If TRUE, only tissue with sample.type 'adjacent' from GTEX_set would be used instead of 'normal'.
#' @param doplot A logical value. If TRUE, draw the heatmap of Spearman correlation coefficients between query samples and subtypes.
#' @param OrgDb The human annotation db for ID convert.
#' @param db.path A character(1) string with the folder path where the background data was saved.
#'
#' @return Profiles of subgoups and control groups for each clinical query sample.
#' @export
#' @import ggplot2

get_subgroup = function(cmat, mat, subtype, k = 10, biopsy = NULL, adjacent = FALSE,
                        doplot = TRUE, OrgDB = NULL, db.path = tempdir()) {

  adj = adjacent
  bio = biopsy
  db.path <- file.path(db.path, 'CPDR_db')
  if (!is.matrix(cmat)) {
    cmat = as.matrix(cmat)
  }
  if (!identical(colnames(subtype$distanceMatrix), colnames(mat))) {
    stop('The patient names do not match in `mat` and `subtype`.')
  }
  if (!k %in% c(3,5,10)) {
    stop('Either k=3, k=5 or k=10.')
  }


  message('... Calculating correlation coefficient ... ')
  cor = apply(cmat, 2, function(x){get_pmat(mat,x,topk = 1500)})
  cor_m = lapply(cor, function(x){
    x = x[["median.cor"]][colnames(subtype$distanceMatrix)]
    cor = c()
    for (i in 1:max(subtype$group)) {
      cor = c(cor, mean(x[which(subtype$group == i)]))
    }
    return(cor)})
  cor_mm = as.matrix(as.data.frame(cor_m))

  cor_n = list()
  for (i in 1:dim(cmat)[2]) {
    group = unlist(lapply(cor_m, function(x){which.max(x)}))[i]
    x = cor[[i]][["median.cor"]][colnames(subtype$distanceMatrix)]
    x = x[which(subtype$group == group)]
    cor_n[[i]] = x
  }
  names(cor_n) = colnames(cmat)


  case_id = lapply(cor_n, function(x){names(sort(x))[1:k]})
  control_id = lapply(case_id, function(x){
    return(get_ref_tissue(x,
                          biopsy=bio,
                          source="octad",
                          control_size=length(x),
                          adjacent=adj))
  })

  if (isTRUE(doplot)) {
    cor_x = stack(cor_m)
    colnames(cor_x) = c('cor','samples')
    cor_x$group = rep(1:max(subtype$group),dim(cmat)[2])
    grid.newpage()
    plot = ggplot(data = cor_x,mapping = aes(x=samples,y=group,fill=cor))+
           geom_tile(color = "white")+
           scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = mean(cor_x$cor), limit = c(min(cor_x$cor),max(cor_x$cor)), space = "Lab",
                           name="Spearman\ncor") +
           theme_minimal()+
           theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                       size = 12, hjust = 1))+
           coord_fixed()+
           geom_text(aes(samples, group, label = round(cor,2)), color = "black", size = 4)
    print(plot)
  }


  message('... Get normalized profiles from octad.db ...')
  tmp = lapply(seq(length(case_id)), function(x){
    tz = suppressMessages(get_profile(case_id[[x]],control_id[[x]],db.path,OrgDB))
    names(tz) = c('case','control')
    return(tz)
  })
  names(tmp) = colnames(cmat)
  return(tmp)
}
