#' Evaluate drug effectiveness on relevant cell line
#'
#' @param mysRGES sRGES data.frame produced by get_reverse_score.
#' @param topline inidividual-related cell lines.
#' @param cell_info cell line information list produced by get_cell.
#'
#' @return A correlation list
#' @export
#' @importFrom data.table dcast.data.table
#' @importFrom reshape2 melt

drugcorTest = function(mysRGES, topline, cell_info){


  auc <- data.table::dcast.data.table(data.table::as.data.table(cell_info$sensp),
                                      V2 ~ cellid, value.var = "auc_published", fun.aggregate = median)
  auc = as.data.frame(auc)

  res = lapply(seq(length(topline)), function(i){

    x = topline[i]
    sRGES = mysRGES[[i]]
    sRGES$pert_iname <- toupper(sRGES$pert_iname)

    auc = subset(auc, select = c("V2", x))
    auc.m = as.data.frame(reshape2::melt(auc, id.vars = "V2"))
    auc.medianauc = aggregate(auc.m[3], by = list(auc.m$V2),
                              FUN = median)
    auc.medianauc$medauc = auc.medianauc$value
    auc.medianauc$value = NULL
    auc.medianauc$DRUGID = auc.medianauc$Group.1
    auc.medianauc$Group.1 = NULL
    auc.medianauc <- auc.medianauc[is.finite(auc.medianauc$medauc),
    ]


    testdf2 <- merge(sRGES, auc.medianauc, by.x = "V2",
                     by.y = "DRUGID")
    AUC.cortest <- cor.test(testdf2$sRGES, testdf2$medauc, method='pearson')
    aucpval <- AUC.cortest$p.value
    aucrho <- AUC.cortest$estimate

    r = data.frame( rho = aucrho,
                    pval = aucpval)
    return(r)
  })

  names(res) = names(sRGES)
  return(res)

}

