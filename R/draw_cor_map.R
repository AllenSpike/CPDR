#' Visualize drug effectiveness test result
#'
#' @inheritParams drugcorTest
#' @export
#' @importFrom data.table dcast.data.table
#' @importFrom reshape2 melt
#' @importFrom ggpubr stat_cor
#' @import ggplot2
#' @importFrom grid grid.newpage pushViewport viewport

draw_cor_map = function(mysRGES, topline, cell_info){

  message('data preprocess......')
  sensp = cell_info$sensp
  auc <- data.table::dcast.data.table(data.table::as.data.table(sensp),
                                      cid ~ cellid, value.var = "auc_published", fun.aggregate = median)
  auc = as.data.frame(auc)

  top = lapply(names(table(topline)), function(j){

    AUC = subset(auc, select = c("cid", j))
    auc.m = as.data.frame(reshape2::melt(AUC, id.vars = "cid"))
    auc.medianauc = aggregate(auc.m[3], by = list(auc.m$cid),
                              FUN = median)
    auc.medianauc$medauc = auc.medianauc$value
    auc.medianauc$value = NULL
    auc.medianauc$DRUGID = auc.medianauc$Group.1
    auc.medianauc$Group.1 = NULL
    auc.medianauc <- auc.medianauc[is.finite(auc.medianauc$medauc),
    ]

    auc.medianauc <- dplyr::mutate(auc.medianauc,
                                   effect  = ifelse(medauc<(median(medauc)-0.5*sd(medauc)), 'effective', 'ineffective'))

    return(list(auc.medianauc))

  })
  names(top) = names(table(topline))
  mysRGES <- merge(mysRGES, res_prism, by.x = "pert_iname", by.y = "LINCS_drug_name")
  testdf <- merge(mysRGES, top[[1]], by.x = "cid",
                  by.y = "DRUGID")

  dot <- ggplot(data = testdf, aes(x = sRGES, y = medauc)) +
    geom_point(aes(color=effect)) +
    geom_smooth(method = "lm", se = T, color = "black", size = 0.5)+
    ggpubr::stat_cor(method = "pearson") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 15, colour='black',face = "bold"),
                   axis.text.y = element_text(size = 15, colour='black',face = "bold"),
                   axis.title.x = element_text(size = 15, colour='black',face = "bold"),
                   axis.title.y = element_text(size = 15, colour='black',face = "bold")) +
    theme(panel.grid = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(size=1, colour = "black")) +
    theme(legend.position="none")

  box <- ggplot(data = testdf, aes(x = sRGES, y = effect)) +
    geom_boxplot(aes(color=effect)) +
    theme(axis.title.y  = element_blank()) +
    ggtitle(paste0('p = ', signif(t.test(medauc~effect,testdf)[["p.value"]],3)))+
    theme_bw() +
    theme(axis.text.x = element_text(size = 15,colour='black',face = "bold"),
                   axis.text.y = element_text(size = 15, colour='black',face = "bold"),
                   axis.title.x = element_text(size = 15, colour='black',face = "bold")) +
    ylab(NULL) +
    theme(panel.grid = element_blank()) +
    theme(panel.border = element_blank()) +
    theme(axis.line = element_line(size=1, colour = "black")) +
    theme(legend.position="none")


  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(4,1)))
  vplayout <- function(x,y){
    grid::viewport(layout.pos.row = x, layout.pos.col = y)
  }
  print(dot, vp = vplayout(1:3,1))
  print(box, vp = vplayout(4,1))
}

