#' Implement transcriptome subtyping
#'
#' @param mat Profiles of selected TCGA data set.
#' @param method Method for feature selection, optional 'MAD', 'VAR' and 'PCA' ('MAD' by default).
#' @param value A numeric value for biological feature selection.
#' If method='MAD' or 'VAR', the top number of value features are selected.
#' If method='PCA', the value of principal component is selected.
#' @param clusterNum Number of subtypes for the data set.
#' @param rank A numeric vector for the estimation of clusterNum. It will be overrided if clusterNum is not NULL.
#' @param nrun Number of runs to perform NMF. A default of 30 runs are performed, allowing the computation of a consensus matrix that is used in selecting the best result for cancer subtypes identification as Consensus Clustering method.
#' @param seed A numeric value is used to seed the random number generator before generating a random starting point.
#' @param doPlot A logical value. If true, draw the heatmap for the distance matrix with samples ordered to form clusters.
#'
#' @return The NMF subtyping result.
#' @export
#' @importFrom NMF nmf consensusmap
#' @importFrom CancerSubtypes silhouette_SimilarityMatrix
#' @importFrom grDevices palette


get_NMF = function(mat, method='MAD', value=1500, clusterNum=NULL, rank=2:6, nrun=30,
                   seed=NULL, doPlot=TRUE) {

  if (!method %in% c('MAD','VAR', 'PCA')){
    stop("Method should be one of MAD, VAR, and PCA")
  }

  message('... Selecting gene feature ... ')

  if (method=='MAD') {
    mads = apply(mat, 1, mad)
    feature_num = length(mads)
    index = sort(mads, decreasing = TRUE, index.return = TRUE)
    if (value < 1) {
      stop("Value should be number of genes with method 'MAD', .")
    }
    if (value > nrow(mat)) {
      value = nrow(mat)
      message("Warning: the feature selection number is beyond the original feature number.")
    }
    index = index$ix[1:value]
    selectmat = mat[index, ]
  }

  if (method=='VAR') {
    vars = apply(mat, 1, var)
    feature_num = length(vars)
    index = sort(vars, decreasing = TRUE, index.return = TRUE)
    if (value < 1) {
      stop("Value should be number of genes with method 'VAR'.")
    }
    if (value > nrow(mat)) {
      value = nrow(mat)
      message("Warning: the feature selection number is beyond the original feature numnber.")
    }
    index = index$ix[1:value]
    selectmat = mat[index, ]
  }

  if (method=='PCA') {
    if (value > 1) {
      stop("Value should be number of the ratio of principal component with method 'PCA'.")
    }
    newmat = t(mat)
    aa = prcomp(newmat, scale = TRUE)
    vars <- apply(aa$x, 2, var)
    props <- vars/sum(vars)
    xx = as.vector(cumsum(props))
    num = which(xx > value)[1]
    coeff = aa$ro
    score = (newmat %*% coeff)[, 1:num]
    selectmat = t(score)
  }

  if (any(selectmat < 0 )) {
    (selectmat <- (selectmat*2)) & (selectmat[selectmat < 0] <- 0)
  }

  if (any(is.na(selectmat))) {
    selectmat[is.na(selectmat)] <- 0
  }

  if (is.null(clusterNum)) {
    message('... Estimating clusterNum ...')
    res <- NMF::nmf(x=selectmat, rank=rank, nrun=nrun, seed=seed)
    num <- rank
    clusterNum <- num[which.max(res[["measures"]][["cophenetic"]])]
  }

  message('... excuting nmf ...')
  if (is.list(selectmat)) {
    temp = NULL
    for (i in 1:length(selectmat)) {
      temp = rbind(temp, selectmat[[i]])
    }
  }
  temp <- selectmat
  data1 <- rbind(pmax(temp, 0), -pmin(temp, 0))
  index <- which(rowSums(data1) == 0)
  data1 <- data1[-index, ]
  res <- NMF::nmf(data1, rank = clusterNum, nrun = nrun)
  distanceMatrix <- slot(res, "consensus")
  attr(distanceMatrix, "class") <- "Similarity"
  group <- as.numeric(as.vector(predict(res)))
  result <- list(group = group,
                 distanceMatrix = distanceMatrix,
                 originalResult = res)
  sil <- CancerSubtypes::silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
  result$sil <- sil

  if (isTRUE(doPlot)) {
    distanceMatrix = result$distanceMatrix
    group = result$group
    attr(distanceMatrix, "class") = NULL
    ind = order(group, -sil[, "sil_width"])
    num = length(unique(group))
    Var1 = c(palette()[2:(num + 1)])
    names(Var1) = sort(unique(group))
    ann_colors = list(group = Var1)
    annotation = data.frame(group = as.factor(group))
    NMF::consensusmap(distanceMatrix, Rowv = ind, Colv = ind,
                      main = "Clustering display", annCol = annotation$group,
                      annColors = ann_colors, labRow = "Sample", labCol = "Sample",
                      scale = "none")
  }
  return(result)
}
