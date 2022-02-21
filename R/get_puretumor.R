#' Purify subgroup profiles
#'
#' @param subgroup The profiles of subgoups and control groups for each clinical query sample.
#' @param pp a GxM matrix, representing the expression profiles whose convex combination form the prior over the purified cancer profile learned.
#' @param minkappa (optional) The minimum value allowed for the strength parameter kappa placed over the reference cancer profile m (see Quon et al, 2013). By default, this is set to 1/min(BB), such that the log likelihood of the model is always finite. However, when the min(BB) is very small, this forces MIN_KAPPA to be very large, and can sometimes cause the reference profile m to look too much like a 'normal profile' (and therefore you may observe the tumor samples having low \% cancer content estimates). If this is the case, you can try setting MIN_KAPPA=1, or some other small value. For reference, for the data presented in Quon et al., 2013, MIN_KAPPA is on the order of 10^5.
#' @param loglevel (optional) A string that gives the logging threshold for futile.logger. The possible options are 'TRACE', 'DEBUG', 'INFO', 'WARN', 'ERROR', 'FATAL'. Currently the messages in ISOpureR are only in the categories 'INFO', 'WARN', and 'FATAL', and the default setting is 'INFO'. Setting a setting for the entire package will over-ride the setting for a particular function.
#'
#' @return The purified profiles of subgoups and control groups for each clinical query sample.
#' @export
#' @importFrom ISOpureR ISOpure.step1.CPE ISOpure.step2.PPE

get_puretumor = function(subgroup, pp = NULL, minkappa = NULL, loglevel = "INFO") {

  pure_analysis = function(i){
    message(paste0(' ... Begin purify',' ',names(subgroup)[i]))
    S1model <- suppressWarnings(ISOpureR::ISOpure.step1.CPE(subgroup[[i]]$case+1, subgroup[[i]]$control+1, PP = pp, MIN_KAPPA = minkappa, logging.level = loglevel))
    S2model <- suppressWarnings(ISOpureR::ISOpure.step2.PPE(subgroup[[i]]$case+1, subgroup[[i]]$control+1, S1model, MIN_KAPPA = minkappa, logging.level = loglevel))
    tumor <- S2model$cc_cancerprofiles
    dimnames(tumor) <- dimnames(subgroup[[i]]$case)
    data = list(tumor,subgroup[[i]]$control)
    names(data) = c('case_puify','control')
    return(data)
  }

  purify_data = lapply(seq(length(subgroup)), function(x){pure_analysis(x)})
  return(purify_data)
}
