#' Download background data
#'
#' @param pset The Pharmacogenomic dataset name. View(pharmaco_sets) for currently available datasets.
#' @param tset The TCGA cancer type name. View(TCGA_sets) for currently available cancer types.
#' @param nset The normal tissue dataset name. Currently available datasets: GTEX.
#' @param saveDir A character(1) string with the folder path where the background data should be saved. Defaults to tempdir(). Will create directory if it does not exist.
#' @param verbose bool Should status messages be printed during download. Defaults to TRUE.
#'
#' @return
#' @export
#' @importFrom PharmacoGx availablePSets
#' @importFrom downloader download
#' @importFrom cBioPortalData downloadStudy

download_db = function( pset, tset, nset,
                        saveDir = tempdir(), verbose = FALSE) {

  saveDir <- file.path(saveDir, 'CPDR_db')

  if(!is.null(pset)){
    opts <- options()
    options(timeout = 600)
    on.exit(options(opts))
    pSetTable <- PharmacoGx::availablePSets()

    saveDir_pset <- paste0(saveDir,'/Pharmacogenomic')
    if (!file.exists(saveDir_pset)) {
      dir.create(saveDir_pset, recursive = TRUE)
    }

    if (any(is.na(match(pset, c('LINCS','CCLE','PRISM'))))) {
      stop("pset not found.Please use View(pharmaco_sets) for available pset.")
    }

    if(any(pset %in% c('CCLE','PRISM'))){
      sens <- pset[pset %in% c('CCLE','PRISM')]
      message(paste0("... Download drug sensitive data (",paste0(sens,collapse = ','),") ... "))
      lapply(sens, function(x){
        if (!file.exists(file.path(saveDir_pset, paste0(x, ".rds")))){
          whichx <- match(x, pSetTable[, "Dataset Name"])
          downloader::download(url = as.character(pSetTable[whichx,
                                                            "Download"]), destfile = file.path(saveDir_pset, paste0(x, ".rds")),
                               quiet = !verbose, mode = "wb")
        }
      })
    } else if(any(pset %in% 'LINCS')){
      message("... Download drug pertubation data (LINCS) ...")
      if (!file.exists(file.path(saveDir_pset, "LINCS.RData"))){
        downloader::download(url = 'https://zenodo.org/record/5880026/files/lincs_signature_new.RData?download=1',
                             destfile = file.path(saveDir_pset,"LINCS.RData"),
                             quiet = !verbose, mode = "wb")
      }
    }
  }


  if(!is.null(tset)){
    saveDir_tset <- paste0(saveDir,'/TCGA')
    TCGATable <- CPDR::TCGA_sets
    if (!file.exists(saveDir_tset)) {
      dir.create(saveDir_tset, recursive = TRUE)
    }
    if (any(is.na(match(tset, TCGATable[, 2])))) {
      stop("tset is not found.Please use View(TCGA_sets) for available tset")
    }
    tset <- TCGATable[match(tset, TCGATable[, 2]),3]
    message(paste0("... Download tcga data (",paste0(tset,collapse = ','),") ... "))
    lapply(tset, function(x){cBioPortalData::downloadStudy(x, use_cache=saveDir_tset)})
  }


  if(!is.null(nset)){
    saveDir_nset <- paste0(saveDir,'/GTEX')
    if (any(is.na(match(nset, 'GTEX')))) {
      stop("The nset is not found.Please use GTEX.")
    }
    if (!file.exists(saveDir_nset)) {
      dir.create(saveDir_nset, recursive = TRUE)
    }
    message(paste0("... Download normal tissue data (",paste0(nset,collapse = ','),") from OCTAD ... "))

    tryCatch(requireNamespace("octad.db"),
             error=function() {
               utils::install.packages("https://chenlab-data-public.s3.amazonaws.com/octad/octad.db_0.99.0.tar.gz%3Fdl%3D0",
                                        method="libcurl",repos=NULL,type="source")
             })
    downloader::download("https://chenlab-data-public.s3-us-west-2.amazonaws.com/octad/octad.counts.and.tpm.h5",
                         destfile = file.path(saveDir_nset, 'octad.counts.and.tpm.h5'), quiet = !verbose, mode = "wb")
  }
}
