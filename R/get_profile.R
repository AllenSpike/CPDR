#' @import dplyr
#' @importFrom rhdf5 h5read H5close
#' @importFrom clusterProfiler bitr
get_profile = function(case_id, control_id, saveDir, OrgDB) {

  file = file.path(saveDir,'GTEX','octad.counts.and.tpm.h5')
  samples = rhdf5::h5read(file,"meta/samples")
  transcripts = rhdf5::h5read(file,"meta/transcripts")

  dz_counts = rhdf5::h5read(file,"data/count",index = list(which(!is.na(transcripts)),
                                                           which(samples %in% case_id)))
  dimnames(dz_counts) = list(transcripts[!is.na(transcripts)],
                             samples[samples %in% case_id])

  ct_counts = rhdf5::h5read(file,"data/count",index = list(which(!is.na(transcripts)),
                                                           which(samples %in% control_id)))
  dimnames(ct_counts) = list(transcripts[!is.na(transcripts)],
                             samples[samples %in% control_id])

  gene <- suppressMessages(suppressWarnings(clusterProfiler::bitr(row.names(dz_counts), fromType='ENSEMBL', toType='SYMBOL',OrgDb = OrgDB)))
  gene <- gene[!duplicated(gene$ENSEMBL),]
  dz_counts <- dz_counts[row.names(dz_counts) %in% gene$ENSEMBL,]
  dz_counts <- dz_counts[gene$ENSEMBL,]
  row.names(dz_counts) <- gene$SYMBOL

  ct_counts <- ct_counts[row.names(ct_counts) %in% gene$ENSEMBL,]
  ct_counts <- ct_counts[gene$ENSEMBL,]
  row.names(ct_counts) <- gene$SYMBOL

  dz_counts = 2^dz_counts-1
  ct_counts = 2^ct_counts-1
  data = cbind(dz_counts, ct_counts)

  dz_counts <- data[,case_id]
  ct_counts <- data[,control_id]
  return(list(dz_counts,ct_counts))
}
