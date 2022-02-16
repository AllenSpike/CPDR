# @importFrom    rhdf5 h5read H5close
get_profile = function(case_id, control_id, saveDir) {

  file = file.path(saveDir,'GTEX')
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

  dz_counts = 2^dz_counts-1
  ct_counts = 2^ct_counts-1
  data = cbind(dz_counts, ct_counts)

  dz_counts <- data[,case_id]
  ct_counts <- data[,control_id]
  return(list(dz_counts,ct_counts))
}
