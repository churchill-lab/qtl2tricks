#' Given a table of LOD peaks and the data used to generate them, calculate the
#' founder allele effects for each peak.
#'
#' @param lod.peaks data.frame containing at least 3 columns named annot.id, marker.id and lod. We used annot.id to find the phenotype in expr and marker.id to find the genoprobs at the QTL.
#' @param genoprobs list of 3 dimensional arrays that contain the allele probabilities for the mice. In qtl2 format.
#' @param pheno matrix of phenotypes with samples in rows and phenotypes in columns. In qtl2 format.
#' @param K list of numeric kinship matrices.
#' @param addcovar numeric matrix of additive covariates, created using model.matrix.
#' @param map list of numeric vectors containing marker locations in qtl2 format.
#' @return data.frame containing the lod.peaks with the eight founder alleles, in columns A through H, appended.
#' @import qtl2
#' @export
add_blup_to_lod_peaks = function(lod.peaks, genoprobs, pheno, K, addcovar, map) {

  # Unlist the map.
  markers = data.frame(marker = unlist(sapply(map, names)), 
                       chr    = rep(names(map), sapply(map, length)),
                       pos    = unlist(map),
                       stringsAsFactors = FALSE)

  # Append founder allele effect columns to the lod.peaks.
  lod.peaks = cbind(lod.peaks, matrix(0, nrow(lod.peaks), 8, dimnames = 
                    list(rownames(lod.peaks), LETTERS[1:8])))

  for(i in 1:nrow(lod.peaks)) {

    mkr = lod.peaks$marker.id[i]
    id  = lod.peaks$annot.id[i]
    chr = markers$chr[mkr == markers$marker]

    gp = genoprobs[,chr]
    gp[[1]] = gp[[1]][,,mkr,drop = FALSE]

    ph = pheno[,id,drop = FALSE]

    blup = scan1blup(genoprobs = gp, pheno = ph, kinship = K[[chr]], 
           addcovar = addcovar)

    lod.peaks[i,LETTERS[1:8]] = blup[1,LETTERS[1:8]]

  } # for(i)

  return(lod.peaks)

} # add_blup_to_lod_peaks()
