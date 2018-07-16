#' Given a list of Ensembl gene IDs, determine which Ensembl version has the highest number of matches.
#' 
#' @param ensids character vector of ensemble gene IDs.
#' 
#' @return data.frame containing the Ensembl version and the proportion of matching IDs.
#' 
#' @importFrom AnnotationHub AnnotationHub
#' @export
determine_ensembl_version = function(ensids) {

  ensids = as.character(ensids)

  hub = AnnotationHub()
  hub = hub[grep("^Mus_musculus.GRCm38\\.[0-9][0-9]\\.gtf$", hub$title)]

  result = data.frame(version = hub$title, prop_match = rep(0, length(hub)))

  for(i in seq_along(hub)) {

    ensembl = hub[[names(hub)[i]]]
    result$prop_match[i] = mean(ensids %in% ensembl$gene_id)
    print(result[i,])

  } # for(i)

} # determine_ensembl_version()



