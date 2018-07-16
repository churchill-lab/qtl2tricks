#' Get ensembl GTF as a GRanges object.
#'
#' @param version integer that is the Ensembl version to get. Default = NULL which gets the most recent version.
#' @return GRanges object taht is the Ensembl GTF.
#'
#' @import AnnotationHub
#' @import GenomicRanges
#' @export
get_ensembl_genes = function(version = NULL) {

    hub = AnnotationHub()
    hub = query(hub, c("ensembl", "gtf", "mus musculus", "GRCm38"))
    hub = hub[grep("^Mus_musculus\\.GRCm38\\.[0-9]+\\.gtf$", hub$title)]

    gtf.title = NULL
    if(is.null(version)) {
      gtf.title = sort(hub$title)[length(hub)]
    } else {
      gtf.title = hub$title[grep(version, hub$title)]
    } # else

    if(length(gtf.title) == 0) {
      stop(paste("Ensembl version", version, "not found on AnnotationHub."))
    }

    ensembl = hub[[names(hub)[hub$title == gtf.title]]]

    return(ensembl[ensembl$type == "gene"])

} # get_ensembl_genes()


#' Estimate the chromosome lengths from the location of the gene with the highest end position.
#' 
#' @param ensembl GRanges object containing the Ensembl GTF.
#' @return numeric vector of chromosome lengths for the autosomes, sex and mitochondrial chromosomes.
#'
#' @import dplyr
#' @import GenomicRanges
#' @export
get_chr_length = function(ensembl) {

  tmp = data.frame(chr = seqnames(ensembl), end = end(ensembl))
  tmp = tmp %>% filter(!(substring(chr, 1, 1) %in% c("G", "J"))) %>%
          group_by(chr) %>% 
          summarize(len = max(end))
  chrlen = tmp$len * 1e-6
  names(chrlen) = tmp$chr

  return(chrlen)

} # get_chr_length()



#' Make a transcriptome map using the ggplot2 framework.
#' 
#' @param data: data.frame (or tibble) with the following columns:
#'       ensembl: (required) character string containing the Ensembl gene ID.
#'       qtl_chr: (required) character string containing QTL chromsome.
#'       qtl_pos: (required) floating point number containing the QTL position 
#'                in Mb.
#'       qtl_lod: (optional) floating point number containing the LOD score.
#'       gene_chr:  (optional) character string containing transcript chromosome.
#'       gene_start: (optional) character string containing transcript start 
#'                 postion in Mb.
#'       gene_end:  (optional) character string containing transcript end
#'                position in Mb.
#' color.points: logical that is TRUE if the points should be colored by LOD.
#' cis.points: logical that is TRUE if the points should be colored if they
#'             are with in cis.
#' cis.radius: numeric value containing the radius in Mb between a gene and a cis-eQTL.
#'             Optional.
#' cis.color: color for cis QTL. Optional.
#'
#' @return plot of the QTL and gene location for each gene.
#' @import ggplot2
#' @import dplyr
#' @export
ggtmap = function(data, color.points = FALSE, cis.points = FALSE, cis.radius = 2, 
         cis.color = "#4286f4") {

  # Check for required column names.
  required.colnames = c("ensembl", "qtl_chr", "qtl_pos")
  
  if(all(!required.colnames %in% colnames(data))) {
    stop(paste("colnames must contain the following columns:", 
         paste(required.colnames, collapse = ",")))
  } # if(!required.colnames %in% colnames(data))

  # Make sure that columns are not factors.
  data$ensembl = as.character(data$ensembl)
  data$qtl_chr = as.character(data$qtl_chr)

  gene.position.colnames = c("gene_chr", "gene_start", "gene_end")
  if(!all(gene.position.colnames %in% colnames(data))) {

    message(paste("Using Ensembl gene locations because optional gene",
            "position columns (", paste0(gene.position.colnames, collapse = ","),
            ") not found."))

	# Get the latest Ensembl GTF.
    ensembl = get_ensembl_genes()
			
    id    = ensembl$gene_id
    chr   = seqnames(ensembl)
    start = start(ensembl) * 1e-6
    end   = end(ensembl)   * 1e-6

    df = data.frame(ensembl = id, gene_chr = chr, gene_start = start, gene_end = end,
         stringsAsFactors = F)
    data = left_join(data, df, by = "ensembl")

  } # if(gene.position.colnames %in% colnames(data))

  # Make sure that columns are not factors.
  data$gene_chr = as.character(data$gene_chr)

  # Get the gene mid-point.
  data = data %>% mutate(gene_pos = (gene_end + gene_start) * 0.5)

  # Fix the factor levels for the chr.
  all.chr = data %>% select(qtl_chr, gene_chr) %>%
            gather(k, v) %>%
            select(v) %>%
            distinct() %>%
            arrange(v)
  all.chr = all.chr$v[!is.na(all.chr$v)]

  if(length(grep("M", all.chr)) > 0) {
    wh = grep("M", all.chr)
    all.chr = all.chr[c(1:(wh-1), (wh+1):length(all.chr), wh)]
  }

  # Remove any NAs.
  data = data %>% na.omit

  data$qtl_chr  = factor(data$qtl_chr,  levels = all.chr[order(as.numeric(all.chr))])
  data$gene_chr = factor(data$gene_chr, levels = rev(all.chr[order(as.numeric(all.chr))]))

  # If we're plotting cis points, then add a cis-QTL column.
  if(cis.points) {
  
    data = data %>% mutate(cis = (gene_chr == qtl_chr) & (abs(gene_start - qtl_pos) <= cis.radius))
	cis.colors = c("black", cis.color)
	names(cis.colors) = c("FALSE", "TRUE")
    print(ggplot(data, aes(x = qtl_pos, y = gene_pos), alpha = 0.5) +
      geom_point(aes(color = cis), alpha = 0.5) + 
      scale_color_manual(values = cis.colors) +
      facet_grid(gene_chr ~ qtl_chr, scales = "free", shrink = TRUE) +
      theme(panel.background = element_blank(),
	        panel.border = element_rect(fill = 0, color = "grey70"),
	        panel.grid.minor = element_blank(),
            panel.spacing = unit(0.05, "lines"),
            axis.text.x = element_text(angle = 90, hjust = 1)))

  } else { 

    print(ggplot(data, aes(x = qtl_pos, y = gene_pos)) +
      geom_point(aes(color = qtl_lod, alpha = 0.5)) + {
        if(color.points) scale_color_continuous(low = "grey50", high = "red") 
      } +
      facet_grid(gene_chr ~ qtl_chr, scales = "free", shrink = TRUE) +
      theme(panel.background = element_blank(),
	  	    panel.border = element_rect(fill = 0, color = "grey70"),
	  	    panel.grid.minor = element_blank(),
            panel.spacing = unit(0.05, "lines"),
            axis.text.x = element_text(angle = 90, hjust = 1)))

   } # else

} # ggtmap()

