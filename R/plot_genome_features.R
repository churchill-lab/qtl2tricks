#' Given a set of genomic features (genes, SNPs, etc.), plot them on a karyogram.
#'
#' @param features GenomicRanges object containing the feature locations to plot.
#' @param col color to use for the features.
#' @return karyogram of the genome with the features plotted.
#' @import ggplot
#' @import ggbio
#' @export
plot_genome_features = function(features, col) {

  autoplot(features, layout = "karyogram", 
           color = col, fill = col) +
           theme(legend.position = "none")

} # plot_genome_features()







