#' Plot the heritability estimates and their confidene intervals.
#' 
#' @param h2 matrix of heritability simulations with simulations in rows and chromosomes in columns.
#' 
#' @import ggplot2
#' @export
herit_ci_plot = function(h2, title = "") {

  obs  = h2[1,]
  sims = h2[-1,]

  mean.h2 = colMeans(sims)
  q95 = apply(sims, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  q95 = t(q95)
  colnames(q95) = c("low", "high")

  # All data.
  df = data.frame(chr = c(1:19, "X"), mean.h2, q95, obs) %>%
       mutate(chr = factor(chr, levels = c(1:19, "X")))

  ggplot(data = df) +
    geom_pointrange(aes(chr, mean.h2, ymin = low, ymax = high)) +
    labs(title = title, y = "h^2")

} # herit_ci_plot()