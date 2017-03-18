#' plot_diff_genes_pattern()
#' This function will plot out the mC level pattern in genes that were identified with diff methylated C sites
#'
#' @param mcobj a methylRawList object
#'
#' @param diff_genes diff genes
#'
#' @param output_figure name of the output figure
#'
#' @return A Granges object that contains a list of genes that have diff methylated C sites
#'
#'
#' @examples
#'   meth_diff_genes <- get_diff_genes(mcobj)
#'
#' @export




plot_diff_genes_pattern <- function(mrobj, diff_genes, output_figure) {
  sample_list <- getSampleID(mrobj)

  # calculate the percentageCs
  for (i in seq_along(sample_list)){
    mrobj[[i]] <- mrobj[[i]] %>% mutate(percentageCs = round(numCs / coverage, digits = 2))
    mrobj[[i]] <- as(mrobj[[i]], 'GRanges')
  }
  # resize the diff_genes (GRanges)
  methdiffgenes_resize.gr <- resize(diff_genes, 20000)
  targets = GRangesList(mrobj)

  test <- ScoreMatrixList(targets,
                          methdiffgenes_resize.gr,
                          bin.num = 20,
                          # using differnt bin.op could be interesting
                          bin.op = 'mean',
                          strand.aware = FALSE,
                          weight.col ='percentageCs',
                          is.noCovNA = TRUE
  )

  test.sub <- intersectScoreMatrixList(test, reorder = FALSE)

  # plot the figure
  pdf(output_figure)
  options(repr.plot.width = 8, repr.plot.height = 5)
  multiHeatMatrix(test.sub,
                  col = topo.colors(10),
                  xcoords=c(0, 20000),
                  matrix.main=names(targets),
                  common.scale = TRUE,
                  xlab = sample_list,
                  legend = TRUE
  )
	dev.off()

	return(test.sub)
}
