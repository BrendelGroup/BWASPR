#' plot_diff_genes_pattern()
#' ??This function will plot out the mC level pattern in genes that were identified with diff methylated C sites
#'
#' @param mcobj a methylRawList object
#' @param diff_genes diff genes
#' @param output_figure dir/file name of the output figure
#'
#' @return A Granges object that contains a list of genes that have diff methylated C sites
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",type="CpGhsm",
#'                        mincvrg=1,assembly="Amel-4.5")
#'   genome <- read_genome_info(myfiles$parameters)
#'   meth_diff <- get_diff_genes(AmHE, genome)
#'   plot_diff_genes_pattern <- function(AmHE, meth_diff)
#'
#' @export

plot_diff_genes_pattern <- function(mrobj, diff_genes, output_figure = 'meth_diff_genes.pdf') {
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
