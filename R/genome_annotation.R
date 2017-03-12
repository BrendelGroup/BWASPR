#' genome_annotation()
#' annotate the methylated sites
#'
#' @param mrobj
#' @param genome_info
#'
#' @return A data frame that store the genome annotation data
#'
#' @importFrom
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles)
#'   genome <- read_genome_info(myfiles$parameters)
#'   genome_annotation_info <- genome_annotation(AmHE, genome, save_output_data = FALSE)
#'   genome_annotation(AmHE, genome, outputdata = 'AmHE_genome_annotation.Rda')
#'
#' @export


genome_annotation <- function(mrobj, genome_info, save_output_data = TRUE, outputdata = 'output_data.Rda'){

    meth <- unite(mrobj, destrand = TRUE)
    matrix <- percMethylation(meth, rowids = FALSE, save.txt = FALSE)

    for (i in names(genome_info)){
      assign(paste(i, 'annot', sep = '_'),
             annotateWithFeature(as(meth, 'GRanges'), genome_info[[i]], strand = TRUE, intersect.chr = FALSE))
    }

    for (i in names(genome_info)){
      meth[i] <- get(paste(i, 'annot', sep = '_'))@members
    }

	  meth_data <- getData(meth)
	  label <- subset(meth_data, select = c('chr', 'start', 'end', 'strand', names(genome_info)))
  	genome_annotation <- cbind(matrix, label)

    if (save_output_data == TRUE){

	    save(genome_annotation, file = outputdata)
    }

  	else {
  	  return(genome_annotation)
  	}
}
