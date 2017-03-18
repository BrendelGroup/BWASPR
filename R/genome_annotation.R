#' genome_annotation()
#' This function annotates the methylated sites with different genomic features (e.g. genes, exon, promoter, etc.)
#'
#' @param mrobj A methylRaw object or a methylRawList object.
#' @param genome_info A list of GRanges objects that contains genome infomation.
#' @param save_output_data If specified other than the default "TRUE", then results
#'   are saved in Rda file that is specified by @param outputdata; otherwise either no Rda file is generated.
#' @param outputdata If specified, then the result is saved in the specified file name; otherwise the result
#'   is saved in '$wd/genome_annotation_data.Rda'
#'
#' @return A table ?that stores the methylation information and genome_annotation
#'
#' @importFrom methylKit percMethylation unite
#' @importFrom genomation annotateWithFeature
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


genome_annotation <- function(mrobj, genome_info, save_output_data = TRUE, outputdata = 'genome_annotation_data.Rda'){
    # unite all the methlRawList object into a single table and calculate the ?'percentage of methylation'
    #
    meth   <- unite(mrobj, destrand = TRUE)
    matrix <- percMethylation(meth, rowids = FALSE, save.txt = FALSE)
    # ?annotate with generic features from genome information
    #
    for (i in names(genome_info)){
        assign(paste(i, 'annot', sep = '_'),
               annotateWithFeature(as(meth, 'GRanges'), genome_info[[i]], strand = TRUE, intersect.chr = FALSE)
              )
    }
    # ?add the annotation to the meth table
    #
    for (i in names(genome_info)){
        meth[i] <- get(paste(i, 'annot', sep = '_'))@members
    }
    #
	meth_data <- getData(meth)
	label <- subset(meth_data, select = c('chr', 'start', 'end', 'strand', names(genome_info)))
  	genome_annotation <- cbind(matrix, label)
    # ?save the table as Rda or return the table
  	#
    if (save_output_data == TRUE){
	    save(genome_annotation, file = outputdata)
    } else {
  	    return(genome_annotation)
  	}
}
