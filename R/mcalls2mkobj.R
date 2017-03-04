#' mcalls2mkobj()
#' This function will select BWASP mcalls files according to specified
#' selection criteria and read the data into methylBase objects.
#'
#' @param inputdf A data frame as returned by setup_BWASPR from the input
#'   data file.  Each row corresponds to a BWASP-generated output, with columns
#'   indicating Species, Study, Sample, Replicate, Type, and File.
#' @param species A character string specifying the Species to be pulled out of
#'   inputdf; default: "all".
#' @param study A list of character strings specifying the Studies to be pulled
#'   out of inputdf; default: "all".
#' @param sample A list of character strings specifying the Samples to be pulled
#'   out of inputdf; default: "all".
#' @param replicate A list of integers specifying the Replicates to be pulled
#'   out of inputdf; default: c(1:20).
#' @param type The type should be CpGhsm or CpGscd, as per BWASP output.
#' @param assembly Label for the underlying genome assembly version; default:
#'   "unknown".
#'
#' @return A methylKit-package methylRaw or methylRawList object.
#'
#' @importFrom methylKit methRead
#' @importFrom utils write.table
#'
#' @examples
#    library("methylKit")
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   mcalls2mkobj(myfiles[[1]],species="all",study="HE",type="CpGhsm",
#'                assembly="Amel-4.5")
#'
#' @export

mcalls2mkobj <- function(inputdf,species="all",study="all",sample="all",
                         replicate=c(1:20),type="CpGhsm",assembly="unknown"){
    message("... loading mc objects ..")
    system("pwd")
    if (species != "all") {
        inputdf <- inputdf[inputdf$Species == species,]
    }
    if (study != "all") {
        inputdf <- inputdf[inputdf$Study %in% study,]
    }
    if (sample != "all") {
        inputdf <- inputdf[inputdf$Sample %in% sample,]
    }
    if (!identical(replicate,c(1:20))) {
        inputdf <- inputdf[inputdf$Replicate %in% replicate,]
    }
    if (type != "all") {
        inputdf <- inputdf[inputdf$Type == type,]
    }

    nrow(inputdf)
    if (!nrow(inputdf)) {
        stop("no matching data")
    }
    else {
        locations <- as.list(inputdf$Source)
        message("   ... matching data sets being loaded are:\n")
        write.table(inputdf,row.names=F,quote=F)
        message("")
    }

    sampleids <- as.list(inputdf$Sample)
    treatmentvec <- c()
    tmp <- rle(inputdf$Sample)
    for (i in 1:length(tmp$lengths)) {
        treatmentvec <- c(treatmentvec,rep(i-1,tmp$lengths[i]))
    }

    mrobj <- methRead(location = locations,
                      sample.id = sampleids,
                      assembly = assembly,
                      header = TRUE,
                      treatment = treatmentvec,
                      mincov = 1
                     )
    message("... done ..")
    return(mrobj)
}
