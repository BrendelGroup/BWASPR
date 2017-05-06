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
#' @param mincov Integer value specifying the mincov parameter for readMeth();
#'   default: 1.
#' @param assembly Label for the underlying genome assembly version; default:
#'   "unknown".
#'
#' @return A methylKit-package methylRaw or methylRawList object (NULL if no
#'   matching data are found).
#'
#' @importFrom methylKit methRead
## @importFrom utils write.table
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",type="CpGhsm",
#'                mincov=1,assembly="Amel-4.5")
#'
#' @export

mcalls2mkobj <- function(inputdf,species="all",study="all",sample="all",
                         replicate=c(1:20),type="CpGhsm",mincov=1,
                         assembly="unknown") {
    message("... reading *.mcalls files into methylRaw objects ...")

    if (species != "all") {
        inputdf <- inputdf[inputdf$Species == species,]
    }
    if (!identical(study,"all")) {
        inputdf <- inputdf[inputdf$Study %in% study,]
    }
    if (!identical(sample,"all")) {
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
        message("NOTE: There are no  matching data for selection ",
		species, " ", study, " ", sample, " replicates=", replicate)
        return(NULL)
    }
    else {
        locations <- as.list(inputdf$Source)
        message("   ... matching data sets being loaded are:\n")
        write.table(inputdf,file=stderr(),row.names=F,quote=F)
        message("")
    }

# ... add replicate numbers to replicates in sampleids:
    s <- as.list(inputdf$Sample)
    r <- as.list(inputdf$Replicate)
    s <- mapply(function(x,y) if (y > 0) {paste(x,y,sep="_")} else {x}, s,r)
    sampleids <- as.list(s)

# ... group distinct samples as "treatments" in methylKit notation:
    treatmentvec <- c()
    tmp <- rle(inputdf$Sample)
    for (i in 1:length(tmp$lengths)) {
        treatmentvec <- c(treatmentvec,rep(i-1,tmp$lengths[i]))
    }

# ... add Study to sampleids if there are multiple studies:
    if (length(unique(inputdf$Study))>1) {
        s <- as.list(inputdf$Study)
        s <- mapply(function(x,y) paste(x,y,sep="_"), sampleids,s)
        sampleids <- as.list(s)
    }

    mrobj <- methRead(location = locations,
                      sample.id = sampleids,
                      assembly = assembly,
                      header = TRUE,
                      treatment = treatmentvec,
                      mincov = mincov
                     )
    message("... done ...")
    return(mrobj)
}
