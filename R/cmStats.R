#' cmStats()
#' This function generates coverage and methylation statistics for the
#'   input data set.
#'
#' @param mrobj A methylRaw object
#' @param covlist a vector of coverage threshols; default: c(10)
#' @param outfile If specified then output is printed to the specified file.
#' @param plotfile If specified other than the default "", then plots
#'   are saved in PDF file "plotfile".pdf; otherwise either no plots are
#'   generated or the R default plot outputs are used.
#'
#' @return the coverage count for the highest threshold in covlist
#'
#' @importFrom methylKit getData getCoverageStats getMethylationStats
#  @importFrom grDevices dev.off pdf
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                        type="CpGhsm", mincov=1,assembly="Amel-4.5")
#'   cmStats(AmHE[[1]],covlist=c(10,20),outfile="Am_HE-statistics.txt",plotfile="myplots")
#'
#' @export

cmStats <- function(mrobj,covlist=c(10),outfile="",plotfile="") {
    if (outfile != "") {
        sink(outfile)
    }
    sampledata <- getData(mrobj)
    slabel <- mrobj@sample.id
    message("... calculating sample-wide coverage and methylation statistics for" ,slabel," ...")
    sampledata$PrcntM <- 100.*sampledata$numCs/sampledata$coverage

    cat( sprintf( "methylKit output for \"%s\"\t- ", slabel ) )
    getCoverageStats(mrobj,plot=F,both.strands=F)
    cat( sprintf( "methylKit output for \"%s\"\t- ", slabel ) )
    getMethylationStats(mrobj,plot=F,both.strands=F)

    if (plotfile!="") {
      pdf(paste(plotfile,"pdf",sep="."))
      getCoverageStats(mrobj,plot=T,both.strands=F,labels=FALSE,
                       cex.main=0.75, cex.sub=0.5, cex.axis=0.5, cex.lab=0.5)
      getMethylationStats(mrobj,plot=T,both.strands=F,labels=FALSE,
                          cex.main=0.75, cex.sub=0.5, cex.axis=0.5, cex.lab=0.5)
      dev.off()
    }

    cat( sprintf(
      "Number of \"%s\" sites with minimal and higher level coverage\n",
      slabel) )
    cat( sprintf(
      "  number of \"%s\" sites with coverage >= %2d: %6d\n",
      slabel, min(sampledata$coverage), length(sampledata$coverage)) )
    for (n in  covlist) {
      cn <- length(sampledata[["coverage"]][sampledata[["coverage"]] >= n])
      cat( sprintf(
        "  number of \"%s\" sites with coverage >= %2d: %6d\n", slabel, n,cn) )
    }
    cat( sprintf( "\n\n") )

    if (outfile != "") {
        sink()
    }
    message("... done ...")
    return(cn)
}
