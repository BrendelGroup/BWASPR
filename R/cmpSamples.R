#' cmpSamples()
#' This function will methylKit::unite() input samples and compare
#' them using methylKit::getCorrelation() and methylKit::PCASamples.
#'
#' @param mrobj A methylRaw object as returned by methylKit::unite()
#' @param destrand methylKit::unite() parameter; default: FALSE.
#'   destrand=TRUE combines CpG methylation calls from both strands.
#' @param mc.cores Integer denoting how many cores should be used for parallel
#'   diffential methylation calculations
#' @param plotfile If specified other than the default "", then plots
#'   are saved in PDF file "plotfile".pdf; otherwise either no plots are
#'   generated or the R default plot outputs are used.
#'
#' @return the data obtained from methylKit::unite()
#'
#' @importFrom methylKit unite getData getCorrelation PCASamples
#' @importFrom grDevices dev.off pdf
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                        sample="all",replicate=c(0),
#'                        type="CpGhsm", mincov=1,
#'                        assembly="Amel-4.5")
#'   cmpSamples(AmHE,destrand=TRUE,mc.cores=4,plotfile="myplots")
#'
#' @export

cmpSamples <- function(mrobj,destrand=FALSE,mc.cores=1,plotfile="") {
    message("... comparing samples ..")

    mbobj <- unite(mrobj,destrand=destrand,mc.cores=mc.cores)
    data <- getData(mbobj)

    if (plotfile != "") {
        pdf(paste(plotfile,"pdf",sep="."))
        getCorrelation(mbobj, plot=T)
        PCASamples(mbobj, adj.lim=c(2, 2))
        dev.off()
    }
    else {
        getCorrelation(mbobj, plot=F)
        PCASamples(mbobj, adj.lim=c(2, 2))
    }

    message("... done ..")
    return(data)
}
