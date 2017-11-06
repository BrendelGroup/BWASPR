#' cmStats()
#' This function generates coverage and methylation statistics for the
#'   input data set.
#'
#' @param mrobj A methylRaw object
#' @param type The type should be CpGhsm or CpGscd, as per BWASP output.
#' @param covlist A vector of coverage threshols; default: c(10)
#' @param locount Low coverage threshold for statistics on a range of sites; default: c(100)
#' @param hicount High coverage threshold for statistics on a range of sites; default: c(1000)
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
#'   cmStats(AmHE[[1]],type="CpGhsm",covlist=c(10,20),locount=100,hicount=1000,
#'           outfile="Am_HE-statistics.txt",plotfile="myplots")
#'
#' @export

cmStats <- function(mrobj,type="CpGhsm",covlist=c(10),locount=100,hicount=1000,outfile="",plotfile="") {
    if (outfile != "") {
        sink(outfile)
    }
    sampledata <- getData(mrobj)
    sampledata$PrcntM <- 100.*sampledata$numCs/sampledata$coverage
    slabel <- mrobj@sample.id

    cat( sprintf(
      "Number of \"%s\" %s-sites with minimal and higher level coverage:\n\n",
      slabel,type) )
    cat( sprintf(
      "  number of \"%s\" %s-sites with coverage >= %2d: %6d\n",
      slabel, type, min(sampledata$coverage), length(sampledata$coverage)) )
    for (n in  covlist) {
      cn <- length(sampledata[["coverage"]][sampledata[["coverage"]] >= n])
      cat( sprintf(
        "  number of \"%s\" %s-sites with coverage >= %2d: %6d\n", slabel, type, n, cn) )
    }
    cat( sprintf( "\n") )


    if (plotfile!="") {
      pdf(paste(plotfile,"pdf",sep="."))
    }

    message("... calculating sample-wide ",type," coverage and methylation statistics for " ,slabel," ...")
    cat( sprintf(
      "\nCoverage and methylation statistics for \"%s\" %s-sites at different levels of minimum coverage:\n",
      slabel,type) )

    covlist <- c(min(sampledata$coverage),covlist)
    for (n in covlist) {
      mrobjFiltered <- filterByCoverage(mrobj,lo.count=n,lo.perc=NULL,hi.count=NULL,hi.perc=NULL)
      sampledata <- getData(mrobjFiltered)

      cat( sprintf( "\n\nmethylKit::getCoverageStats output for \"%s\" %s-sites (#: %d) at minimum coverage %d\t- ",
                    slabel,type,length(sampledata$coverage),n ) )
      getCoverageStats(mrobjFiltered,plot=F,both.strands=F)
      cat( sprintf( "methylKit::getMethylationStats output for \"%s\" %s-sites (#: %d) at minimum coverage %d\t- ",
                    slabel,type,length(sampledata$coverage),n ) )
      getMethylationStats(mrobjFiltered,plot=F,both.strands=F)

      if (plotfile!="") {
        subtitle <- paste(slabel,type,"with coverage at least",n,
                          " (number of sites: ",length(sampledata$coverage),")",sep=" ")
        getCoverageStats(mrobjFiltered,plot=T,sub=subtitle,both.strands=F,labels=FALSE,
                         cex.main=0.75, cex.sub=0.5, cex.axis=0.5, cex.lab=0.5, breaks=20)
        getMethylationStats(mrobjFiltered,plot=T,sub=subtitle,both.strands=F,labels=FALSE,
                            cex.main=0.75, cex.sub=0.5, cex.axis=0.5, cex.lab=0.5)
      }
    }

    mrobjFiltered <- filterByCoverage(mrobj,lo.count=locount,lo.perc=NULL,hi.count=hicount,hi.perc=NULL)
    sampledata <- getData(mrobjFiltered)

    cat( sprintf( "\n\nmethylKit::getCoverageStats output for \"%s\" %s-sites (#: %d) in coverage range [%d-%d]\t- ",
                  slabel,type,length(sampledata$coverage),locount,hicount ) )
    getCoverageStats(mrobjFiltered,plot=F,both.strands=F)
    cat( sprintf( "methylKit::getMethylationStats output for \"%s\" %s-sites (#: %d) in coverage range [%d-%d]\t- ",
                  slabel,type,length(sampledata$coverage),locount,hicount ) )
    getMethylationStats(mrobjFiltered,plot=F,both.strands=F)

    if (plotfile!="") {
      subtitle <- paste(slabel,type,"coverage range [",locount,"-",hicount,"]",
                        " (number of sites: ",length(sampledata$coverage),")",sep=" ")
      getCoverageStats(mrobjFiltered,plot=T,sub=subtitle,both.strands=F,labels=FALSE,
                       cex.main=0.75, cex.sub=0.5, cex.axis=0.5, cex.lab=0.5, breaks=20)
      getMethylationStats(mrobjFiltered,plot=T,sub=subtitle,both.strands=F,labels=FALSE,
                          cex.main=0.75, cex.sub=0.5, cex.axis=0.5, cex.lab=0.5)
    }

    if (plotfile!="") {
      dev.off()
    }
    if (outfile != "") {
        sink()
    }
    message("... done ...")
    return(cn)
}
