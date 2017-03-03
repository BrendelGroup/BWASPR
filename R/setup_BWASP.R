#' setup_BWASPR()
#' Setting up a BWASPR analysis structure
#'
#' @param datafile A BWASP data description file consisting of rows of
#'   tab-deliminted entries with consecutive labesl for Species, Study,
#'   Sample, Replicate, Type, and BWASP mcalls file location.
#' @param parfile  A BWASPR parameter file with parameter settings used in
#'   in the analysis.
#'
#' @return A list consisting of a data frame with the data input files and a 
#'         data frame with parameters
#'
#' @examples
#'   setup_BWASRP(datefile="./Pc.dat",parfile="Pc.par")
#'
#' @export

setup_BWASPR <- function(datafile,parfile){
    message("... setting up BWASPR based on your data and parameter files ...")
    dfiles <- read.csv(datafile, sep="\t", skip=0, header=FALSE,
                       comment.char="#",
                       stringsAsFactors=FALSE
                      )
    colnames(dfiles) <- c("Species", "Study", "Sample", "Replicate", "Type",
                          "Source")
    pfile<- read.csv(parfile, sep="\t", skip=0, header=FALSE,
                     comment.char="#",
                     stringsAsFactors=FALSE
                    )
    colnames(pfile) <- c("Variable", "Value")
    files <- list(dfiles,pfile)
    message("... done ...")
    return(files)
}
