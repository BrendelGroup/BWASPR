#' Setting up a BWASPR analysis structure
#'
#' @param datafile, parfile
#' setup_BWASPR()
#' @return a list consisting of a data frame with the data input files) and a 
#'         data frame with parameters
#' @export

setup_BWASPR <- function(datafile,parfile){
    message("Setting up")
    dfiles <- read.csv(datafile, sep="\t", skip=0, header=FALSE,
                       comment.char="#",
                       stringsAsFactors=FALSE
                      )
    colnames(dfiles) <- c("Species", "Study", "Sample", "Replicate", "Type", "Source")
    pfile<- read.csv(parfile, sep="\t", skip=0, header=FALSE,
                     comment.char="#",
                     stringsAsFactors=FALSE
                    )
    colnames(pfile) <- c("Variable", "Value")
    files <- list(dfiles,pfile)
    return(files)
}
