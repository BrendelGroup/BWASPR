#' Setting up a BWASPR analysis structure
#'
#' @param datafile
#' mcalls2mkobj()
#' @return a methylBase object
#'
#' @importFrom methylKit methRead
#'
#' @export

mcalls2mkobj <- function(inputdf,species,study,type,assembly){
    message("... loading mc objects ..")
    if (species != "all") {
        inputdf <- inputdf[inputdf$Species == species,]
    }
    if (study != "all") {
        inputdf <- inputdf[inputdf$Study == study,]
    }
    if (type != "all") {
        inputdf <- inputdf[inputdf$Type == type,]
    }

    locations <- as.list(inputdf$Source)

    sampleids <- as.list(inputdf$Sample)
    treatmentvec <- c()
    tmp <- rle(inputdf$Sample)
    for (i in 1:length(tmp$lengths)) {
        treatmentvec <- c(treatmentvec,rep(i-1,tmp$lengths[i]))
    }

    mcobj <- methRead(location = locations,
                      sample.id = sampleids,
                      assembly = assembly,
                      header = TRUE,
                      treatment = treatmentvec,
                      mincov = 1
                      )
    return(mcobj)
}
