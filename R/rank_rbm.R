#' rank_rbm() 
#'   This function subsets the provided MethylRaw objects by the specified GRange
#'   and returns a list of dataframes containing the GRange regions ranked by
#'   methylation site density.
#'
#' @param mrobj A methyRaw/methRawList object
#' @param region.gr A GRanges object (typically based on the relevant genome annotation)
#' @param rlabel A string to identify the type of region analyzed
#' @param withglink Either NCBIgene or "" to indicate inculsion of a gene link column in the output
#' @param outflabel A string to identify the output dataframe 
#' 
#' @return A list of dataframes containing the GRange regions ranked by methylation
#'   site density
#' 
#' @importFrom GenomicRanges findOverlaps
#' @importFrom utils write.table 
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom dplyr group_by %>% summarise
#' @importFrom ggplot2 ggplot aes geom_col
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                        type="CpGhsm",mincov=1,assembly="Amel-4.5")
#'   genome_ann <- get_genome_annotation(myfiles$parameters)
#'   summaries <- rank_rbm(AmHE,region.gr=genome_ann$gene,rlabel="sir",
#'                         withglink="NCBIgene",outflabel="Am_HE_gene")
#'
#' @export

rank_rbm <- function(mrobj,region.gr,rlabel="sir",withglink="",outflabel="") {
    message('... rank_rbm ...')
    # read basic information
    #
    sample_list <- getSampleID(mrobj)
    message('   ... subset individual samples ...')
    rnk_summaries <- lapply(sample_list, function(sample) {
        message(paste('      ... rank ',sample,' regions of interest ...',sep=''))
        # subset the mrobj ...
        #
        sites             <- reorganize(mrobj,
                                        sample.ids=list(sample),
                                        treatment=c(0))[[1]]
        sites.gr          <- as(sites,'GRanges')
        sites.gr$perc_meth <- (sites.gr$numCs/sites.gr$coverage) * 100

        # identify the msites within the specified regions ...
        #
        match              <- suppressWarnings(findOverlaps(sites.gr,region.gr,ignore.strand=TRUE))
        sites.gr           <- sites.gr[queryHits(match)]
        region.gr          <- region.gr[subjectHits(match)]
        sites.df            <- as.data.frame(sites.gr)
        region.df           <- as.data.frame(region.gr)
        colnames(region.df) <- lapply(colnames(region.df),
                                      function(i) paste('region',i,sep='_'))
        sites_region        <- cbind(sites.df,region.df)

        wtoutfile <- paste(rlabel,outflabel,sep="-")
        wtoutfile <- paste(wtoutfile,sample,sep="_")
        wtoutfile <- paste(wtoutfile,"txt",sep=".")
        write.table(sites_region, wtoutfile, sep='\t', row.names=FALSE, quote=FALSE)
      
        # calculate site densities for each region ...
        #
        if (withglink == "NCBIgene"){
            ss_summary <- sites_region %>% group_by(region_ID) %>%
                summarise(rwidth = round(mean(region_width),2),
                          nbrsites = n(),
                          nbrper10kb = round((nbrsites/rwidth)*10000,2),
                          pmsum = round(sum(perc_meth),2),
                          pmpersite = round(pmsum/nbrsites,2),
                          pmpernucl = round(pmsum/rwidth,2),
                          pglink = paste("https://www.ncbi.nlm.nih.gov/gene/?term",
                                         unique(region_Name),sep="=")
                          )
        }
        else {
            ss_summary<- sites_region %>% group_by(region_ID) %>%
                summarise(rwidth = round(mean(region_width),2),
                          nbrsites = n(),
                          nbrper10kb = round((nbrsites/rwidth)*10000,2),
                          pmsum = round(sum(perc_meth),2),
                          pmpersite = round(pmsum/nbrsites,2),
                          pmpernucl = round(pmsum/rwidth,2)
                          )
        }
        # order the regions by nbrper10kb ...
        #
        ss_summary <- ss_summary[order(- ss_summary$nbrper10kb),]
        ss_summary <- subset(ss_summary, select = -c(pmsum))

        outfile <- paste("rnk",rlabel,outflabel,sep="-")
        outfile <- paste(outfile,sample,sep="_")
        wtoutfile <- paste(outfile,"txt",sep=".")
        write.table(ss_summary, wtoutfile, sep='\t',
                    row.names=FALSE, quote=FALSE)

        ptoutfile <- paste(outfile,"pdf",sep=".")
        pdf(ptoutfile)
        ss_summary$region_ID <- factor(ss_summary$region_ID,
                                       levels=unique(as.character(ss_summary$region_ID)))
        print(ggplot(ss_summary, aes(x=region_ID,y=nbrper10kb)) + geom_col())
        print(ggplot(head(ss_summary,25), aes(x=region_ID,y=nbrper10kb)) + geom_col())

        dev.off()
        return(ss_summary)
     })

    message('... rank_rbm finished ...')
    names(rnk_summaries) = sample_list
    return(rnk_summaries)
}
