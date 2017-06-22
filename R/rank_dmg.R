#' rank_dmg() 
#' This function subsets the mrobj by the GRanges that the user provides 
#' and returns a List of dataframes containing the msites info summries and
#' save the dataframes as tab delimited files
#'
#' @param mrobj A methyRaw/methRawList object or a methyRawList object
#' @param region.gr A Granges object that the user provieds 
#' @param outflabel A string to identify the output dataframe 
#' 
#' @return A list of data frames that contains the msites info within the GRanges the user provided
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
#'   summaries <- subset_mrobj(AmHE,region.gr=genome_ann$promoter,
#'                             outflabel="Am_HE_promoter")
#'
#' @export

rank_dmg <- function(mrobj,region.gr,withglink="NCBIgene",
                         outflabel="") {
    message('... rank_dmg ...')
    message('... \'id\'& \'gene_name\' is required in region.gr ...')
    # read basic information
    #
    sample_list     <- getSampleID(mrobj)
    # for each sample, subset the msites info within the region.gr
    #
    message('   ... subset individual sample ...')
    sr_summaries <- lapply(sample_list, function(sample) {
        message(paste('      ... rank ',sample,' in interested region...',sep=''))
        # subset the mrobj
        #
        sites             <- reorganize(mrobj,
                                        sample.ids=list(sample),
                                        treatment=c(0))[[1]]
        sites.gr          <- as(sites,'GRanges')
        # calc the methylation level in percentage
        #
        sites.gr$perc_meth <- (sites.gr$numCs/sites.gr$coverage) * 100
        # annotate the msites with genes
        #
        match              <- findOverlaps(sites.gr,region.gr)
        sites.gr           <- sites.gr[queryHits(match)]
        region.gr          <- region.gr[subjectHits(match)]
        # combine
        sites.df            <- as.data.frame(sites.gr)
        region.df           <- as.data.frame(region.gr)
        colnames(region.df) <- lapply(colnames(region.df),
                                    function(i) paste('region',i,sep='_'))

        sites_region        <- cbind(sites.df,region.df)

        wtoutfile <- paste(outflabel,sample,"txt",sep=".")
        write.table(sites_region, wtoutfile, sep='\t',
                    row.names=FALSE, quote=FALSE)
      
        # calc a set of parameters for each gene
        #
        if (withglink == "NCBIgene"){
            ss_summary <- sites_region %>% group_by(region_ID) %>%
                summarise(rwidth = round(mean(region_width),2),
                          nbrsites = n(),
                          nbrper10kb = round((nbrsites/rwidth)*10000,2),
                                              pmsum = round(sum(perc_meth),2),
                          pmpersite = round(pmsum/nbrsites,2),
                          mpernucl = round(pmsum/rwidth,2),
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
        # order the regions by nbrper10kb 
        #
        ss_summary <- ss_summary[order(- ss_summary$nbrper10kb),]
        ss_summary <- subset(ss_summary, select = -c(pmsum))

        outfile <- paste("rnk",outflabel,sep="-")
        outfile <- paste(outfile,sample,sep="_")
        ptoutfile <- paste(outfile,"pdf",sep=".")
        wtoutfile <- paste(outfile,"txt",sep=".")
        write.table(ss_summary, wtoutfile, sep='\t',
                    row.names=FALSE, quote=FALSE)
        pdf(ptoutfile)
        ss_summary$region_ID <- factor(ss_summary$region_ID,
                                       levels=unique(as.character(ss_summary$region_ID)))
        print(ggplot(ss_summary, aes(x=region_ID,y=nbrper10kb)) + geom_col())
        print(ggplot(head(ss_summary,50), aes(x=region_ID,y=nbrper10kb)) + geom_col())

        dev.off()
        return(ss_summary)
     })
 
       
    message('... subset_mrobj finished ...')
    names(sr_summaries) = sample_list
    return(sr_summaries)
}
