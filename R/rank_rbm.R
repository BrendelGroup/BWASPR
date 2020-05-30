#' rank_rbm() 
#'   This function subsets the provided methylRaw objects by the specified GRange
#'   and returns a list of dataframes containing the GRange regions ranked by
#'   methylation site density.
#'
#' @param mrobjscd A methylKit methylRaw or methylRawList object for scd sites
#' @param mrobjhsm A methylKit methylRaw or methylRawList object for hsm sites
#' @param region.gr A GRanges object (typically based on the relevant genome annotation)
#' @param rlabel A string to identify the type of region analyzed
#' @param withglink Either NCBIgene or "" to indicate inclusion of a gene link column
#'   in the output if available
#' @param outflabel A string to identify the output dataframe 
#' 
#' @return A list of lists of dataframes containing (for each sample) the GRange
#'   regions ranked by overall methylation percentage (over all scd sites) or
#'   methylation hsm site density
#' 
#' @importFrom GenomicRanges findOverlaps
#' @importFrom utils write.table 
#' @importFrom S4Vectors subjectHits queryHits
#' @import     dplyr
#' @importFrom ggplot2 ggplot aes geom_col
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHEscd <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                           type="CpGscd",mincov=1,assembly="Amel-4.5")
#'   AmHEhsm <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                           type="CpGhsm",mincov=1,assembly="Amel-4.5")
#'   genome_ann <- get_genome_annotation(myfiles$parameters)
#'   rgL <- rank_rbm(AmHEscd,AmHEhsm,region.gr=genome_ann$gene,rlabel="genes",
#'                   withglink="NCBIgene",outflabel="Am_HE")
#'
#' @export

rank_rbm <- function(mrobjscd,mrobjhsm,region.gr,rlabel="",withglink="",outflabel="") {
    message('... rank_rbm ...')
    # read basic information:
    #
    sample_list <- getSampleID(mrobjscd)

    # get regional counts of methylation calls:
    #
    rcL <- regionCounts(mrobjscd,region.gr)

    # determinine methylation (percent of calls) in specified regions ...
    #
    message('   ... subset individual samples ...')
    #
    MpRegion <- lapply(1:length(getSampleID(mrobjscd)), function(i) {
        rcd           <- getData(rcL[[i]])
        rcd           <- cbind(rcd,round(100.*rcd[,6]/rcd[,5],2))
	names(rcd)[8] <- "prcntM"
	rcg   <- makeGRangesFromDataFrame(rcd,keep.extra.columns=TRUE)
	match <- findOverlaps(rcg,region.gr,type="equal")
	r     <- as.data.frame(rcg[queryHits(match)])
	a     <- as.data.frame(region.gr[subjectHits(match)])
        colnames(a) <- lapply(colnames(a), function(i) paste('region',i,sep='_'))
	ra    <- cbind(r,a)
	c     <- c("region_ID","seqnames","start","end","width","strand",
	           "coverage","numCs","numTs","prcntM")
	ra    <- ra[c]
        wtoutfile <- paste("ranked",rlabel,"byPrcntM",outflabel,sep="-")
        wtoutfile <- paste(wtoutfile,sample_list[i],sep="_")
        wtoutfile <- paste(wtoutfile,"txt",sep=".")
        write.table(ra[order(- ra$prcntM),], wtoutfile, sep='\t', row.names=FALSE, quote=FALSE)
	return(ra)
    })

    # determinine methylation site density in specified regions ...
    #
    SpRegion <- lapply(sample_list, function(sample) {
        message(paste('      ... rank ',sample,' regions of interest ...',sep=''))
        # subset mrobjhsm ...
        #
        sites               <- reorganize(mrobjhsm,
                                          sample.ids=list(sample),
                                          treatment=c(0)
					 )[[1]]
        sites.gr            <- as(sites,'GRanges')
        sites.gr$perc_meth  <- (sites.gr$numCs/sites.gr$coverage) * 100

        # identify the msites within the specified regions ...
        #
        match               <- suppressWarnings(findOverlaps(sites.gr,region.gr,ignore.strand=TRUE))
        sites.gr            <- sites.gr[queryHits(match)]
        region.gr           <- region.gr[subjectHits(match)]
        sites.df            <- as.data.frame(sites.gr)
        region.df           <- as.data.frame(region.gr)
        colnames(region.df) <- lapply(colnames(region.df),
                                      function(i) paste('region',i,sep='_'))
        sites_region        <- cbind(sites.df,region.df)

        wtoutfile <- paste("sites-in",rlabel,outflabel,sep="-")
        wtoutfile <- paste(wtoutfile,sample,sep="_")
        wtoutfile <- paste(wtoutfile,"txt",sep=".")
        write.table(sites_region, wtoutfile, sep='\t', row.names=FALSE, quote=FALSE)
      
        # calculate site densities for each region ...
        #
        if (withglink == "NCBIgene"){
            rSstats <- sites_region %>% group_by(region_ID) %>%
                summarise(rwidth = round(mean(region_width),2),
                          nbrsites = dplyr::n(),
                          nbrper10kb = round((nbrsites/rwidth)*10000,2),
                          pmsum = round(sum(perc_meth),2),
                          pmpersite = round(pmsum/nbrsites,2),
                          pmpernucl = round(pmsum/rwidth,2),
                          pglink = paste("https://www.ncbi.nlm.nih.gov/gene/?term",
                                         unique(region_Name),sep="=")
                         )
        }
        else {
            rSstats <- sites_region %>% group_by(region_ID) %>%
                summarise(rwidth = round(mean(region_width),2),
                          nbrsites = dplyr::n(),
                          nbrper10kb = round((nbrsites/rwidth)*10000,2),
                          pmsum = round(sum(perc_meth),2),
                          pmpersite = round(pmsum/nbrsites,2),
                          pmpernucl = round(pmsum/rwidth,2)
                         )
        }
        # order the regions by nbrper10kb ...
        #
        rSstats <- rSstats[order(- rSstats$nbrper10kb),]
        rSstats <- subset(rSstats, select = -c(pmsum))

	# prepare the regional Site & Methylation statistics data frame:
	#
 	rSMdf <- merge(rSstats,MpRegion[match(sample,sample_list)],by="region_ID")
	rSMdf <- rSMdf[order(- rSMdf$nbrper10kb),]
	rSMdf$width <- NULL	# ... remove the duplicated width information

        outfile   <- paste("ranked",rlabel,"bySiteDensity",outflabel,sep="-")
        outfile   <- paste(outfile,sample,sep="_")
        wtoutfile <- paste(outfile,"txt",sep=".")
        write.table(rSMdf, wtoutfile, sep='\t',
                    row.names=FALSE, quote=FALSE)

        outfile   <- paste("plot",rlabel,"prcntM-vs-SiteDensity",outflabel,sep="-")
        outfile   <- paste(outfile,sample,sep="_")
        ptoutfile <- paste(outfile,"pdf",sep=".")
        pdf(ptoutfile)
        print(ggplot(head(rSMdf[rSMdf$nbrper10kb<200,],50), aes(x=nbrper10kb,y=prcntM)) +
	      labs(title=outfile) +
	      geom_point(size=2,shape=23,color="red") +
	      geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95)
             )

        dev.off()
        return(rSstats)
     })

    message('... rank_rbm finished ...')
    names(SpRegion) = sample_list
    names(MpRegion) = sample_list
    return(list(MpRegion,SpRegion))
}
