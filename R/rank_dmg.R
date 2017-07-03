#' rank_dmg()
#' This function sorts dmgenes by ADMpNucl and plots bar graph for demonstration 
#'
#' @param explore_dmsg_summaries A list returned by explore_dmsg()
#' @param outflabel A string to identify the study in the output file
#' 
#' @return A list of data frames summarizing genes sorted by ADMpNucl 
#' 
#' @importFrom utils write.table 
#' @importFrom ggplot2 ggplot
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                        type="CpGhsm",mincov=1,assembly="Amel-4.5")
#'   genome_ann <- get_genome_annotation(myfiles$parameters)
#'   dmsgList <- det_dmsg(AmHE,genome_ann,
#'                        threshold=25.0,qvalue=0.01,mc.cores=4,
#'                        destrand=TRUE,
#'                        outfile1="AmHE-dmsites.txt", 
#'                        outfile2="AmHE-dmgenes.txt")
#'   dmgprp <- show_dmsg(AmHE,dmsgList,min.nsites=10,destrand=TRUE,
#'                       outflabel="Am_HE")
#'   explore_dmsg_summaries <- explore_dmsg(AmHE,genome_ann,dmgprp,withglink="NCBIgene",
#'                             outflabel="Am_HE")
#'   rnk_pwcomparison_summaries <- rank_dmg(explore_dmsg_summaries,outflabel="Am_HE")
#'
#' @export
 
rank_dmg <- function(explore_dmsg_summaries,outflabel){
    message('... rank_dmg ...')
    ## read explore_dmsg_summaries
    #
    pw_summaries    <- explore_dmsg_summaries[[2]] 
    comparison_list <- names(pw_summaries)
    ## sort the gene ID by ADMpNucl 
    # and plot the distribution
    #
    rnk_summaries <- lapply(comparison_list, function(comparison){
        data     <- pw_summaries[[comparison]]
        wtoutfile <- paste("rnk_dmg",outflabel,sep="-")
        wtoutfile <- paste(wtoutfile,comparison,sep="_")
        wtoutfile <- paste(wtoutfile,"txt",sep=".")
        data <- data[order(- data$ADMpNucl),]
        write.table(data,wtoutfile,sep="\t",row.names=FALSE,quote=FALSE)

        ## plot the distribution of ADMpNucl ...
        #
        message('   ... sorted by ADMpNucl ...')
        pdf(paste("rnk_dmg_",comparison,'.pdf',sep=''))
        data$gene_ID <- factor(data$gene_ID,levels=data$gene_ID[order(data$ADMpNucl,decreasing=TRUE)])
        print(ggplot(data, aes(x=gene_ID,y=ADMpNucl)) + geom_col())
        print(ggplot(head(data,25), aes(x=gene_ID,y=ADMpNucl)) + geom_col())
        dev.off()

        return(data)
    })
    names(rnk_summaries) <- comparison_list
    return(rnk_summaries)
}
