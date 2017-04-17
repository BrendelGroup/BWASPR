#' show_dmsg() 
#'   This function generates will plot heatmaps of methylation level of sites in dmgenes
#'
#' @param mrobj A methyRaw object or a methyRawList object
#' @param dmsg A list containing GRanges objects dmsites and dmgenes returned by
#'   det_dmsg()
#' @param min.nsites Minimal number of msites per genes
#' 
#' @return A data frame
#' 
#' @importFrom methylKit percMethylation reorganize unite
#' @importFrom GenomicRanges findOverlaps values
#' @importFrom gplots heatmap.2 greenred
#' @importFrom utils write.table 
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom dplyr group_by_ %>%
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                        type="CpGhsm",mincov=1,assembly="Amel-4.5")
#'   genome_ann <- get_genome_annotation(myfiles$parameters)
#'   meth_diff <- det_dmsg(AmHE,genome_ann,
#'                         threshold=25.0,qvalue=0.01,mc.cores=4,
#'                         outfile1="AmHE-dmsites.txt", 
#'                         outfile2="AmHE-dmgenes.txt")
#'   show_dmsg <- show_dmsg(AmHE,meth_diff)
#'
#' @export

show_dmsg <- function(mrobj,dmsg,min.nsites=2){
    message('... show_dmsg() ...')
    # load dmsites and dmgenes and sample_match_list
    dmsites.gr          <- do.call("c", dmsg$dmsites)
    dmgenes.gr          <- do.call("c", dmsg$dmgenes)
    sample_match_list   <- as.list(unique(as.character(dmgenes.gr$comparison)))
    # analyze each sample_match    
    show_dmsg.df <- lapply(sample_match_list, function(sample_match){
        sample1         <- unlist(strsplit(sample_match,'\\.'))[1]
        sample2         <- unlist(strsplit(sample_match,'\\.'))[3]
        message(paste('... comparing ',sample1,' & ',sample2,' ...',sep=''))
        # subset the dmsites.gr & dmgenes.gr with this sample_match
        #
        pair_dmsites.gr <- dmsites.gr[GenomicRanges::values(dmsites.gr)$comparison%in%sample_match]
        pair_dmgenes.gr <- dmgenes.gr[GenomicRanges::values(dmgenes.gr)$comparison%in%sample_match]
        # subset the mrobj with current sample_match 
        #
        pair_mrobj      <- reorganize(mrobj,sample.ids=list(sample1,sample2),
                                      treatment=c(0,1))
        pair_meth       <- unite(pair_mrobj,destrand=TRUE)
        # calc methylation level
        #
        p_meth          <- round(percMethylation(pair_meth,rowids=FALSE,
                                                 save.txt=FALSE),2)
        pair_p_meth     <- cbind(pair_meth,p_meth)
        pair_p_meth.gr  <- as(pair_p_meth,'GRanges')
        # identify scd sites in each gene
        #
        match                 <- findOverlaps(pair_dmgenes.gr,pair_p_meth.gr)
        sub_pair_p_meth.gr    <- pair_p_meth.gr[subjectHits(match)]
        sub_pair_dmgenes.gr   <- pair_dmgenes.gr[queryHits(match)]
        # identify dmsites in scd sites
        #
        match2                <- findOverlaps(sub_pair_p_meth.gr,pair_dmsites.gr)
        pair_dmsites_index    <- queryHits(match2)
        # transform GRanges objects to dataframes and combine them
        #
        sub_pair_p_meth            <- as.data.frame(sub_pair_p_meth.gr)
        sub_pair_dmgenes           <- as.data.frame(sub_pair_dmgenes.gr)
        colnames(sub_pair_dmgenes) <- lapply(colnames(sub_pair_dmgenes),
                                             function(i) paste('gene',i,sep='_'))
        meth_dmg_comb              <- cbind(sub_pair_p_meth,
                                            sub_pair_dmgenes)
        # label each scd if it is a dmsite
        #
        meth_dmg_comb['is.dm']                    <- FALSE
        meth_dmg_comb[pair_dmsites_index,'is.dm'] <- TRUE
        # save
        #
        meth_dmg_comb <- meth_dmg_comb[colSums(! is.na(meth_dmg_comb))>0]
        write.table(meth_dmg_comb,
                    file=(paste(sample_match,'show_dmsg.txt',sep='_')),
                    sep='\t',row.names=FALSE,quote=FALSE)
        # split the dataframe
        #
        splitter     <- c('gene_ID','gene_Name','gene_gene')
        splitter     <- splitter[splitter%in%names(meth_dmg_comb)][1]
        grouped      <- meth_dmg_comb%>%group_by_(.dots=splitter)
        out          <- split(grouped,grouped[splitter])
        # plot heatmap for each dmgene
        #
        pdf(paste(sample_match,'.pdf',sep=''))
        for (i in out) {
            if (dim(i)[1] >= min.nsites) {
                plot <- as.matrix(i[,c(sample1,sample2)])
                # make sure that there are difference to show in the heatmap:
                if (! all(plot[1] == plot)) {
                    heatmap.2(plot, 
                              margins=c(10,10),
                              dendrogram='none',
                              Rowv=FALSE,
                              col=greenred(10),
                              trace='none',
                              main=paste('msites@gene',unique(i[splitter]),sep=''),
                              srtCol=45,
                              RowSideColors=as.character(as.numeric(i$is.dm)))
                }
            }
        }
        dev.off()
        return(meth_dmg_comb)
    })
    names(show_dmsg.df) <- sample_match_list
    message('... show_dmsg() finished ...')
    return(show_dmsg.df)
}
