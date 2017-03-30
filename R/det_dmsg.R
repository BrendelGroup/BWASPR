#' det_dmsg()
#' This function will get a list of genes with diff methylated CpG sites
#'
#' @param mrobj A methylRaw object or a methylRawList object.
#' @param genome_ann Genome annotation returned by get_genome_annotation()
#' @param threshold cutoff for percent methylation difference, default threshold = 0.25
#' @param qvalue cutoff for q-value, defualt q-value = 0.05
#' @param mc.cores Integer denoting how many cores should be used for parallel
#'   diffential methylation calculations
#' @param outfile1 File name to which diff sites are written
#' @param outfile2 File name to which diff genes are written
#'
#' @return A Granges object that contains a list of genes that have diff methylated C sites
#'
#' @import methylKit
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",type="CpGhsm",
#'                        mincov=1,assembly="Amel-4.5")
#'   genome_ann <- get_genome_annotation(myfiles$parameters)
#'   meth_diff <- det_dmsg(AmHE,genome_ann,threshold=0.25,qvalue=0.05,mc.cores=4,
#'                         outfile1="AmHE-dmsites.txt", 
#'                         outfile2="AmHE-dmgenes.txt") 
#'
#' @export

det_dmsg <- function(mrobj,genome_ann,threshold=0.25,qvalue=0.05,mc.cores=1,
                     outfile1="methyl_dmsites.txt",
                     outfile2="methyl_dmgenes.txt"){
    sample_list <- getSampleID(mrobj)
    treatment_list <- getTreatment(mrobj)
    genes <- genome_ann$gene
    nbtreatments <- length(unique(treatment_list))
    unique_treatment_list <- unique(treatment_list)
# If mrobj consists of only one treatment, compare replicates:
    if (nbtreatments == 1) {
        mrobj <- reorganize(mrobj,treatment=0 : (length(treatment_list) - 1))
        meth <- unite(mrobj,destrand=TRUE,mc.cores=mc.cores)
        diff <- calculateDiffMeth(meth,mc.cores=mc.cores)
        diff_th <- getMethylDiff(diff,difference=threshold,qvalue=qvalue)
        dmsites.gr <- as(diff_th,'GRanges')
        dmgenes.gr <- subsetByOverlaps(gene.gr, as(diff_th, 'GRanges'))
    }

# If mrobj consists of multiple treatments, compare corresponding samples:
    if (nbtreatments > 1) {
        pairs <- combn(unique_treatment_list, 2, simplify=FALSE)
        for (pair in pairs) {
            pair_samples <- list(sample_list[pair[1]+1],sample_list[pair[2]+1])
            pair_treatments <- c(pair[1],pair[2])
            pair_mrobj <- reorganize(mrobj,
                                     sample.ids=pair_samples,
                                     treatment=pair_treatments
                                    )
            meth <- unite(pair_mrobj,destrand=TRUE)
            meth@treatment=c(0,1)
            pairdiff <- calculateDiffMeth(meth,mc.cores=mc.cores)
            assign(paste(sample_list[pair[1]+1],sample_list[pair[2]+1],sep=".vs."),
                   getMethylDiff(pairdiff,difference=threshold,qvalue=qvalue)
                  )
        }

        pairnames <- lapply(pairs, function(p) paste(sample_list[p[1]+1],sample_list[p[2]+1],sep=".vs."))
        dmsites.gr <- sapply(pairnames, function(p) {
                               gr <- as(get(p),'GRanges')
                               S4Vectors::mcols(gr)$comparison <- p
                               return(gr)
                               }
                            )
        dmsites <- unlist(GRangesList(unlist(dmsites.gr)))
        dmgenes.gr <- sapply(pairnames, function(p) {
                               gr <- subsetByOverlaps(genes,
                                                      as(get(p),'GRanges')
                                                     )
                               S4Vectors::mcols(gr)$comparison <- p
                               return(gr)
                               }
                            )
        dmgenes <- unlist(GRangesList(unlist(dmgenes.gr)))
    }
    if (outfile1 != ''){
        dmsites$qvalue <- round(dmsites$qvalue,3)
        dmsites$meth.diff <- round(dmsites$meth.diff,2)
        write.table(dmsites,file=outfile1,sep='\t',row.names=FALSE,quote=FALSE)
    }
    if (outfile2 != ''){
        write.table(dmgenes,file=outfile2,sep='\t',row.names=FALSE,quote=FALSE)
    }
    return(list('dmgenes' = dmgenes.gr,
                'dmsites' = dmsites.gr))
}
