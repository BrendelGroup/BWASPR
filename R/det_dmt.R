#' det_dmt()
#'   This function generates list of differentially methylated CpG tiles and the
#'   corresponding genes, as determined by application of methylKit functions.
#'
#' @param mrobj A methylKit methylRaw or methylRawList object
#' @param genome_ann Genome annotation returned by get_genome_annotation()
#' @param wsize The window size for the sliding windows (tiles), default wsize = 1000
#' @param stepsize The step size for the sliding window starts, default stepsize = 1000
#' @param threshold Cutoff for percent methylation difference, default threshold = 25.0
#' @param qvalue Cutoff for q-value, default q-value = 0.01
#' @param mc.cores Integer denoting how many cores should be used for parallel
#'   diffential methylation calculations
#' @param destrand methylKit::unite() parameter; default: FALSE.
#'   destrand=TRUE combines CpG methylation calls from both strands.
#' @param outfile1 File name to which differentially methylated tiles are written
#' @param outfile2 File name to which differentially methylated genes are written
#'
#' @return A Granges object that contains a list of genes that have diff methylated C tiles
#'
#' @import methylKit
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BiocGenerics as.data.frame
#' @importFrom S4Vectors mcols
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",type="CpGscd",
#'                        mincov=1,assembly="Amel-4.5")
#'   genome_ann <- get_genome_annotation(myfiles$parameters)
#'   mTlist <- det_dmt(AmHE,genome_ann,threshold=25.0,qvalue=0.01,mc.cores=4,
#'                     outfile1="AmHE-dmtiles.txt", 
#'                     outfile2="AmHE-dmgenes.txt") 
#'
#' @export

det_dmt <- function(mrobj,genome_ann,wsize=1000,stepsize=1000,
		     threshold=25.0,qvalue=0.01,mc.cores=1,
                     destrand=FALSE,
                     outfile1="dmt.txt",
                     outfile2="dmg.txt") {
    message('... det_dmt() ...')
    sample_list <- getSampleID(mrobj)
    treatment_list <- getTreatment(mrobj)
    genes <- genome_ann$gene
    nbtreatments <- length(unique(treatment_list))

    tiles <- tileMethylCounts(mrobj,win.size=wsize,step.size=stepsize)

    if (nbtreatments < 2) {
 	return(tiles)
    }

# If mrobj consists of multiple treatments, compare corresponding samples:
    unique_treatment_list <- unique(treatment_list)
    pairs <- combn(unique_treatment_list, 2, simplify=FALSE)
    dmtiles.gr <- lapply(pairs, function(pair) {
        pair_samples <- list(sample_list[pair[1]+1],sample_list[pair[2]+1])
        pair_treatments <- c(pair[1],pair[2])
        pair_tiles <- reorganize(tiles,
                                 sample.ids=pair_samples,
                                 treatment=pair_treatments
                                )
        pair_tiles@treatment=c(0,1)
        utiles <- unite(pair_tiles,destrand=destrand)
#getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
        pairname <- paste(sample_list[pair[1]+1],
                          sample_list[pair[2]+1],sep=".vs.")
        message(paste("... comparison: ",pairname," ...",sep=""))
        pairdiff <- calculateDiffMeth(utiles,mc.cores=mc.cores)
        diftiles <- getMethylDiff(pairdiff,difference=threshold,qvalue=qvalue)
        gr <- as(diftiles,"GRanges")
        if (length(gr) > 0) {
            S4Vectors::mcols(gr)$comparison <- pairname
        }
        message("... done ...")
        return(gr)
       }
       )
    if (length(dmtiles.gr) == 0) {
        message("No differentially methylated tiles are found.")
        return(list('dmgenes' = GRanges(),
                    'dmtiles' = GRanges()))
    }
    if (outfile1 != '') {
        dmtiles <- unlist(GRangesList(unlist(dmtiles.gr)))
        dmtiles$qvalue <- round(dmtiles$qvalue,3)
        dmtiles$meth.diff <- round(dmtiles$meth.diff,2)
        dmtiles.df <- BiocGenerics::as.data.frame(dmtiles)
        write.table(dmtiles.df,file=outfile1,sep='\t',row.names=FALSE,quote=FALSE)
    }

    dmgenes.gr <- lapply(seq_along(pairs), function(i) {
        gr <- subsetByOverlaps(genes,dmtiles.gr[[i]],ignore.strand=TRUE)
        pairname <- paste(sample_list[pairs[[i]][1]+1],
                          sample_list[pairs[[i]][2]+1],sep=".vs.")
        if (length(gr) > 0) {
            S4Vectors::mcols(gr)$comparison <- pairname
        }
        return(gr)
       }
       )
    if (length(dmgenes.gr) == 0) {
        message("No differentially methylated genes are found.")
        return(list('dmtiles' = dmtiles.gr,
		    'dmgenes' = GRanges()  ))
    }
    if (outfile2 != '') {
        dmgenes <- unlist(GRangesList(unlist(dmgenes.gr)))
        dmgenes.df <- BiocGenerics::as.data.frame(dmgenes)
        write.table(dmgenes.df,file=outfile2,sep='\t',row.names=FALSE,quote=FALSE)
    }

    message('... det_dmt() finished ...')
    return(list('tiles'   = tiles,
                'dmtiles' = dmtiles.gr,
                'dmgenes' = dmgenes.gr))
}
