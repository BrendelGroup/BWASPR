#' get_diff_genes()
#' This function will get a list of genes with diff methylated CpG sites
#'
#' @param mrobj A methylRaw object or a methylRawList object.
#' @param genome_ann Genome annotation returned by get_genome_annotation()
#' @param threshold cutoff for percent methylation difference, default threshold = 0.25
#' @param qvalue cutoff for q-value, defualt q-value = 0.05
#'
#' @return A Granges object that contains a list of genes that have diff methylated C sites
#'
#' @import methylKit
#' @importFrom GenomicRanges GRangesList
#' @importFrom IRanges subsetByOverlaps
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",type="CpGhsm",
#'                        mincov=1,assembly="Amel-4.5")
#'   genome_ann <- get_genome_annotation(myfiles$parameters)
#'   meth_diff <- get_diff_genes(AmHE,genome_ann,threshold=0.25,qvalue=0.05)
#'
#' @export

get_diff_genes <- function(mrobj,genome_ann,threshold=0.25,qvalue=0.05){
    sample_list <- getSampleID(mrobj)
    treatment_list <- getTreatment(mrobj)
    genes <- genome_ann$gene
    nbtreatments <- length(unique(treatment_list))
    unique_treatment_list <- unique(treatment_list)
    # if number of treatment == 1, compare intra a caste
    if (nbtreatments == 1) {
        mrobj <- reorganize(mrobj,treatment=0 : (length(treatment_list) - 1))
        meth <- unite(mrobj, destrand = TRUE, mc.cores = 8)
        diff <- calculateDiffMeth(meth, mc.cores = 8)
        diff_th <- getMethylDiff(diff, difference = threshold, qvalue = qvalue)
        methygenes <- subsetByOverlaps(gene.gr, as(diff_th, 'GRanges'))
    }

    # if number of treatment > 1, compare inter castes
    if (nbtreatments > 1) {
        pairs <- combn(unique_treatment_list, 2, simplify=FALSE)
        for (pair in pairs) {
            pair_sample_list <- list(sample_list[pair[1]+1],sample_list[pair[2]+1])
            pair_treatment_list <- c(pair[1],pair[2])
            pair_mrobj <- reorganize(mrobj,
                                     sample.ids=pair_sample_list,
                                     treatment=c(pair[1],pair[2])
				    )
            meth <- unite(pair_mrobj,destrand=TRUE)
            meth@treatment=c(0,1)
            diff <- calculateDiffMeth(meth)
            assign(paste(pair, collapse = ''),
                   getMethylDiff(diff,difference=threshold,qvalue=qvalue)
                  )
        }

        all_diff_combination_result_names <- lapply(pairs, function(i) paste(i, collapse = ''))
        all_diff_sites <- sapply(all_diff_combination_result_names, function(i) as(get(i),'GRanges'))

        methygenes <- sapply(all_diff_combination_result_names, function(i) subsetByOverlaps(genes, as(get(i), 'GRanges')))
        diff_genes <- unlist(GRangesList(unlist(methygenes)))
    }
    return(list('diff_genes' = diff_genes,
                'diff_sites' = all_diff_sites))
}
