#' get_diff_genes()
#' ??This function will get a list of genes with diff methylated CpG sites
#'
#' @param mrobj A methylRaw object or a methylRawList object.
#' @param threshold cutoff for percent methylation difference
#' @param qvalue cutoff for q-value
#'
#' @return A Granges object that contains a list of genes that have diff methylated C sites
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",type="CpGhsm",
#'                        mincvrg=1,assembly="Amel-4.5")
#'   genome <- read_genome_info(myfiles$parameters)
#'   meth_diff <- get_diff_genes(AmHE, genome)
#'
#' @export

get_diff_genes<- function(mrobj, genome, threshold = 0, qvalue = 1){
	# fetch sample_list and treatment_list
	sample_list <- getSampleID(mrobj)
	treatment_list <- getTreatment(mrobj)
	# get the gene GRanges
	genes <- genome[['gene']]
	# get the number of treatments
	# if number of treatment == 1, compare intra a casteqaz123
	# if number of treatment > 1, compare inter castes
	number_of_treatment <- length(unique(treatment_list))
	unique_treatment_list <- unique(treatment_list)
	# if number of treatment == 1, compare intra a caste
	if (number_of_treatment == 1) {
		# change the treatment list to a c(0:length(treatment_list))
		mrobj <- reorganize(mrobj, treatment = 0 : (length(treatment_list) - 1))
		# unite the mrobj
		meth <- unite(mrobj, destrand = TRUE, mc.cores = 8)
		# calculate the diff meth
		diff <- calculateDiffMeth(meth, mc.cores = 8)
		# set the threshold
		diff_th <- getMethylDiff(diff, difference = threshold, qvalue = qvalue)
		# get the genes with diff meth sites
		methygenes <- subsetByOverlaps(gene.gr, as(diff_th, 'GRanges'))
	}

	# if number of treatment > 1, compare inter castes
	if (number_of_treatment > 1) {
		# generate combinations between different castes
		pair_combination <- combn(unique_treatment_list, 2)
		pair_combination <- lapply(seq_len(ncol(pair_combination)), function(i) pair_combination[, i])
		# go through each pair in all the pair combinations
		for (pair in pair_combination) {
			# get the sample ID list and treatment list from the 2 castes for comparison
		  # very clumsy here, need modification
		  # modified with grepl
		  pair_mask <- grepl(paste(pair, collapse = '|'), treatment_list)
			pair_sample_list <- sample_list[pair_mask]
			pair_treatment_list <- treatment_list[pair_mask]
			# subset the mrobj that contains all the mcfiles
			# to create a mrobj that only contains the mcfiles from the 2 castes for comparision
			pair_mrobj <- reorganize(mrobj,
									 sample.ids = pair_sample_list,
									 treatment = pair_treatment_list)
			# unite the mrobj in the mrobjRawList for calculate the meth diff
			meth <- unite(mrobj, destrand = TRUE)
			# calculate the diff meth
			diff <- calculateDiffMeth(meth)
			# set the threshold
			assign(paste(pair, collapse = ''),getMethylDiff(diff, difference = threshold, qvalue = qvalue))
		}
		# combine all the meth diff results together
		all_diff_combination_result_names <- lapply(pair_combination, function(i) paste(i, collapse = ''))
		all_diff_sites <- sapply(all_diff_combination_result_names, function(i) as(get(i),'GRanges'))
		# get the genes with diff meth sites
		methygenes <- sapply(all_diff_combination_result_names, function(i) subsetByOverlaps(genes, as(get(i), 'GRanges')))
		# get the final methdiffgenes
		diff_genes <- unlist(GRangesList(unlist(methygenes)))
	}
	return(diff_genes)
}

