#' cmpSites()
#' This function compares (potential) methylation sites between two samples.
#'
#' @param sample1hsm A methylRaw object for hsm sites of sample 1.
#' @param sample1scd A methylRaw object for scd sites of sample 1.
#' @param sample1label A string serving as a label for sample 1.
#' @param sample2hsm A methylRaw object for hsm sites of sample 2.
#' @param sample2scd A methylRaw object for scd sites of stample 2.
#' @param sample2label A string serving as a label for sample 2.
#' @param nbrpms Integer representing the number of potential methylation
#'   sites in the genome; typically derived from the input *.par file.
#'
#' @return A list of data frames containing data on unique and common
#'   sites comparing the two samples.
#'
#' @importFrom methylKit getData
#' @import sqldf
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   nbrpms  <- as.numeric(myfiles$parameters[
#'     myfiles$parameters$Variable == "TOTALNBRPMSITES",2])
#'   AmHEhsm <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                           type="CpGhsm", mincov=1,assembly="Amel-4.5")
#'   AmHEscd <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                           type="CpGscd", mincov=1,assembly="Amel-4.5")
#'   s1hsm <- methylKit::getData(AmHEhsm[[1]])
#'   s1scd <- methylKit::getData(AmHEscd[[1]])
#'   s2hsm <- methylKit::getData(AmHEhsm[[2]])
#'   s2scd <- methylKit::getData(AmHEscd[[2]])
#'   mydflist <- cmpSites(s1hsm,s1scd,"Am_HE_fr",s2hsm,s2scd,"Am_HE_rn",nbrpms)
#'
#' @export

cmpSites <- function(sample1hsm,sample1scd,sample1label,
                     sample2hsm,sample2scd,sample2label,nbrpms) {
    message("... calculating site statistics ...")

    # ... adding a unique key to the data frame rows:
    sample1hsm$SeqPos <- paste(sample1hsm$chr,sample1hsm$start,sep=".")
    sample1hsm$PrcntM <- 100.*sample1hsm$numCs/sample1hsm$coverage
    label1hsm <- paste(sample1label,"hsm",sep="_")
    
    sample1scd$SeqPos <- paste(sample1scd$chr,sample1scd$start,sep=".")
    sample1scd$PrcntM <- 100.*sample1scd$numCs/sample1scd$coverage
    label1scd <- paste(sample1label,"scd",sep="_")
    
    sample2hsm$SeqPos <- paste(sample2hsm$chr,sample2hsm$start,sep=".")
    sample2hsm$PrcntM <- 100.*sample2hsm$numCs/sample2hsm$coverage
    label2hsm <- paste(sample2label,"hsm",sep="_")
    
    sample2scd$SeqPos <- paste(sample2scd$chr,sample2scd$start,sep=".")
    sample2scd$PrcntM <- 100.*sample2scd$numCs/sample2scd$coverage
    label2scd <- paste(sample2label,"scd",sep="_")
    
    # ... accounting for common and unique sites comparing the two samples:
    #
    cat( sprintf( "Numbers of common and distinct sites comparing %s versus %s\n",sample1label,sample2label) )
    cat( sprintf( "================================================================================\n\n" ) )
    
    # ? Which sites are common between samples 1 and 2, and which sites are unique in each sample?
    #
    commonHSM <- sqldf('SELECT sample1hsm.SeqPos as SeqPos, sample1hsm.chr as chr, sample1hsm.start as start, sample1hsm.end as end, sample1hsm.strand as strand, sample1hsm.coverage as covHSM1, sample1hsm.numCs as hsm1numCs, sample1hsm.numTs as hsm1numTs, sample1hsm.PrcntM as s1PrcntM, sample2hsm.coverage as covHSM2, sample2hsm.numCs as hsm2numCs, sample2hsm.numTs as hsm2numTs, sample2hsm.PrcntM as s2PrcntM from sample1hsm INNER JOIN sample2hsm ON sample1hsm.SeqPos = sample2hsm.SeqPos')
    NBRcommonHSM <- dim(commonHSM)[1]
    
    unique1HSM<- sqldf('SELECT sample1hsm.SeqPos, sample1hsm.chr, sample1hsm.start, sample1hsm.end, sample1hsm.coverage, sample1hsm.numCs, sample1hsm.numTs from sample1hsm LEFT JOIN sample2hsm ON sample1hsm.SeqPos = sample2hsm.SeqPos WHERE sample2hsm.SeqPos IS NULL')
    NBRunique1HSM <- dim(unique1HSM)[1]
    
    unique2HSM <- sqldf('SELECT sample2hsm.SeqPos, sample2hsm.chr, sample2hsm.start, sample2hsm.end, sample2hsm.coverage, sample2hsm.numCs, sample2hsm.numTs from sample2hsm LEFT JOIN sample1hsm ON sample2hsm.SeqPos = sample1hsm.SeqPos WHERE sample1hsm.SeqPos IS NULL')
    NBRunique2HSM <- dim(unique2HSM)[1]
    
    
    # ? For a more refined analysis, we find the numbers of common and unique sites only within
    #   the set of sites that were sufficiently covered (i.e., detectable) in both samples:
    #
    commonSCD <- sqldf('SELECT sample1scd.SeqPos, sample1scd.chr, sample1scd.start, sample1scd.end, sample1scd.strand, sample1scd.coverage as covSCD1, sample1scd.numCs as scd1numCs, sample1scd.numTs as scd1numTs, sample2scd.coverage as covSCD2, sample2scd.numCs as scd2numCs, sample2scd.numTs as scd2numTs from sample1scd INNER JOIN sample2scd ON sample1scd.SeqPos = sample2scd.SeqPos')
    NBRcommonSCD <- dim(commonSCD)[1]
    
    unique1SCD <- sqldf('SELECT sample1scd.SeqPos, sample1scd.chr, sample1scd.start, sample1scd.end, sample1scd.coverage, sample1scd.numCs, sample1scd.numTs from sample1scd LEFT JOIN sample2scd ON sample1scd.SeqPos = sample2scd.SeqPos WHERE sample2scd.SeqPos IS NULL')
    NBRunique1SCD <- dim(unique1SCD)[1]
    
    unique2SCD <- sqldf('SELECT sample2scd.SeqPos, sample2scd.chr, sample2scd.start, sample2scd.end, sample2scd.coverage, sample2scd.numCs, sample2scd.numTs from sample2scd LEFT JOIN sample1scd ON sample2scd.SeqPos = sample1scd.SeqPos WHERE sample1scd.SeqPos IS NULL')
    NBRunique2SCD <- dim(unique2SCD)[1]
    
    hsm1cSCD <- sqldf('SELECT sample1hsm.SeqPos, sample1hsm.chr, sample1hsm.start, sample1hsm.end, sample1hsm.strand, sample1hsm.coverage, sample1hsm.numCs, sample1hsm.numTs, sample1hsm.PrcntM, commonSCD.covSCD1, commonSCD.covSCD2 from sample1hsm INNER JOIN commonSCD ON sample1hsm.SeqPos = commonSCD.SeqPos')
    NBRhsm1cSCD <- dim(hsm1cSCD)[1]
    
    unique1HSMn2SCD <- sqldf('SELECT sample1hsm.SeqPos, sample1hsm.chr, sample1hsm.start, sample1hsm.end, sample1hsm.coverage, sample1hsm.numCs, sample1hsm.numTs from sample1hsm LEFT JOIN commonSCD ON sample1hsm.SeqPos = commonSCD.SeqPos WHERE commonSCD.SeqPos IS NULL')
    NBRunique1HSMn2SCD <- dim(unique1HSMn2SCD)[1]
    
    hsm2cSCD <- sqldf('SELECT sample2hsm.SeqPos, sample2hsm.chr, sample2hsm.start, sample2hsm.end, sample2hsm.strand, sample2hsm.coverage, sample2hsm.numCs, sample2hsm.numTs, sample2hsm.PrcntM, commonSCD.covSCD1, commonSCD.covSCD2 from sample2hsm INNER JOIN commonSCD ON sample2hsm.SeqPos = commonSCD.SeqPos')
    NBRhsm2cSCD <- dim(hsm2cSCD)[1]
    
    unique2HSMn1SCD <- sqldf('SELECT sample2hsm.SeqPos, sample2hsm.chr, sample2hsm.start, sample2hsm.end, sample2hsm.coverage, sample2hsm.numCs, sample2hsm.numTs from sample2hsm LEFT JOIN commonSCD ON sample2hsm.SeqPos = commonSCD.SeqPos WHERE commonSCD.SeqPos IS NULL')
    NBRunique2HSMn1SCD <- dim(unique2HSMn1SCD)[1]
    
    unique1HSM2SCD <- sqldf('SELECT hsm1cSCD.* from hsm1cSCD LEFT JOIN commonHSM ON hsm1cSCD.SeqPos = commonHSM.SeqPos WHERE commonHSM.SeqPos IS NULL')
    NBRunique1HSM2SCD <- dim(unique1HSM2SCD)[1]
    
    unique2HSM1SCD <- sqldf('SELECT hsm2cSCD.* from hsm2cSCD LEFT JOIN commonHSM ON hsm2cSCD.SeqPos = commonHSM.SeqPos WHERE commonHSM.SeqPos IS NULL')
    NBRunique2HSM1SCD <- dim(unique2HSM1SCD)[1]
    
    # ! Printing out the results (with expected values):
    cat( sprintf( "  total number of potential sites:\t%9d\n", nbrpms ) )
    cat( sprintf( "\n" ) )
    
    cat( sprintf( "  number of \"%s\"         sites:\t%6d\n", label1hsm, length(sample1hsm$coverage) ) )
    cat( sprintf( "  number of \"%s\"         sites:\t%6d\n", label2hsm, length(sample2hsm$coverage) ) )
    cat( sprintf( "    number of \"%s\"-unique sites:\t%6d\n", label1hsm, NBRunique1HSM ) )
    cat( sprintf( "    number of common sites:       \t%6d\n", NBRcommonHSM ) )
    cat( sprintf( "    number of \"%s\"-unique sites:\t%6d\n", label2hsm, NBRunique2HSM ) )
    cat( sprintf( "  total number of \"%s+%s\"-sites observed:\t%6d\n", label1hsm, label2hsm, NBRunique1HSM+NBRcommonHSM+NBRunique2HSM ) )
    cat( sprintf( "\n" ) )
    
    cat( sprintf( "  number of \"%s\"        sites:\t%9d (%6.2f%% of total)\n", label1scd, length(sample1scd$coverage), 100.*length(sample1scd$coverage)/nbrpms ) )
    cat( sprintf( "  number of \"%s\"        sites:\t%9d (%6.2f%% of total)\n", label2scd, length(sample2scd$coverage), 100.*length(sample2scd$coverage)/nbrpms ) )
    eNBRcommonSCD <- round(length(sample2scd$coverage) * (length(sample1scd$coverage) / nbrpms), 0)
    cat( sprintf( "    number of \"%s\"-unique sites:\t%9d\n", label1scd, NBRunique1SCD ) )
    cat( sprintf( "    number of sites in common:    \t%9d\t(Expected: %9.0f; O/E: %9.1f)\n", NBRcommonSCD, eNBRcommonSCD, NBRcommonSCD/eNBRcommonSCD) )
    cat( sprintf( "    number of \"%s\"-unique sites:\t%9d\n", label2scd, NBRunique2SCD ) )
    cat( sprintf( "  total number of \"%s+%s\"-sites observed:\t%9d\n", label1scd, label2scd, NBRunique1SCD+NBRcommonSCD+NBRunique2SCD ) )
    cat( sprintf( "\n" ) )
    
    cat( sprintf( "  number of sites in \"%s\" that are not  detectable in \"%s\":\t%6d\n",label1hsm,label2hsm,NBRunique1HSMn2SCD ) )
    cat( sprintf( "  number of sites in \"%s\" that are also detectable in \"%s\":\t%6d\n",label1hsm,label2hsm,NBRhsm1cSCD ) )
    cat( sprintf( "    number of sites unique to \"%s\" although detectable in \"%s\":\t\t%6d\n",label1hsm,label2hsm,NBRunique1HSM2SCD ) )
    eNBRcommonHSM <- round(NBRhsm2cSCD * (NBRhsm1cSCD / NBRcommonSCD), 0)
    cat( sprintf( "    number of \"%s\" / \"%s\" common sites                        :\t\t%6d\t(Expected: %6.0f; O/E: %9.1f)\n",label1hsm,label2hsm,NBRcommonHSM,eNBRcommonHSM,NBRcommonHSM/eNBRcommonHSM ) )
    cat( sprintf( "\n" ) )
    
    cat( sprintf( "  number of sites in \"%s\" that are not  detectable in \"%s\":\t%6d\n",label2hsm,label1hsm,NBRunique2HSMn1SCD ) )
    cat( sprintf( "  number of sites in \"%s\" that are also detectable in \"%s\":\t%6d\n",label2hsm,label1hsm,NBRhsm2cSCD ) )
    cat( sprintf( "    number of sites unique to \"%s\" although detectable in \"%s\":\t\t%6d\n",label2hsm,label1hsm,NBRunique2HSM1SCD ) )
    cat( sprintf( "    number of \"%s\" / \"%s\" common sites                        :\t\t%6d\t(Expected: %6.0f; O/E: %9.1f)\n",label2hsm,label1hsm,NBRcommonHSM,eNBRcommonHSM,NBRcommonHSM/eNBRcommonHSM ) )
    cat( sprintf( "\n" ) )
    
    # Reporting the (relative) overlap index (= common over total sites ratio, normalized by maximum):
    #
    if (NBRunique1HSM2SCD <= NBRunique2HSM1SCD) {
      oidx <- (NBRcommonHSM / (NBRunique1HSM2SCD + NBRcommonHSM + NBRunique2HSM1SCD)) / ((NBRunique1HSM2SCD + NBRcommonHSM)/(NBRcommonHSM + NBRunique2HSM1SCD))
    } else {
      oidx <- (NBRcommonHSM / (NBRunique1HSM2SCD + NBRcommonHSM + NBRunique2HSM1SCD)) / ((NBRunique2HSM1SCD + NBRcommonHSM)/(NBRcommonHSM + NBRunique1HSM2SCD))
    }
    cat( sprintf( "  Overlap index of \"%s\" with \"%s\":\t\t%5.3f\n",label1hsm,label2hsm,oidx ) )
    
    eps <- NBRhsm1cSCD * (NBRhsm2cSCD / NBRcommonHSM)
    cat( sprintf( "  Estimated number of \"%s\" = \"%s\" sites (assuming sampling from one population):\t\t%9.0f\n",label1hsm,label2hsm,eps ) )
    cat( sprintf( "  Adjusted population size of \"%s\" = \"%s\" sites (assuming all sites detectable):\t\t%9.0f (%6.2fx of observed)\n",label1hsm,label2hsm,eps*nbrpms/NBRcommonSCD,(eps*nbrpms/NBRcommonSCD)/(NBRunique1HSM+NBRcommonHSM+NBRunique2HSM) ) )
    
    cat( sprintf( "\n" ) )
    cat( sprintf( "\n" ) )
    
    message("... done ..")
    return(list(commonHSM, unique1HSM, unique2HSM,
                commonSCD, unique1SCD, unique2SCD,
                hsm1cSCD, unique1HSMn2SCD,
                hsm2cSCD, unique2HSMn1SCD,
                unique1HSM2SCD, unique2HSM1SCD)
          )
}
