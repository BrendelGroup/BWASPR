#' cmpSites()
#' This function will 
#'
#' @param sample1hsm A methylRaw object
#' @param sample1scd A methylRaw object
#' @param sample1label A string
#' @param sample2hsm A methylRaw object
#' @param sample3scd A methylRaw object
#' @param sample2label A string
#' @param nbrpms integer
#'
#' @return something
#'
#' @imporFrom methylKit getData
#' @import sqldf
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   nbrpms  <- myfiles$parameters[myfiles$parameters$Variable == "TOTALNBRPMSITES",2]
#'   AmHEhsm <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                           type="CpGhsm", mincov=1,assembly="Amel-4.5")
#'   AmHEscd <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                           type="CpGscd", mincov=1,assembly="Amel-4.5")
#'   cmpSites(getData(AmHEhsm[[1]]),getData(AmHEscd[[1]]),"Am_HE_fr",
#'            getData(AmHEhsm[[2]]),getData(AmHEscd[[2]]),"Am_HE_rn",nbrpms)
#'
#' @export

cmpSites <- function(sample1hsm,sample1scd,sample1label,
		     sample2hsm,sample2scd,sample2label,nbrpms) {
    message("... calculating site statistics ...")
verbose <- 1

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
    cat( sprintf( "================================================================================\n" ) )
    cat( sprintf( "2. Numbers of common and distinct sites comparing the samples\n") )
    cat( sprintf( "================================================================================\n\n" ) )
    
    # ? Which sites are common between samples 1 and 2, and which sites are unique in each sample?
    #
    
    common12 <- sqldf('SELECT sample1hsm.SeqPos as SeqPos, sample1hsm.chr as chr, sample1hsm.start as start, sample1hsm.end as end, sample1hsm.strand as strand, sample1hsm.coverage as s1coverage, sample1hsm.numCs as s1numCs, sample1hsm.numTs as s1numTs, sample1hsm.PrcntM as s1PrcntM, sample2hsm.coverage as s2coverage, sample2hsm.numCs as s2numCs, sample2hsm.numTs as s2numTs, sample2hsm.PrcntM as s2PrcntM from sample1hsm INNER JOIN sample2hsm ON sample1hsm.SeqPos = sample2hsm.SeqPos')
    if ( verbose ) {
      print( head(common12) )
    }
    nbrc12 <- dim(common12)[1]
    
    unique1 <- sqldf('SELECT sample1hsm.SeqPos, sample1hsm.chr, sample1hsm.start, sample1hsm.end, sample1hsm.coverage, sample1hsm.numCs, sample1hsm.numTs from sample1hsm LEFT JOIN sample2hsm ON sample1hsm.SeqPos = sample2hsm.SeqPos WHERE sample2hsm.SeqPos IS NULL')
    if ( verbose ) {
      print( head(unique1) )
    }
    nbru1 <- dim(unique1)[1]
    
    unique2 <- sqldf('SELECT sample2hsm.SeqPos, sample2hsm.chr, sample2hsm.start, sample2hsm.end, sample2hsm.coverage, sample2hsm.numCs, sample2hsm.numTs from sample2hsm LEFT JOIN sample1hsm ON sample2hsm.SeqPos = sample1hsm.SeqPos WHERE sample1hsm.SeqPos IS NULL')
    if ( verbose ) {
      print( head(unique2) )
    }
    nbru2 <- dim(unique2)[1]
    
    
    # ? For a more refined analysis, we find the numbers of common and unique sites only within
    #   the set of sites that were sufficiently covered (i.e., detectable) in both samples:
    #
    
    common34 <- sqldf('SELECT sample1scd.SeqPos, sample1scd.chr, sample1scd.start, sample1scd.end, sample1scd.strand, sample1scd.coverage as s3coverage, sample1scd.numCs as s3numCs, sample1scd.numTs as s3numTs, sample2scd.coverage as s4coverage, sample2scd.numCs as s4numCs, sample2scd.numTs as s4numTs from sample1scd INNER JOIN sample2scd ON sample1scd.SeqPos = sample2scd.SeqPos')
    if ( verbose ) {
      print( head(common34) )
    }
    nbrc34 <- dim(common34)[1]
    
    unique3 <- sqldf('SELECT sample1scd.SeqPos, sample1scd.chr, sample1scd.start, sample1scd.end, sample1scd.coverage, sample1scd.numCs, sample1scd.numTs from sample1scd LEFT JOIN sample2scd ON sample1scd.SeqPos = sample2scd.SeqPos WHERE sample2scd.SeqPos IS NULL')
    if ( verbose ) {
      print( head(unique3) )
    }
    nbru3 <- dim(unique3)[1]
    
    unique4 <- sqldf('SELECT sample2scd.SeqPos, sample2scd.chr, sample2scd.start, sample2scd.end, sample2scd.coverage, sample2scd.numCs, sample2scd.numTs from sample2scd LEFT JOIN sample1scd ON sample2scd.SeqPos = sample1scd.SeqPos WHERE sample1scd.SeqPos IS NULL')
    if ( verbose ) {
      print( head(unique4) )
    }
    nbru4 <- dim(unique4)[1]
    
    common134 <- sqldf('SELECT sample1hsm.SeqPos, sample1hsm.chr, sample1hsm.start, sample1hsm.end, sample1hsm.strand, sample1hsm.coverage, sample1hsm.numCs, sample1hsm.numTs, sample1hsm.PrcntM, common34.s3coverage, common34.s4coverage from sample1hsm INNER JOIN common34 ON sample1hsm.SeqPos = common34.SeqPos')
    if ( verbose ) {
      print( head(common134) )
    }
    nbrc134 <- dim(common134)[1]
    
    unique134 <- sqldf('SELECT sample1hsm.SeqPos, sample1hsm.chr, sample1hsm.start, sample1hsm.end, sample1hsm.coverage, sample1hsm.numCs, sample1hsm.numTs from sample1hsm LEFT JOIN common34 ON sample1hsm.SeqPos = common34.SeqPos WHERE common34.SeqPos IS NULL')
    if ( verbose ) {
      print( head(unique134) )
    }
    nbru134 <- dim(unique134)[1]
    
    common234 <- sqldf('SELECT sample2hsm.SeqPos, sample2hsm.chr, sample2hsm.start, sample2hsm.end, sample2hsm.strand, sample2hsm.coverage, sample2hsm.numCs, sample2hsm.numTs, sample2hsm.PrcntM, common34.s3coverage, common34.s4coverage from sample2hsm INNER JOIN common34 ON sample2hsm.SeqPos = common34.SeqPos')
    if ( verbose ) {
      print( head(common234) )
    }
    nbrc234 <- dim(common234)[1]
    
    unique234 <- sqldf('SELECT sample2hsm.SeqPos, sample2hsm.chr, sample2hsm.start, sample2hsm.end, sample2hsm.coverage, sample2hsm.numCs, sample2hsm.numTs from sample2hsm LEFT JOIN common34 ON sample2hsm.SeqPos = common34.SeqPos WHERE common34.SeqPos IS NULL')
    if ( verbose ) {
      print( head(unique234) )
    }
    nbru234 <- dim(unique234)[1]
    
    unique112 <- sqldf('SELECT common134.* from common134 LEFT JOIN common12 ON common134.SeqPos = common12.SeqPos WHERE common12.SeqPos IS NULL')
    if ( verbose ) {
      print( head(unique112) )
    }
    nbru112 <- dim(unique112)[1]
    
    unique212 <- sqldf('SELECT common234.* from common234 LEFT JOIN common12 ON common234.SeqPos = common12.SeqPos WHERE common12.SeqPos IS NULL')
    if ( verbose ) {
      print( head(unique212) )
    }
    nbru212 <- dim(unique212)[1]
    
    # ! Printing out the results (with expected values):
    cat( sprintf( "  total number of potential sites:\t%9d\n", nbrpms ) )
    cat( sprintf( "\n" ) )
    
    cat( sprintf( "  number of \"%s\"         sites:\t%6d\n", label1hsm, length(sample1hsm$coverage) ) )
    cat( sprintf( "  number of \"%s\"         sites:\t%6d\n", label2hsm, length(sample2hsm$coverage) ) )
    cat( sprintf( "    number of \"%s\"-unique sites:\t%6d\n", label1hsm, nbru1 ) )
    cat( sprintf( "    number of common sites:       \t%6d\n", nbrc12 ) )
    cat( sprintf( "    number of \"%s\"-unique sites:\t%6d\n", label2hsm, nbru2 ) )
    cat( sprintf( "  total number of \"%s+%s\"-sites observed:\t%6d\n", label1hsm, label2hsm, nbru1+nbrc12+nbru2 ) )
    cat( sprintf( "\n" ) )
    
    cat( sprintf( "  number of \"%s\"        sites:\t%9d (%6.2f%% of total)\n", label1scd, length(sample1scd$coverage), 100.*length(sample1scd$coverage)/nbrpms ) )
    cat( sprintf( "  number of \"%s\"        sites:\t%9d (%6.2f%% of total)\n", label2scd, length(sample2scd$coverage), 100.*length(sample2scd$coverage)/nbrpms ) )
    enbrc34 <- round(length(sample2scd$coverage) * (length(sample1scd$coverage) / nbrpms), 0)
    cat( sprintf( "    number of \"%s\"-unique sites:\t%9d\n", label1scd, nbru3 ) )
    cat( sprintf( "    number of sites in common:    \t%9d\t(Expected: %9.0f; O/E: %9.1f)\n", nbrc34, enbrc34, nbrc34/enbrc34) )
    cat( sprintf( "    number of \"%s\"-unique sites:\t%9d\n", label2scd, nbru4 ) )
    cat( sprintf( "  total number of \"%s+%s\"-sites observed:\t%9d\n", label1scd, label2scd, nbru3+nbrc34+nbru4 ) )
    cat( sprintf( "\n" ) )
    
    cat( sprintf( "  number of sites in \"%s\" that are not  detectable in \"%s\":\t%6d\n",label1hsm,label2hsm,nbru134 ) )
    cat( sprintf( "  number of sites in \"%s\" that are also detectable in \"%s\":\t%6d\n",label1hsm,label2hsm,nbrc134 ) )
    cat( sprintf( "    number of sites unique to \"%s\" although detectable in \"%s\":\t\t%6d\n",label1hsm,label2hsm,nbru112 ) )
    enbrc12 <- round(nbrc234 * (nbrc134 / nbrc34), 0)
    cat( sprintf( "    number of \"%s\" / \"%s\" common sites                        :\t\t%6d\t(Expected: %6.0f; O/E: %9.1f)\n",label1hsm,label2hsm,nbrc12,enbrc12,nbrc12/enbrc12 ) )
    cat( sprintf( "\n" ) )
    
    cat( sprintf( "  number of sites in \"%s\" that are not  detectable in \"%s\":\t%6d\n",label2hsm,label1hsm,nbru234 ) )
    cat( sprintf( "  number of sites in \"%s\" that are also detectable in \"%s\":\t%6d\n",label2hsm,label1hsm,nbrc234 ) )
    cat( sprintf( "    number of sites unique to \"%s\" although detectable in \"%s\":\t\t%6d\n",label2hsm,label1hsm,nbru212 ) )
    enbrc21 <- round(nbrc134 * (nbrc234 / nbrc34), 0)
    cat( sprintf( "    number of \"%s\" / \"%s\" common sites                        :\t\t%6d\t(Expected: %6.0f; O/E: %9.1f)\n",label2hsm,label1hsm,nbrc12,enbrc21,nbrc12/enbrc21 ) )
    cat( sprintf( "\n" ) )
    
    # Reporting the (relative) overlap index (= common over total sites ratio, normalized by maximum):
    #
    if (nbru112 <= nbru212) {
      oidx <- (nbrc12 / (nbru112 + nbrc12 + nbru212)) / ((nbru112 + nbrc12)/(nbrc12 + nbru212))
    } else {
      oidx <- (nbrc12 / (nbru112 + nbrc12 + nbru212)) / ((nbru212 + nbrc12)/(nbrc12 + nbru112))
    }
    cat( sprintf( "  Overlap index of \"%s\" with \"%s\":\t\t%5.3f\n",label1hsm,label2hsm,oidx ) )
    
    eps <- nbrc134 * (nbrc234 / nbrc12)
    cat( sprintf( "  Estimated number of \"%s\" = \"%s\" sites (assuming sampling from one population):\t\t%9.0f\n",label1hsm,label2hsm,eps ) )
    cat( sprintf( "  Adjusted population size of \"%s\" = \"%s\" sites (assuming all sites detectable):\t\t%9.0f (%6.2fx of observed)\n",label1hsm,label2hsm,eps*nbrpms/nbrc34,(eps*nbrpms/nbrc34)/(nbru1+nbrc12+nbru2) ) )
    
    cat( sprintf( "\n" ) )
    cat( sprintf( "\n" ) )
    
    message("... done ..")
    return(eps)
}
