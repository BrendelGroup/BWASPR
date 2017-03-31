#' map_methylome()
#' This function maps the scd and hsm sites onto different genomic features
#'   (e.g. genes, exon, promoter, etc.)
#'
#' @param studymk methylRaw object representing an hsm sample
#' @param slabel Label for the sample
#' @param studymc methylRaw object representing the scd control
#' @param clabel Label for the control sample
#' @param genome_ann A list of GRanges objects that contains genome annotation.
#' @param species Label for the species being analyzed
#' @param gnmsize Genome size
#' @param UTRflag Numerical, indicating whether or not the annotation included UTRs (1 or 0)
#' @param outfile If specified, then the result is saved in the specified file name.
#'
#' @return Should return something useful
#'
#' @importFrom BiocGenerics width
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   infiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   asmblv <- infiles$parameters[infiles$parameters$Variable == "ASSEMBLYVERSION",2]
#'   gnmsize <- as.numeric(infiles$parameters[infiles$parameters$Variable == "GENOMESIZE",2])
#'   UTRflag <- as.numeric(infiles$parameters[infiles$parameters$Variable == "UTRFLAGSET",2])
#'   AmHEhsm <- mcalls2mkobj(infiles$datafiles,species="Am",study="HE",
#'                           sample=list("forager","nurse"),replicate=c(0),
#'                           type="CpGhsm",mincov=1,assembly=asmblv
#'                          )
#'   AmHEscd <- mcalls2mkobj(infiles$datafiles,species="Am",study="HE",
#'                           sample=list("forager","nurse"),replicate=c(0),
#'                           type="CpGscd",mincov=1,assembly=asmblv
#'                          )
#'   ginfo <- get_genome_annotation(infiles$parameters)
#'   map_methylome(AmHEhsm[[1]],"forager_hsm",AmHEscd[[1]],"forager_scd",ginfo,
#'                 species="Am",gnmsize,UTRflag=UTRflag,
#'                 outfile="AmHE-methylome-map.txt"
#'                )
#'
#' @export

################################################################################
map_methylome <- function(studymk,slabel,studymc,clabel,
                          genome_ann,species,gnmsize,UTRflag,
                          outfile=""){

    if (outfile != "") {
        sink(outfile)
    }
    gene.gr                <- genome_ann$gene
    exon.gr                <- genome_ann$exon
    pcexon.gr              <- genome_ann$pcexon
    promoter.gr            <- genome_ann$promoter
    CDS.gr                 <- genome_ann$CDS
    fiveprimeUTR.gr        <- genome_ann$fiveprimeUTR
    threeprimeUTR.gr       <- genome_ann$threeprimeUTR
    fiveprimeUTRunique.gr  <- genome_ann$fiveprimeUTRunique
    threeprimeUTRunique.gr <- genome_ann$threeprimeUTRunique
    ncexon.gr              <- genome_ann$ncexon

    # Calculate the fraction of genic and intergenic regions in the genome:
    #
    gene.width                     <- sum(width(gene.gr))
      exon.width                   <- sum(width(exon.gr))
      pcexon.width                 <- sum(width(pcexon.gr))
      ncexon.width                 <- exon.width-pcexon.width
        CDS.width                  <- sum(width(CDS.gr))
        if (UTRflag == 1) {
            fiveprimeUTR.width     <- sum(width(fiveprimeUTR.gr))
            threeprimeUTR.width    <- sum(width(threeprimeUTR.gr))
        }
      intron.width                 <- gene.width - exon.width
    intergenic.width               <- gnmsize - gene.width
      promoter.width               <- sum(width(promoter.gr))
      intergenicsanspromoter.width <- intergenic.width - promoter.width
    
    gene.fraction                     <- gene.width/gnmsize
      exon.fraction                   <- exon.width/gnmsize
      pcexon.fraction                 <- pcexon.width/exon.width
      ncexon.fraction                 <- ncexon.width/exon.width
         CDS.fraction                 <- CDS.width/exon.width
         if (UTRflag == 1) {
            fiveprimeUTR.fraction     <- fiveprimeUTR.width/exon.width
            threeprimeUTR.fraction    <- threeprimeUTR.width/exon.width
         }
      intron.fraction                 <- intron.width/gnmsize
    intergenic.fraction               <- intergenic.width/gnmsize
      promoter.fraction               <- promoter.width/gnmsize
      intergenicsanspromoter.fraction <- intergenicsanspromoter.width/gnmsize
    
    
    ################################################################################
    cat( sprintf( "================================================================================\n" ) )
    cat( sprintf( "Genomic composition in terms of feature regions\n") )
    cat( sprintf( "================================================================================\n\n" ) )
    
    ################################################################################
    
    cat( sprintf("%s genome size:                       %9d bp\n",species,gnmsize) )
    cat( sprintf("  %s genic region size:                 %9d bp (%5.1f%%)\n",species,gene.width,round(100*gene.fraction,1)) )
    cat( sprintf("    %s exon region size:                  %9d bp (%5.1f%%)\n",species,exon.width,round(100*exon.fraction,1)) )
    cat( sprintf("      %s CDS region size:                   %9d bp (%5.1f%%)\n",species,CDS.width,round(100*CDS.fraction,1)) )
    if (UTRflag == 1) {
        cat( sprintf("      %s five-prime UTR region size:        %9d bp (%5.1f%%)\n",species,fiveprimeUTR.width,round(100*fiveprimeUTR.fraction,1)) )
        cat( sprintf("      %s three-prime UTR region size:       %9d bp (%5.1f%%)\n",species,threeprimeUTR.width,round(100*threeprimeUTR.fraction,1)) )
    }
    cat( sprintf("      %s other exon region size:            %9d bp (%5.1f%%)\n",species,ncexon.width,round(100*ncexon.fraction,1)) )
    cat( sprintf("    %s intron region size:                %9d bp (%5.1f%%)\n",species,intron.width,round(100*intron.fraction,1)) )
    cat( sprintf("  %s intergenic region size:            %9d bp (%5.1f%%)\n",species,intergenic.width,round(100*intergenic.fraction,1)) )
    cat( sprintf("    %s promoter region size:              %9d bp (%5.1f%%)\n",species,promoter.width,round(100*promoter.fraction,1)) )
    cat( sprintf("    %s other intergenic region size:      %9d bp (%5.1f%%)\n",species,intergenicsanspromoter.width,round(100*intergenicsanspromoter.fraction,1)) )
    cat( sprintf( "\n\n" ) )
    
    
    ################################################################################
    cat( sprintf( "================================================================================\n" ) )
    cat( sprintf( "Methylation sites in genomic feature regions \n") )
    cat( sprintf( "================================================================================\n\n" ) )
    
    ################################################################################
    
    clabel.gr <- as(studymc,"GRanges")
    clabel.geneoverlap.gr                       <- suppressWarnings((subsetByOverlaps(clabel.gr,gene.gr,ignore.strand=TRUE)))
    clabel.exonoverlap.gr                       <- suppressWarnings((subsetByOverlaps(clabel.gr,exon.gr,ignore.strand=TRUE)))
    clabel.ncexonoverlap.gr                               <- suppressWarnings((subsetByOverlaps(clabel.gr,ncexon.gr,ignore.strand=TRUE)))
    clabel.CDSoverlap.gr                        <- suppressWarnings((subsetByOverlaps(clabel.gr,CDS.gr,ignore.strand=TRUE)))
    if (UTRflag == 1) {
        clabel.fiveprimeUTRoverlap.gr               <- suppressWarnings((subsetByOverlaps(clabel.gr,fiveprimeUTRunique.gr,ignore.strand=TRUE)))
        clabel.threeprimeUTRoverlap.gr              <- suppressWarnings((subsetByOverlaps(clabel.gr,threeprimeUTRunique.gr,ignore.strand=TRUE)))
    }
    clabel.promoteroverlap.gr                   <- suppressWarnings((subsetByOverlaps(clabel.gr,promoter.gr,ignore.strand=TRUE)))
    clabel.SitesInGenicRegions.number           <- length(clabel.geneoverlap.gr) 
    clabel.SitesInExonRegions.number            <- length(clabel.exonoverlap.gr)
    clabel.SitesInNcExonRegions.number            <- length(clabel.ncexonoverlap.gr)
    clabel.SitesInCDSRegions.number             <- length(clabel.CDSoverlap.gr)
    if (UTRflag == 1) {
        clabel.SitesInfiveprimeUTRRegions.number    <- length(clabel.fiveprimeUTRoverlap.gr)
        clabel.SitesInthreeprimeUTRRegions.number   <- length(clabel.threeprimeUTRoverlap.gr)
    }
    clabel.SitesInIntronRegions.number          <- length(clabel.geneoverlap.gr) - length(clabel.exonoverlap.gr)
    clabel.SitesInPromoterRegions.number        <- length(clabel.promoteroverlap.gr)
    clabel.SitesInIntergenicRegions.number      <- length(clabel.gr)- length(clabel.geneoverlap.gr) - length(clabel.promoteroverlap.gr)
    clabel.SitesInFullIntergenicRegions.number  <- clabel.SitesInIntergenicRegions.number + clabel.SitesInPromoterRegions.number
    clabel.SitesInGenicRegions.fraction          <- clabel.SitesInGenicRegions.number/length(clabel.gr)
    clabel.SitesInExonRegions.fraction           <- clabel.SitesInExonRegions.number/length(clabel.gr) 
    clabel.SitesInNcExonRegions.fraction         <- clabel.SitesInNcExonRegions.number/clabel.SitesInExonRegions.number
    clabel.SitesInCDSRegions.fraction            <- clabel.SitesInCDSRegions.number/clabel.SitesInExonRegions.number
    if (UTRflag == 1) {
        clabel.SitesInfiveprimeUTRRegions.fraction   <- clabel.SitesInfiveprimeUTRRegions.number/clabel.SitesInExonRegions.number
        clabel.SitesInthreeprimeUTRRegions.fraction  <- clabel.SitesInthreeprimeUTRRegions.number/clabel.SitesInExonRegions.number
    }
    clabel.SitesInIntronRegions.fraction         <- clabel.SitesInIntronRegions.number/length(clabel.gr)
    clabel.SitesInPromoterRegions.fraction       <- clabel.SitesInPromoterRegions.number/length(clabel.gr)
    clabel.SitesInIntergenicRegions.fraction     <- clabel.SitesInIntergenicRegions.number/length(clabel.gr)
    clabel.SitesInFullIntergenicRegions.fraction <- clabel.SitesInFullIntergenicRegions.number/length(clabel.gr)
    
    clabel.density                              <- ((length(clabel.gr) * 10000) / gnmsize)
    clabel.SitesInGenicRegions.density          <- ((clabel.SitesInGenicRegions.number * 10000) / gene.width)
    clabel.SitesInExonRegions.density           <- ((clabel.SitesInExonRegions.number  * 10000) / exon.width)
    clabel.SitesInNcExonRegions.density         <- ((clabel.SitesInNcExonRegions.number  * 10000) / ncexon.width)
    clabel.SitesInCDSRegions.density            <- ((clabel.SitesInCDSRegions.number   *10000)/CDS.width)
    if (UTRflag ==1) {
        clabel.SitesInfiveprimeUTRRegions.density   <- ((clabel.SitesInfiveprimeUTRRegions.number *10000)/ fiveprimeUTR.width)
        clabel.SitesInthreeprimeUTRRegions.density  <- ((clabel.SitesInthreeprimeUTRRegions.number *10000)/ threeprimeUTR.width)
    }
    clabel.SitesInIntronRegions.density         <- ((clabel.SitesInIntronRegions.number * 10000) /intron.width)
    clabel.SitesInFullIntergenicRegions.density <- ((clabel.SitesInFullIntergenicRegions.number * 10000) / intergenic.width)
    clabel.SitesInPromoterRegions.density       <- ((clabel.SitesInPromoterRegions.number * 10000) / promoter.width)
    clabel.SitesInIntergenicRegions.density     <- ((clabel.SitesInIntergenicRegions.number * 10000) / intergenicsanspromoter.width)
    
    cat( sprintf( "Number and density (per 10kb) of sites in %s\n\n",clabel) )
    cat( sprintf("Number of sites identified in %s %s:                            %9d     (%8.2f)\n",species,clabel,length(clabel.gr),clabel.density) )
    cat( sprintf("Number of sites identified in %s %s genic regions:              %9d     (%8.2f)                      %6.2f%%     (%5.2f O/E)\n",species,clabel,clabel.SitesInGenicRegions.number,clabel.SitesInGenicRegions.density,100*clabel.SitesInGenicRegions.fraction,clabel.SitesInGenicRegions.fraction/gene.fraction) )   
    cat( sprintf("Number of sites identified in %s %s   exon regions:               %9d   (%8.2f)                        %6.2f%%   (%5.2f O/E)\n",species,clabel,clabel.SitesInExonRegions.number,clabel.SitesInExonRegions.density,100*clabel.SitesInExonRegions.fraction,clabel.SitesInExonRegions.fraction/exon.fraction) )
    cat( sprintf("Number of sites identified in %s %s     CDS regions:                %9d (%8.2f)                          %6.2f%% (%5.2f O/E)\n",species,clabel,clabel.SitesInCDSRegions.number,clabel.SitesInCDSRegions.density,100*clabel.SitesInCDSRegions.fraction,clabel.SitesInCDSRegions.fraction/CDS.fraction) )
    if (UTRflag == 1) {
        cat( sprintf("Number of sites identified in %s %s     five-prime UTR regions:     %9d (%8.2f)                          %6.2f%% (%5.2f O/E)\n",species,clabel,clabel.SitesInfiveprimeUTRRegions.number,clabel.SitesInfiveprimeUTRRegions.density,100*clabel.SitesInfiveprimeUTRRegions.fraction,clabel.SitesInfiveprimeUTRRegions.fraction/fiveprimeUTR.fraction) )
        cat( sprintf("Number of sites identified in %s %s     three-prime UTR regions:    %9d (%8.2f)                          %6.2f%% (%5.2f O/E)\n",species,clabel,clabel.SitesInthreeprimeUTRRegions.number,clabel.SitesInthreeprimeUTRRegions.density,100*clabel.SitesInthreeprimeUTRRegions.fraction,clabel.SitesInthreeprimeUTRRegions.fraction/threeprimeUTR.fraction) )
    }
    cat( sprintf("Number of sites identified in %s %s     other exon regions:         %9d (%8.2f)                          %6.2f%% (%5.2f O/E)\n",species,clabel,clabel.SitesInNcExonRegions.number,clabel.SitesInNcExonRegions.density,100*clabel.SitesInNcExonRegions.fraction,clabel.SitesInNcExonRegions.fraction/ncexon.fraction) )
    cat( sprintf("Number of sites identified in %s %s   intron regions:             %9d   (%8.2f)                        %6.2f%%   (%5.2f O/E)\n",species,clabel,clabel.SitesInIntronRegions.number,clabel.SitesInIntronRegions.density,100*clabel.SitesInIntronRegions.fraction,clabel.SitesInIntronRegions.fraction/intron.fraction) )
    cat( sprintf("Number of sites identified in %s %s intergenic regions:         %9d     (%8.2f)                      %6.2f%%     (%5.2f O/E)\n",species,clabel,clabel.SitesInFullIntergenicRegions.number,clabel.SitesInFullIntergenicRegions.density,100*clabel.SitesInFullIntergenicRegions.fraction,clabel.SitesInFullIntergenicRegions.fraction/intergenic.fraction) )
    cat( sprintf("Number of sites identified in %s %s   promoter regions:           %9d   (%8.2f)                        %6.2f%%   (%5.2f O/E)\n",species,clabel,clabel.SitesInPromoterRegions.number,clabel.SitesInPromoterRegions.density,100*clabel.SitesInPromoterRegions.fraction,clabel.SitesInPromoterRegions.fraction/promoter.fraction) )
    cat( sprintf("Number of sites identified in %s %s   other intergenic regions:   %9d   (%8.2f)                        %6.2f%%   (%5.2f O/E)\n",species,clabel,clabel.SitesInIntergenicRegions.number,clabel.SitesInIntergenicRegions.density,100*clabel.SitesInIntergenicRegions.fraction,clabel.SitesInIntergenicRegions.fraction/intergenicsanspromoter.fraction) )   
    cat( sprintf( "\n\n" ) )
    
    
    slabel.gr <- as(studymk,"GRanges")
    slabel.geneoverlap.gr                        <- suppressWarnings((subsetByOverlaps(slabel.gr,gene.gr,ignore.strand=TRUE)))
    slabel.exonoverlap.gr                        <- suppressWarnings((subsetByOverlaps(slabel.gr,exon.gr,ignore.strand=TRUE)))
    slabel.ncexonoverlap.gr                      <- suppressWarnings((subsetByOverlaps(slabel.gr,ncexon.gr,ignore.strand=TRUE)))
    slabel.CDSoverlap.gr                         <- suppressWarnings((subsetByOverlaps(slabel.gr,CDS.gr,ignore.strand=TRUE)))
    if (UTRflag == 1) {
        slabel.fiveprimeUTRoverlap.gr                <- suppressWarnings((subsetByOverlaps(slabel.gr,fiveprimeUTRunique.gr,ignore.strand=TRUE)))
        slabel.threeprimeUTRoverlap.gr               <- suppressWarnings((subsetByOverlaps(slabel.gr,threeprimeUTRunique.gr,ignore.strand=TRUE)))
    }
    slabel.promoteroverlap.gr                    <- suppressWarnings((subsetByOverlaps(slabel.gr,promoter.gr,ignore.strand=TRUE)))
    slabel.SitesInGenicRegions.number            <- length(slabel.geneoverlap.gr) 
    slabel.SitesInExonRegions.number             <- length(slabel.exonoverlap.gr)
    slabel.SitesInNcExonRegions.number           <- length(slabel.ncexonoverlap.gr)
    slabel.SitesInCDSRegions.number              <- length(slabel.CDSoverlap.gr)
    if (UTRflag == 1){
        slabel.SitesInfiveprimeUTRRegions.number     <- length(slabel.fiveprimeUTRoverlap.gr)
        slabel.SitesInthreeprimeUTRRegions.number    <- length(slabel.threeprimeUTRoverlap.gr)
    }
    slabel.SitesInIntronRegions.number           <- length(slabel.geneoverlap.gr) - length(slabel.exonoverlap.gr)
    slabel.SitesInPromoterRegions.number         <- length(slabel.promoteroverlap.gr)
    slabel.SitesInIntergenicRegions.number       <- length(slabel.gr)- length(slabel.geneoverlap.gr) - length(slabel.promoteroverlap.gr)
    slabel.SitesInFullIntergenicRegions.number   <- slabel.SitesInIntergenicRegions.number + slabel.SitesInPromoterRegions.number
    slabel.SitesInGenicRegions.fraction          <- slabel.SitesInGenicRegions.number/length(slabel.gr)
    slabel.SitesInExonRegions.fraction           <- slabel.SitesInExonRegions.number/length(slabel.gr) 
    slabel.SitesInNcExonRegions.fraction         <- slabel.SitesInNcExonRegions.number/slabel.SitesInExonRegions.number
    slabel.SitesInCDSRegions.fraction            <- slabel.SitesInCDSRegions.number/slabel.SitesInExonRegions.number
    if (UTRflag == 1) {
        slabel.SitesInfiveprimeUTRRegions.fraction   <- slabel.SitesInfiveprimeUTRRegions.number/slabel.SitesInExonRegions.number
        slabel.SitesInthreeprimeUTRRegions.fraction  <- slabel.SitesInthreeprimeUTRRegions.number/slabel.SitesInExonRegions.number
    }
    slabel.SitesInIntronRegions.fraction         <- slabel.SitesInIntronRegions.number/length(slabel.gr)
    slabel.SitesInPromoterRegions.fraction       <- slabel.SitesInPromoterRegions.number/length(slabel.gr)
    slabel.SitesInIntergenicRegions.fraction     <- slabel.SitesInIntergenicRegions.number/length(slabel.gr)
    slabel.SitesInFullIntergenicRegions.fraction <- slabel.SitesInFullIntergenicRegions.number/length(slabel.gr)
    
    slabel.density                              <- ((length(slabel.gr) * 10000) / gnmsize)
    slabel.SitesInGenicRegions.density          <- ((slabel.SitesInGenicRegions.number * 10000) / gene.width)
    slabel.SitesInExonRegions.density           <- ((slabel.SitesInExonRegions.number  * 10000) / exon.width)
    slabel.SitesInNcExonRegions.density         <- ((slabel.SitesInNcExonRegions.number  * 10000) / ncexon.width)
    slabel.SitesInCDSRegions.density            <- ((slabel.SitesInCDSRegions.number   *10000)/CDS.width)
    if (UTRflag ==1) {
        slabel.SitesInfiveprimeUTRRegions.density   <- ((slabel.SitesInfiveprimeUTRRegions.number *10000)/ fiveprimeUTR.width)
        slabel.SitesInthreeprimeUTRRegions.density  <- ((slabel.SitesInthreeprimeUTRRegions.number *10000)/ threeprimeUTR.width)
    }
    slabel.SitesInIntronRegions.density         <- ((slabel.SitesInIntronRegions.number * 10000) /intron.width)
    slabel.SitesInFullIntergenicRegions.density <- ((slabel.SitesInFullIntergenicRegions.number * 10000) / intergenic.width)
    slabel.SitesInPromoterRegions.density       <- ((slabel.SitesInPromoterRegions.number * 10000) / promoter.width)
    slabel.SitesInIntergenicRegions.density     <- ((slabel.SitesInIntergenicRegions.number * 10000) / intergenicsanspromoter.width)
    
    cat( sprintf( "Number and density (per 10kb) of sites in %s\n\n",slabel) )
    cat( sprintf("Number of sites identified in %s %s:                            %9d     (%8.2f) (%6.2f%% of control)\n",species,slabel,length(slabel.gr),slabel.density,100*length(slabel.gr)/length(clabel.gr)) )
    cat( sprintf("Number of sites identified in %s %s genic regions:              %9d     (%8.2f) (%6.2f%% of control) %6.2f%%     (%5.2f O/E)\n",species,slabel,slabel.SitesInGenicRegions.number,slabel.SitesInGenicRegions.density,100*slabel.SitesInGenicRegions.number/clabel.SitesInGenicRegions.number,100*slabel.SitesInGenicRegions.fraction,slabel.SitesInGenicRegions.fraction/gene.fraction) )   
    cat( sprintf("Number of sites identified in %s %s   exon regions:               %9d   (%8.2f) (%6.2f%% of control)   %6.2f%%   (%5.2f O/E)\n",species,slabel,slabel.SitesInExonRegions.number,slabel.SitesInExonRegions.density,100*slabel.SitesInExonRegions.number/clabel.SitesInExonRegions.number,100*slabel.SitesInExonRegions.fraction,slabel.SitesInExonRegions.fraction/exon.fraction) )
    cat( sprintf("Number of sites identified in %s %s     CDS regions:                %9d (%8.2f) (%6.2f%% of control)     %6.2f%% (%5.2f O/E)\n",species,slabel,slabel.SitesInCDSRegions.number,slabel.SitesInCDSRegions.density,100*slabel.SitesInCDSRegions.number/clabel.SitesInCDSRegions.number,100*slabel.SitesInCDSRegions.fraction,slabel.SitesInCDSRegions.fraction/CDS.fraction) )
    if (UTRflag == 1) {
    cat( sprintf("Number of sites identified in %s %s     five-prime UTR regions:     %9d (%8.2f) (%6.2f%% of control)     %6.2f%% (%5.2f O/E)\n",species,slabel,slabel.SitesInfiveprimeUTRRegions.number,slabel.SitesInfiveprimeUTRRegions.density,100*slabel.SitesInfiveprimeUTRRegions.number/clabel.SitesInfiveprimeUTRRegions.number,100*slabel.SitesInfiveprimeUTRRegions.fraction,slabel.SitesInfiveprimeUTRRegions.fraction/fiveprimeUTR.fraction) )
    cat( sprintf("Number of sites identified in %s %s     three-prime UTR regions:    %9d (%8.2f) (%6.2f%% of control)     %6.2f%% (%5.2f O/E)\n",species,slabel,slabel.SitesInthreeprimeUTRRegions.number,slabel.SitesInthreeprimeUTRRegions.density,100*slabel.SitesInthreeprimeUTRRegions.number/clabel.SitesInthreeprimeUTRRegions.number,100*slabel.SitesInthreeprimeUTRRegions.fraction,slabel.SitesInthreeprimeUTRRegions.fraction/threeprimeUTR.fraction) )
    }
    cat( sprintf("Number of sites identified in %s %s     other exon regions:         %9d (%8.2f) (%6.2f%% of control)     %6.2f%% (%5.2f O/E)\n",species,slabel,slabel.SitesInNcExonRegions.number,slabel.SitesInNcExonRegions.density,100*slabel.SitesInNcExonRegions.number/clabel.SitesInNcExonRegions.number,100*slabel.SitesInNcExonRegions.fraction,slabel.SitesInNcExonRegions.fraction/ncexon.fraction) )
    cat( sprintf("Number of sites identified in %s %s   intron regions:             %9d   (%8.2f) (%6.2f%% of control)   %6.2f%%   (%5.2f O/E)\n",species,slabel,slabel.SitesInIntronRegions.number,slabel.SitesInIntronRegions.density,100*slabel.SitesInIntronRegions.number/clabel.SitesInIntronRegions.number,100*slabel.SitesInIntronRegions.fraction,slabel.SitesInIntronRegions.fraction/intron.fraction) )
    cat( sprintf("Number of sites identified in %s %s intergenic regions:         %9d     (%8.2f) (%6.2f%% of control) %6.2f%%     (%5.2f O/E)\n",species,slabel,slabel.SitesInFullIntergenicRegions.number,slabel.SitesInFullIntergenicRegions.density,100*slabel.SitesInFullIntergenicRegions.number/clabel.SitesInFullIntergenicRegions.number,100*slabel.SitesInFullIntergenicRegions.fraction,slabel.SitesInFullIntergenicRegions.fraction/intergenic.fraction) )
    cat( sprintf("Number of sites identified in %s %s   promoter regions:           %9d   (%8.2f) (%6.2f%% of control)   %6.2f%%   (%5.2f O/E)\n",species,slabel,slabel.SitesInPromoterRegions.number,slabel.SitesInPromoterRegions.density,100*slabel.SitesInPromoterRegions.number/clabel.SitesInPromoterRegions.number,100*slabel.SitesInPromoterRegions.fraction,slabel.SitesInPromoterRegions.fraction/promoter.fraction) )
    cat( sprintf("Number of sites identified in %s %s   other intergenic regions:   %9d   (%8.2f) (%6.2f%% of control)   %6.2f%%   (%5.2f O/E)\n",species,slabel,slabel.SitesInIntergenicRegions.number,slabel.SitesInIntergenicRegions.density,100*slabel.SitesInIntergenicRegions.number/clabel.SitesInIntergenicRegions.number,100*slabel.SitesInIntergenicRegions.fraction,slabel.SitesInIntergenicRegions.fraction/intergenicsanspromoter.fraction) )   
    
    cat( sprintf( "\n\n" ) )
    if (outfile != "") {
        sink()
    }
    return(1)
}
