#' map_mrpr()
#' This function maps the scd and hsm sites onto different genomic features
#'   (e.g. genes, exon, promoter, etc.)
#'
#' @param hsmrL List of data frames with hsm rich and poor regions, as put out by
#    det_mrpr()
#' @param species Label for the species being analyzed
#' @param slabel Label for the sample
#' @param genome_ann A list of GRanges objects that contains genome annotation.
#' @param gnmsize Genome size
#' @param UTRflag Numerical, indicating whether or not the annotation included UTRs (1 or 0)
#' @param outfile If specified, then the result is saved in the specified file name.
#'
#' @return Data frames describing genes overlapping with methylation-rich and poor regions
#'
#' @import GenomicRanges
#' @import IRanges
#' @importFrom S4Vectors subjectHits queryHits
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   infiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   asmblv <- infiles$parameters[infiles$parameters$Variable == "ASSEMBLYVERSION",2]
#'   gnmsize <- as.numeric(infiles$parameters[infiles$parameters$Variable == "GENOMESIZE",2])
#'   UTRflag <- as.numeric(infiles$parameters[infiles$parameters$Variable == "UTRFLAGSET",2])
#'   AmHE <- mcalls2mkobj(infiles$datafiles,species="Am",study="HE",
#'                        sample=list("forager","nurse"),replicate=c(0),
#'                        type="CpGhsm",mincov=1,assembly=asmblv
#'                       )
#'   hsmrL <- det_mrpr(AmHE[[1]],"Am_HE_fr",ddset=c(1,5),outfile="dst-Am_HE_fr.txt",nr2d=10L,doplots=TRUE)
#'   ginfo <- get_genome_annotation(infiles$parameters)
#'   map_mrpr(hsmrL,species="Am",slabel="forager_hsm",ginfo,gnmsize,UTRflag=UTRflag,
#'            outfile="AmHE-region.map.txt"
#'           )
#'
#' @export

################################################################################
map_mrpr <- function(hsmrL,species,slabel,genome_ann,gnmsize,UTRflag,outfile) {

    if (outfile != "") {
        sink(outfile)
    }
    gene.gr                <- genome_ann$gene
    exon.gr                <- genome_ann$exon
    pcexon.gr              <- genome_ann$pcexon
    promoter.gr            <- genome_ann$promoter
    CDS.gr                 <- genome_ann$CDS
    fiveprimeUTR.gr        <- genome_ann$fpUTR
    threeprimeUTR.gr       <- genome_ann$tpUTR
    fiveprimeUTRunique.gr  <- genome_ann$fpUTRnotCDS
    threeprimeUTRunique.gr <- genome_ann$tpUTRunique
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
    
    
    #hsm-rich regions:
    #
    tmpdf <- hsmrL$hsmrR
    hsmrR.gr <- with(tmpdf,{GRanges(tmpdf$SeqID,
                                    IRanges(tmpdf$From,tmpdf$To),
                                    NbrSites=(tmpdf$NbrSites),                                               
                                    Sdnsty=(tmpdf$Sdnsty))
                           }
                           )
    hsmrR.geneoverlap.gr           <- suppressWarnings((intersect(hsmrR.gr,gene.gr,ignore.strand=TRUE)))
    hsmrR.exonoverlap.gr           <- suppressWarnings((intersect(hsmrR.gr,exon.gr,ignore.strand=TRUE)))
    hsmrR.promoteroverlap.gr       <- suppressWarnings((intersect(hsmrR.gr,promoter.gr,ignore.strand=TRUE)))

    hsmrR.In.GenicRegions          <- sum(width(hsmrR.geneoverlap.gr))
    hsmrR.In.ExonRegions           <- sum(width(hsmrR.exonoverlap.gr))
    hsmrR.In.IntronRegions         <- sum(width(hsmrR.geneoverlap.gr))-sum(width(hsmrR.exonoverlap.gr))

    hsmrR.In.FullIntergenicRegions <- sum(width(hsmrR.gr)) - sum(width(hsmrR.geneoverlap.gr))
    hsmrR.promotergeneoverlap.gr   <- setdiff(hsmrR.promoteroverlap.gr,hsmrR.geneoverlap.gr,ignore.strand=TRUE)
    hsmrR.In.PromoterRegions       <- sum(width(hsmrR.promotergeneoverlap.gr))
    hsmrR.In.IntergenicRegions     <- hsmrR.In.FullIntergenicRegions - hsmrR.In.PromoterRegions

    hsmrR.pct          <- round(100*sum(width(hsmrR.gr))/gnmsize,2)
    hsmrR.slices1      <- c(hsmrR.In.GenicRegions,hsmrR.In.FullIntergenicRegions)
    hsmrR.pct1         <- round(100*hsmrR.slices1/sum(hsmrR.slices1),2)
    hsmrR.slices2      <- c(hsmrR.In.ExonRegions,hsmrR.In.IntronRegions,hsmrR.In.PromoterRegions,hsmrR.In.IntergenicRegions)
    hsmrR.pct2         <- round(100*hsmrR.slices2/sum(hsmrR.slices2),2)

    cat( sprintf("Overlap of methylation-rich regions with genome features for %s sample %s\n\n", species, slabel) )
    cat( sprintf("Total size of the methylation-rich regions of %s:                                  %9d bp (%5.1f%% of genome)\n",slabel,sum(width(hsmrR.gr)),hsmrR.pct) )      
    cat( sprintf("Size of overlap of methylation-rich regions of %s %s with genic regions:              %9d bp (%5.1f%% of total size)   (%5.2f O/E)\n",species,slabel,hsmrR.In.GenicRegions,hsmrR.pct1[1],hsmrR.pct1[1]/(100*gene.fraction)) )
    cat( sprintf("Size of overlap of methylation-rich regions of %s %s with   exon regions:               %9d bp (%5.1f%% of total size) (%5.2f O/E)\n",species,slabel,hsmrR.In.ExonRegions,hsmrR.pct2[1],hsmrR.pct2[2]/(100*exon.fraction)) )
    cat( sprintf("Size of overlap of methylation-rich regions of %s %s with   intron regions:             %9d bp (%5.1f%% of total size) (%5.2f O/E)\n",species,slabel,hsmrR.In.IntronRegions,hsmrR.pct2[2],hsmrR.pct2[2]/(100*intron.fraction)) )
    cat( sprintf("Size of overlap of methylation-rich regions of %s %s with intergenic regions:         %9d bp (%5.1f%% of total size)   (%5.2f O/E)\n",species,slabel,hsmrR.In.FullIntergenicRegions,hsmrR.pct1[2],hsmrR.pct1[1]/(100*intergenic.fraction)) )
    cat( sprintf("Size of overlap of methylation-rich regions of %s %s with   promoter regions:           %9d bp (%5.1f%% of total size) (%5.2f O/E)\n",species,slabel,hsmrR.In.PromoterRegions,hsmrR.pct2[3],hsmrR.pct2[3]/(100*promoter.fraction)) )
    cat( sprintf("Size of overlap of methylation-rich regions of %s %s with   other intergenic regions:   %9d bp (%5.1f%% of total size) (%5.2f O/E)\n",species,slabel,hsmrR.In.IntergenicRegions,hsmrR.pct2[4],hsmrR.pct2[4]/(100*intergenicsanspromoter.fraction)) )
    cat( sprintf( "\n\n" ) )

    hsmrRinGenes  <- suppressWarnings(findOverlaps(hsmrR.gr,gene.gr,ignore.strand=TRUE))
    gwith    <- gene.gr[subjectHits(hsmrRinGenes)]
    rwith    <- hsmrR.gr[queryHits(hsmrRinGenes)]

    gwithdf <- as.data.frame(gwith)
    rwithdf <- as.data.frame(rwith)
    gwithrRdf <- cbind(gwithdf,rwithdf)
    gwrfile <- paste("gwr",slabel,sep="-")
    gwrfile <- paste(gwrfile,"txt",sep=".")
    write.table(gwithrRdf, file=gwrfile, sep="\t", row.names=FALSE, quote=FALSE)


    #hsm-poor regions:
    #
    tmpdf <- hsmrL$hsmrP
    hsmpR.gr <- with(tmpdf,{GRanges(tmpdf$SeqID,
                                    IRanges(tmpdf$From,tmpdf$To),
                                    NbrSites=(tmpdf$NbrSites),                                               
                                    Sdnsty=(tmpdf$Sdnsty))
                           }
                           )
    hsmpR.geneoverlap.gr           <- suppressWarnings((intersect(hsmpR.gr,gene.gr,ignore.strand=TRUE)))
    hsmpR.exonoverlap.gr           <- suppressWarnings((intersect(hsmpR.gr,exon.gr,ignore.strand=TRUE)))
    hsmpR.promoteroverlap.gr       <- suppressWarnings((intersect(hsmpR.gr,promoter.gr,ignore.strand=TRUE)))

    hsmpR.In.GenicRegions          <- sum(width(hsmpR.geneoverlap.gr))
    hsmpR.In.ExonRegions           <- sum(width(hsmpR.exonoverlap.gr))
    hsmpR.In.IntronRegions         <- sum(width(hsmpR.geneoverlap.gr))-sum(width(hsmpR.exonoverlap.gr))

    hsmpR.In.FullIntergenicRegions <- sum(width(hsmpR.gr)) - sum(width(hsmpR.geneoverlap.gr))
    hsmpR.promotergeneoverlap.gr   <- setdiff(hsmpR.promoteroverlap.gr,hsmpR.geneoverlap.gr,ignore.strand=TRUE)
    hsmpR.In.PromoterRegions       <- sum(width(hsmpR.promotergeneoverlap.gr))
    hsmpR.In.IntergenicRegions     <- hsmpR.In.FullIntergenicRegions - hsmpR.In.PromoterRegions

    hsmpR.pct          <- round(100*sum(width(hsmpR.gr))/gnmsize,2)
    hsmpR.slices1      <- c(hsmpR.In.GenicRegions,hsmpR.In.FullIntergenicRegions)
    hsmpR.pct1         <- round(100*hsmpR.slices1/sum(hsmpR.slices1),2)
    hsmpR.slices2      <- c(hsmpR.In.ExonRegions,hsmpR.In.IntronRegions,hsmpR.In.PromoterRegions,hsmpR.In.IntergenicRegions)
    hsmpR.pct2         <- round(100*hsmpR.slices2/sum(hsmpR.slices2),2)

    cat( sprintf("Overlap of methylation-poor regions with genome features for %s sample %s\n\n", species, slabel) )
    cat( sprintf("Total size of the methylation-poor regions of %s:                                  %9d bp (%5.1f%% of genome)\n",slabel,sum(width(hsmpR.gr)),hsmpR.pct) )      
    cat( sprintf("Size of overlap of methylation-poor regions of %s %s with genic regions:              %9d bp (%5.1f%% of total size)   (%5.2f O/E)\n",species,slabel,hsmpR.In.GenicRegions,hsmpR.pct1[1],hsmpR.pct1[1]/(100*gene.fraction)) )
    cat( sprintf("Size of overlap of methylation-poor regions of %s %s with   exon regions:               %9d bp (%5.1f%% of total size) (%5.2f O/E)\n",species,slabel,hsmpR.In.ExonRegions,hsmpR.pct2[1],hsmpR.pct2[2]/(100*exon.fraction)) )
    cat( sprintf("Size of overlap of methylation-poor regions of %s %s with   intron regions:             %9d bp (%5.1f%% of total size) (%5.2f O/E)\n",species,slabel,hsmpR.In.IntronRegions,hsmpR.pct2[2],hsmpR.pct2[2]/(100*intron.fraction)) )
    cat( sprintf("Size of overlap of methylation-poor regions of %s %s with intergenic regions:         %9d bp (%5.1f%% of total size)   (%5.2f O/E)\n",species,slabel,hsmpR.In.FullIntergenicRegions,hsmpR.pct1[2],hsmpR.pct1[1]/(100*intergenic.fraction)) )
    cat( sprintf("Size of overlap of methylation-poor regions of %s %s with   promoter regions:           %9d bp (%5.1f%% of total size) (%5.2f O/E)\n",species,slabel,hsmpR.In.PromoterRegions,hsmpR.pct2[3],hsmpR.pct2[3]/(100*promoter.fraction)) )
    cat( sprintf("Size of overlap of methylation-poor regions of %s %s with   other intergenic regions:   %9d bp (%5.1f%% of total size) (%5.2f O/E)\n",species,slabel,hsmpR.In.IntergenicRegions,hsmpR.pct2[4],hsmpR.pct2[4]/(100*intergenicsanspromoter.fraction)) )
    cat( sprintf( "\n\n" ) )

    hsmpRinGenes  <- suppressWarnings(findOverlaps(hsmpR.gr,gene.gr,ignore.strand=TRUE))
    gwith    <- gene.gr[subjectHits(hsmpRinGenes)]
    rwith    <- hsmpR.gr[queryHits(hsmpRinGenes)]

    gwithdf <- as.data.frame(gwith)
    rwithdf <- as.data.frame(rwith)
    gwithpRdf <- cbind(gwithdf,rwithdf)
    gwrfile <- paste("gwp",slabel,sep="-")
    gwrfile <- paste(gwrfile,"txt",sep=".")
    write.table(gwithpRdf, file=gwrfile, sep="\t", row.names=FALSE, quote=FALSE)


    cat( sprintf( "\n\n" ) )
    if (outfile != "") {
        sink()
    }
    return(list(gwithrRdf,gwithpRdf))
}
