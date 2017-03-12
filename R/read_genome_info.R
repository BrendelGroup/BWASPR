#' read_genome_info()
#' This function will read the genome information
#'
#' @param inputdf A data frame as returned by setup_BWASPR from the input
#'   data file.  Each row corresponds to a BWASP-generated output, with columns
#'   indicating Species, Study, Sample, Replicate, Type, and File.
#'
#' @return A list consisting of genome information including
                             #' gene.gr,
                             #' exon.gr,
                             #' pcexon.gr,
                             #' promoter.gr,
                             #' CDS.gr,
                             #' fiveprimeUTR.gr,
                             #' threeprimeUTR.gr,
                             #' fiveprimeUTR_unique.gr,
                             #' threeprimeUTRnotCDS.gr,
                             #' threeprimeUTR_unique.gr,
                             #' ncexon.grl
#'
#' @importFrom GRanges
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   genome <- read_genome_info(myfiles$parameters)
#'


# reading genome information
read_genome_info <- function(inputdf){
    # read the directory of genome info from inputdf
    GFF3DIR               <- inputdf[inputdf$Variable == 'SPECIESGFF3DIR', "Value"]
    # testing:
    GFF3DIR <- '../Am/genome/GFF3DIR/'
    genelist              <- paste(GFF3DIR,inputdf[inputdf$Variable == 'GENELISTGFF3', "Value"],sep="/")
    exonlist              <- paste(GFF3DIR,inputdf[inputdf$Variable == 'EXONLISTGFF3', "Value"],sep="/")
    proteincodingexonlist <- paste(GFF3DIR,inputdf[inputdf$Variable == 'PCGEXNLISTGFF3', "Value"],sep="/")
    promoterlist          <- paste(GFF3DIR,inputdf[inputdf$Variable == 'PROMOTRLISTGFF3', "Value"],sep="/")
    cdslist               <- paste(GFF3DIR,inputdf[inputdf$Variable == 'CDSLISTGFF3', "Value"],sep="/")
    UTRflag               <- inputdf[inputdf$Variable == 'UTRFLAGSET', "Value"]
    if (UTRflag == 1) {
      fiveprimeUTRlist    <- paste(GFF3DIR,inputdf[inputdf$Variable == '5UTRLISTGFF3', "Value"],sep="/")
      threeprimeUTRlist   <- paste(GFF3DIR,inputdf[inputdf$Variable == '3UTRLISTGFF3', "Value"],sep="/")
    }

    # Read genome information
    #
    gene <- read.table(genelist, sep="\t", quote="", fill=FALSE, strip.white=TRUE)

    gene.SeqID  <- gene$V1
    gene.start  <- gene$V4
    gene.end    <- gene$V5
    gene.strand <- gene$V7
    gene.ID     <- gene$V9

    # and create a GRanges object:
    #
    gene.gr <- with(gene,{
       GRanges(gene.SeqID,
               IRanges(gene.start,gene.end),
               GC= (gene.ID),
               strand=(gene.strand)
              )
       }
    )

    # Read the exon file ...
    #
    exon <- read.table(exonlist, sep="\t", quote="", fill=FALSE, strip.white=TRUE)

    exon.SeqID  <- exon$V1
    exon.start  <- exon$V4
    exon.end    <- exon$V5
    exon.strand <- exon$V7
    exon.ID     <- exon$V9

    # and create a GRanges object:
    #
    exon.gr <- with(exon,{
       GRanges(exon.SeqID,
               IRanges(exon.start,exon.end),
               GC = (exon.ID),
               strand=(exon.strand)
              )
       }
    )

    #Read the protein coding exon file ...

    pcexon <-read.table(proteincodingexonlist, sep="\t", quote="", fill=FALSE, strip.white=TRUE)

    pcexon.SeqID  <- pcexon$V1
    pcexon.start  <- pcexon$V4
    pcexon.end    <- pcexon$V5
    pcexon.strand <- pcexon$V7
    pcexon.ID     <- pcexon$V9


    # and create a GRanges object:
    #

    pcexon.gr <- with(pcexon,{
       GRanges(pcexon.SeqID,
               IRanges(pcexon.start,pcexon.end),
               GC = (pcexon.ID),
               strand=(pcexon.strand)
              )
       }
    )


    # Read the promoter file ...
    #
    promoter <- read.table(promoterlist, sep="\t", quote="", fill=FALSE, strip.white=TRUE)

    promoter.SeqID  <- promoter$V1
    promoter.start  <- promoter$V4
    promoter.end    <- promoter$V5
    promoter.strand <- promoter$V7
    promoter.ID     <- promoter$V9

    # and create a GRanges object:
    #
    promoter.gr <- with(promoter,{
       GRanges(promoter.SeqID,
               IRanges(promoter.start,promoter.end),
               GC = (promoter.ID),
               strand=(promoter.strand)
              )
       }
    )

    # Read the CDS file ...
    #

    cds <- read.table(cdslist, sep="\t", quote="", fill=FALSE, strip.white=TRUE)

    cds.SeqID  <- cds$V1
    cds.start  <- cds$V4
    cds.end    <- cds$V5
    cds.strand <- cds$V7
    cds.ID     <- cds$V9

    # and create a GRanges object:
    #
    CDS.gr <- with(cds,{
       GRanges(cds.SeqID,
               IRanges(cds.start,cds.end),
               GC = (cds.ID),
               strand=(cds.strand)
              )
       }
    )

    if (UTRflag == 1) {

    # Read the 5'UTR file ...
    #

    fiveprimeUTR <- read.table(fiveprimeUTRlist, sep="\t", quote="", fill=FALSE, strip.white=TRUE)

    fiveprimeUTR.SeqID  <- fiveprimeUTR$V1
    fiveprimeUTR.start  <- fiveprimeUTR$V4
    fiveprimeUTR.end    <- fiveprimeUTR$V5
    fiveprimeUTR.strand <- fiveprimeUTR$V7
    fiveprimeUTR.ID     <- fiveprimeUTR$V9

    # and create a GRanges object:
    #
    fiveprimeUTR.gr <- with(fiveprimeUTR,{
       GRanges(fiveprimeUTR.SeqID,
               IRanges(fiveprimeUTR.start,fiveprimeUTR.end),
               GC = (fiveprimeUTR.ID),
               strand=(fiveprimeUTR.strand)
              )
       }
    )


    # Read the 3'UTR file ...
    #

    threeprimeUTR <- read.table(threeprimeUTRlist, sep="\t", quote="", fill=FALSE, strip.white=TRUE)

    threeprimeUTR.SeqID  <- threeprimeUTR$V1
    threeprimeUTR.start  <- threeprimeUTR$V4
    threeprimeUTR.end    <- threeprimeUTR$V5
    threeprimeUTR.strand <- threeprimeUTR$V7
    threeprimeUTR.ID     <- threeprimeUTR$V9

    # and create a GRanges object:
    #
    threeprimeUTR.gr <- with(threeprimeUTR,{
       GRanges(threeprimeUTR.SeqID,
               IRanges(threeprimeUTR.start,threeprimeUTR.end),
               GC = (threeprimeUTR.ID),
               strand=(threeprimeUTR.strand)
              )
       }
    )


    # To obtain the five-prime UTR that does not overlap with CDS
    #
    #
    if (length(fiveprimeUTR) > 0) {
      fiveprimeUTR_unique.gr <- suppressWarnings(GenomicRanges::setdiff(fiveprimeUTR.gr,CDS.gr, ignore.strand = TRUE))
      } else {
      fiveprimeUTR_unique.gr <- fiveprimeUTR.gr;
      }

    # # To obtain the three-prime UTR that does not overlap with CDS
    # #
    if (length(threeprimeUTR.gr) > 0) {
      threeprimeUTRnotCDS.gr <- suppressWarnings(GenomicRanges::setdiff(threeprimeUTR.gr,CDS.gr, ignore.strand = TRUE))
      } else {
      threeprimeUTRnotCDS.gr <- threeprimeUTR.gr;
      }

    if (length(threeprimeUTRnotCDS.gr) > 0) {
      threeprimeUTR_unique.gr <- suppressWarnings(GenomicRanges::setdiff(threeprimeUTRnotCDS.gr,fiveprimeUTR_unique.gr, ignore.strand = TRUE))
      } else {
      threeprimeUTR_unique.gr <- threeprimeUTRnotCDS.gr;
    }
    #
    }
    #
    # # Calculate the non-coding regions of the exons
    # #
    #
    if (length(exon.gr) > 0) {
      ncexon.gr <- suppressWarnings(GenomicRanges::setdiff(exon.gr,pcexon.gr, ignore.strand = TRUE))
      } else {
      ncexon.gr <- exon.gr
    }

    return(list('gene' = gene.gr,
                'exon' = exon.gr,
                'pcexon' = pcexon.gr,
                'promoter' = promoter.gr,
                'CDS' = CDS.gr,
                'fiveprimeUTR' = fiveprimeUTR.gr,
                'threeprimeUTR' = threeprimeUTR.gr,
                'fiveprimeUTR' = fiveprimeUTR_unique.gr,
                'threeprimeUTRnotCDS' = threeprimeUTRnotCDS.gr,
                'threeprimeUTR_unique' = threeprimeUTR_unique.gr,
                'ncexon' = ncexon.gr
                ))
}
