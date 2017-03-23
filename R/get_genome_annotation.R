#' get_genome_annotation()
#' This function derives GRanges objects from the GFF3-formatted input file data.
#' The GFF3 input files should be specified by setup_BWASPR().
#'
#' @param inputdf A data frame as returned by setup_BWASPR from the input data file.
#'   The data frame contains the directories to specific generic features of DNA and the UTRflag.
#'
#' @return A list consisting of GRanges objects that describe generic features of DNA including:
#'   gene,exon,pcexon,promoter,CDS,fiveprimeUTR,threeprimeUTR,
#'   fiveprimeUTRunique,threeprimeUTRnotCDS,threeprimeUTRunique,ncexon
#'
#' @importFrom genomation gffToGRanges
#' @importFrom GenomicRanges setdiff
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   genome <- get_genome_annotation(myfiles$parameters)
#'
#' @export


get_genome_annotation <- function(inputdf){
    # read the directory of genome info from inputdf
    GFF3DIR               <- inputdf[inputdf$Variable == 'SPECIESGFF3DIR', "Value"]
    genelist              <- paste(GFF3DIR,inputdf[inputdf$Variable == 'GENELISTGFF3', "Value"],sep="/")
    exonlist              <- paste(GFF3DIR,inputdf[inputdf$Variable == 'EXONLISTGFF3', "Value"],sep="/")
    proteincodingexonlist <- paste(GFF3DIR,inputdf[inputdf$Variable == 'PCGEXNLISTGFF3', "Value"],sep="/")
    promoterlist          <- paste(GFF3DIR,inputdf[inputdf$Variable == 'PROMOTRLISTGFF3', "Value"],sep="/")
    cdslist               <- paste(GFF3DIR,inputdf[inputdf$Variable == 'CDSLISTGFF3', "Value"],sep="/")
    UTRflag               <- inputdf[inputdf$Variable == 'UTRFLAGSET', "Value"]

    if (UTRflag == 1) {
        fiveprimeUTRlist  <- paste(GFF3DIR,inputdf[inputdf$Variable == '5UTRLISTGFF3', "Value"],sep="/")
        threeprimeUTRlist <- paste(GFF3DIR,inputdf[inputdf$Variable == '3UTRLISTGFF3', "Value"],sep="/")
    }

    gene.gr     <- gffToGRanges(genelist)
    exon.gr     <- gffToGRanges(exonlist)
    pcexon.gr   <- gffToGRanges(proteincodingexonlist)
    promoter.gr <- gffToGRanges(promoterlist)
    CDS.gr      <- gffToGRanges(cdslist)

    if (UTRflag == 1) {
        fiveprimeUTR.gr  <- gffToGRanges(fiveprimeUTRlist)
        threeprimeUTR.gr <- gffToGRanges(threeprimeUTRlist)
        # To obtain the five-prime UTR that does not overlap with CDS
        #
        if (length(fiveprimeUTR.gr) > 0) {
            fiveprimeUTRunique.gr <- suppressWarnings(setdiff(fiveprimeUTR.gr,CDS.gr,
                                                               ignore.strand=TRUE
                                                              )
                                                      )
        } else {
            fiveprimeUTRunique.gr <- fiveprimeUTR.gr;
        }
        # To obtain the three-prime UTR that does not overlap with CDS
        #
        if (length(threeprimeUTR.gr) > 0) {
            threeprimeUTRnotCDS.gr <- suppressWarnings(setdiff(threeprimeUTR.gr,CDS.gr,
                                                               ignore.strand=TRUE
                                                              )
                                                      )
        } else {
            threeprimeUTRnotCDS.gr <- threeprimeUTR.gr;
        }
        if (length(threeprimeUTRnotCDS.gr) > 0) {
            threeprimeUTRunique.gr <- suppressWarnings(setdiff(threeprimeUTRnotCDS.gr,
                                                                fiveprimeUTRunique.gr,
                                                                ignore.strand=TRUE
                                                               )
                                                       )
        } else {
            threeprimeUTRunique.gr <- threeprimeUTRnotCDS.gr;
        }
    }

    # Calculate the non-coding regions of the exons
    if (length(exon.gr) > 0) {
        ncexon.gr <- suppressWarnings(setdiff(exon.gr,pcexon.gr, ignore.strand=TRUE))
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
                'fiveprimeUTRunique' = fiveprimeUTRunique.gr,
                'threeprimeUTRnotCDS' = threeprimeUTRnotCDS.gr,
                'threeprimeUTRunique' = threeprimeUTRunique.gr,
                'ncexon' = ncexon.gr
               )
          )
}
