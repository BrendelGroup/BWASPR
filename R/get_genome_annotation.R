#' get_genome_annotation()
#' This function derives GRanges objects from the GFF3-formatted input file data.
#' The GFF3 input files should be specified by setup_BWASPR().
#'
#' @param inputdf A data frame as returned by setup_BWASPR from the input data file.
#'   The data frame contains the directories to specific generic features of DNA and the UTRflag.
#'
#' @return A list consisting of GRanges objects that describe generic features of DNA including:
#'   gene,exon,pcexon,promoter,CDS,fpUTR,tpUTR,
#'   fpUTRnotCDS,tpUTRnotCDS,tpUTRunique,ncexon
#'
#' @importFrom genomation gffToGRanges
#' @importFrom GenomicRanges GRanges setdiff
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   genome <- get_genome_annotation(myfiles$parameters)
#'
#' @export


get_genome_annotation <- function(inputdf) {
    # read the directory of genome info from inputdf
    GFF3DIR               <- inputdf[inputdf$Variable == 'SPECIESGFF3DIR', "Value"]
    genelist              <- paste(GFF3DIR,inputdf[inputdf$Variable == 'GENELISTGFF3', "Value"],sep="/")
    exonlist              <- paste(GFF3DIR,inputdf[inputdf$Variable == 'EXONLISTGFF3', "Value"],sep="/")
    proteincodingexonlist <- paste(GFF3DIR,inputdf[inputdf$Variable == 'PCGEXNLISTGFF3', "Value"],sep="/")
    promoterlist          <- paste(GFF3DIR,inputdf[inputdf$Variable == 'PROMOTRLISTGFF3', "Value"],sep="/")
    cdslist               <- paste(GFF3DIR,inputdf[inputdf$Variable == 'CDSLISTGFF3', "Value"],sep="/")
    UTRflag               <- inputdf[inputdf$Variable == 'UTRFLAGSET', "Value"]

    if (UTRflag == 1) {
        fpUTRlist <- paste(GFF3DIR,inputdf[inputdf$Variable == '5UTRLISTGFF3', "Value"],sep="/")
        tpUTRlist <- paste(GFF3DIR,inputdf[inputdf$Variable == '3UTRLISTGFF3', "Value"],sep="/")
    }

    gene.gr     <- gffToGRanges(genelist)
    exon.gr     <- gffToGRanges(exonlist)
    pcexon.gr   <- gffToGRanges(proteincodingexonlist)
    promoter.gr <- gffToGRanges(promoterlist)
    CDS.gr      <- gffToGRanges(cdslist)

    if (UTRflag == 1) {
        fpUTR.gr  <- gffToGRanges(fpUTRlist)
        tpUTR.gr <- gffToGRanges(tpUTRlist)
        # To obtain the five-prime UTR that does not overlap with CDS
        #
        if (length(fpUTR.gr) > 0) {
            fpUTRnotCDS.gr <- suppressWarnings(setdiff(fpUTR.gr,CDS.gr,
                                                               ignore.strand=TRUE
                                                              )
                                                      )
        } else {
            fpUTRnotCDS.gr <- fpUTR.gr;
        }
        # To obtain the three-prime UTR that does not overlap with CDS
        #
        if (length(tpUTR.gr) > 0) {
            tpUTRnotCDS.gr <- suppressWarnings(setdiff(tpUTR.gr,CDS.gr,
                                                               ignore.strand=TRUE
                                                              )
                                                      )
        } else {
            tpUTRnotCDS.gr <- tpUTR.gr;
        }
        if (length(tpUTRnotCDS.gr) > 0) {
            tpUTRunique.gr <- suppressWarnings(setdiff(tpUTRnotCDS.gr,
                                                                fpUTRnotCDS.gr,
                                                                ignore.strand=TRUE
                                                               )
                                                       )
        } else {
            tpUTRunique.gr <- tpUTRnotCDS.gr;
        }
    }
    else {	# return empty GRanges if there are no UTR annotations ..,
        fpUTR.gr <- GRanges()
        tpUTR.gr <- GRanges()
        fpUTRnotCDS.gr <- GRanges()
        tpUTRnotCDS.gr <- GRanges()
        tpUTRunique.gr <- GRanges()
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
                'fpUTR' = fpUTR.gr,
                'tpUTR' = tpUTR.gr,
                'fpUTRnotCDS' = fpUTRnotCDS.gr,
                'tpUTRnotCDS' = tpUTRnotCDS.gr,
                'tpUTRunique' = tpUTRunique.gr,
                'ncexon' = ncexon.gr
               )
          )
}
