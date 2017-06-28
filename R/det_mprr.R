#' det_mprr()
#' This function generates coverage and methylation statistics for the
#'   input data set.
#'
#' @param mrobj A methylRaw object
#' @param covlist a vector of coverage threshols; default: c(10)
#' @param outfstem If specified then output is printed to the specified file.
#' @param doplots Logical; if true, then show plots of distances
#'
#' @return data frame
#'
#' @importFrom methylKit getData getCoverageStats getMethylationStats
#' @importFrom pastecs stat.pen
#' @import ggplot2
#  @import dplyr
#' @importFrom dplyr arrange
#' @import grid
#' @import gridExtra
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                        type="CpGhsm", mincov=1,assembly="Amel-4.5")
#'   det_mprr(AmHE[[1]],ddset=c(1,5),nr2d=10L,outfstem="Am_HE-regions",doplots=TRUE)
#'
#' @export

det_mprr <- function(mrobj,ddset=c(1,5),nr2d=10L,outfstem="",doplots=TRUE) {
# ... set checkflag=1 to print out out loads of detail to stdout ...:
    checkflag <- 0
    sampledata <- getData(mrobj)
    slabel <- mrobj@sample.id
    message("... determining methylation site poor and rich regions for " ,slabel," ...")

    DFhsmr <- data.frame( Rtype = as.character(), SeqId = as.character(), From = as.numeric(), To = as.numeric(), Rlgth = as.numeric(), NbrSites = as.numeric(), Sdnsty = as.numeric(), Fsite = as.numeric(), Tsite = as.numeric(), Dvalue = as.numeric() )
    
    for ( d in ddset ) {
      cat( sprintf("\n") )
      plotme <- list()
    
        xlabel <- sprintf("%s (%2d-distances)",slabel,d)
        sampledata <- getData(mrobj)
      
        cat( sprintf("Analysis of %2d-distances for sample \"%s\":\n\n",d,slabel) )
      
    # gdvector = vector of d-distances over the entire data set:
    #
        gdvector <- c()
        gdbegpnt <- c()
        gdendpnt <- c()
      
        seqid <- sampledata$chr[1]
        nsites= 0
        icnt= 0
    # dvector = vector of d-distances over sequence seqid
    #
        dvector <- c()
      
        for ( i in 1:nrow(sampledata) ) {
          if ( !identical(sampledata$chr[i],seqid) ) {
            if (checkflag) {
              cat( sprintf("%s\thas %4d sites and %4d %2d-distances\n", seqid, nsites, length(dvector), d) )
            }
            seqid <- sampledata$chr[i]
            nsites= 1
            icnt= 1
            gdvector <- c(gdvector,dvector)
            dvector  <- c()
            if ( i == nrow(sampledata)  &&  checkflag ) {
              cat( sprintf("%s\thas %4d sites and %4d %2d-distances\n", sampledata$chr[i], 1, 0, d) )
            }
          } else if ( i == nrow(sampledata) ) {
            icnt <- icnt + 1
            nsites <- nsites + 1
            if ( icnt >= d ) {
              dstnc = sampledata$start[i] - sampledata$start[i-d]
              if (checkflag) {
                cat( sprintf("\t%s\t%9d %9d %d-dstnc:\t%6d\n", sampledata$chr[i], sampledata$start[i-d], sampledata$start[i], d, dstnc) ) 
              }
              dvector <- c(dvector, dstnc)
              gdbegpnt <- c(gdbegpnt, i-d)
              gdendpnt <- c(gdendpnt, i)
            }
            gdvector <- c(gdvector,dvector)
            if (checkflag) {
              cat( sprintf("%s\thas %4d sites and %4d %2d-distances\n", seqid, nsites, length(dvector), d) )
            }
          } else {
            icnt <- icnt + 1
            if ( icnt >  d ) {
              dstnc = sampledata$start[i] - sampledata$start[i-d]
              if (checkflag) {
                cat( sprintf("\t%s\t%9d %9d %d-dstnc:\t%6d\n", sampledata$chr[i], sampledata$start[i-d], sampledata$start[i], d, dstnc) ) 
              }
              dvector <- c(dvector, dstnc)
              gdbegpnt <- c(gdbegpnt, i-d)
              gdendpnt <- c(gdendpnt, i)
            }
            nsites <- nsites + 1
          }
        }
        cat( sprintf("  Total number of %2d-distances:\t%6d\n", d, length(gdvector)) )
        ddlist <- ddstats(gdvector,xlabel,"black")
        dstats <- ddlist$dstats
        plotme[[1]] <- ddlist$plot1
        plotme[[2]] <- ddlist$plot2
        cat( sprintf("  Median: %7.1f Mean: %7.1f Std: %7.1f\n", dstats["median"], dstats["mean"], dstats["std.dev"]) )
        cat( sprintf("  Quantiles:\n") )
        q <- quantile(gdvector,probs=seq(0,1,0.05))
        print( q )
        cat( sprintf("\n") )
        cat( sprintf("  Low  density region distance cutoff (90-percentile): %6.0f\n", q["90%"]) )
        cat( sprintf("  High density region distance cutoff (10-percentile): %6.0f\n", q["10%"]) )
        cat( sprintf("\n") )
    
        dhsmr.name <- ""
        phsmr.name <- ""
    
        for ( i in 1:length(gdvector) ) {
          if ( d > 1  &&  gdvector[i] <= unname(q["10%"]) ) {
            if (dhsmr.name == "") {
              dhsmr.name <- as.character(sampledata$chr[gdbegpnt[i]])
              dhsmr.begi <- gdbegpnt[i]
              dhsmr.begpnt <- sampledata$start[gdbegpnt[i]]
              dhsmr.endi <- gdendpnt[i]
              dhsmr.endpnt <- sampledata$start[gdendpnt[i]]
            } else if (as.character(sampledata$chr[gdbegpnt[i]]) == dhsmr.name  &&  sampledata$start[gdbegpnt[i]] <= dhsmr.endpnt) {
              dhsmr.endi <- gdendpnt[i]
              dhsmr.endpnt <- as.numeric(sampledata$start[gdendpnt[i]])
            } else {
              hsmr.record <- data.frame( Rtype = "Rich", Sample = slabel, SeqID = dhsmr.name, From = dhsmr.begpnt, To = dhsmr.endpnt, Rlgth = dhsmr.endpnt-dhsmr.begpnt+1, NbrSites = dhsmr.endi-dhsmr.begi+1, Sdnsty = 1000.*(dhsmr.endi-dhsmr.begi+1)/(dhsmr.endpnt-dhsmr.begpnt), Fsite = dhsmr.begi, Tsite = dhsmr.endi, Dvalue = d )
              DFhsmr <- rbind(DFhsmr,hsmr.record)
              if (checkflag) {
                cat( sprintf("Rich \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\t:number of sites, %4d, per-kb-density: %6.2f\n", slabel, dhsmr.name, dhsmr.begpnt, dhsmr.begi, dhsmr.endpnt, dhsmr.endi, d, dhsmr.endpnt-dhsmr.begpnt, dhsmr.endi-dhsmr.begi+1, 1000.*(dhsmr.endi-dhsmr.begi+1)/(dhsmr.endpnt-dhsmr.begpnt) ) )
                for ( j in dhsmr.begi:dhsmr.endi ) {
                  print( sampledata[j,] )
                }
              }
              dhsmr.name <- as.character(sampledata$chr[gdbegpnt[i]])
              dhsmr.begi <- gdbegpnt[i]
              dhsmr.begpnt <- sampledata$start[gdbegpnt[i]]
              dhsmr.endi <- gdendpnt[i]
              dhsmr.endpnt <- sampledata$start[gdendpnt[i]]
            }
            if (checkflag) {
              cat( sprintf("Rich hsm region %s from %d to %d, %2d-distance:\t%6d\n",sampledata$chr[gdbegpnt[i]], sampledata$start[gdbegpnt[i]], sampledata$start[gdendpnt[i]], d, sampledata$start[gdendpnt[i]] - sampledata$start[gdbegpnt[i]] ) )
              for ( j in gdbegpnt[i]:gdendpnt[i] ) {
                print( sampledata[j,] )
              }
            }
          } else if ( d == 1  &&  gdvector[i] > unname(q["90%"]) ) {
            if (phsmr.name == "") {
              phsmr.name <- as.character(sampledata$chr[gdbegpnt[i]])
              phsmr.begi <- gdbegpnt[i]
              phsmr.begpnt <- sampledata$start[gdbegpnt[i]]
              phsmr.endi <- gdendpnt[i]
              phsmr.endpnt <- sampledata$start[gdendpnt[i]]
            } else if (as.character(sampledata$chr[gdbegpnt[i]]) == phsmr.name  &&  sampledata$start[gdbegpnt[i]] <= phsmr.endpnt) {
              phsmr.endi <- gdendpnt[i]
              phsmr.endpnt <- as.numeric(sampledata$start[gdendpnt[i]])
            } else {
              hsmr.record <- data.frame( Rtype = "Poor", Sample = slabel, SeqID = phsmr.name, From = phsmr.begpnt, To = phsmr.endpnt, Rlgth = phsmr.endpnt-phsmr.begpnt+1, NbrSites = phsmr.endi-phsmr.begi+1, Sdnsty = 1000.*(phsmr.endi-phsmr.begi+1)/(phsmr.endpnt-phsmr.begpnt), Fsite = phsmr.begi, Tsite = phsmr.endi, Dvalue = d )
              DFhsmr <- rbind(DFhsmr,hsmr.record)
              if (checkflag) {
                cat( sprintf("Poor \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\t:number of sites, %4d, per-kb-density: %6.2f\n", slabel, phsmr.name, phsmr.begpnt, phsmr.begi, phsmr.endpnt, phsmr.endi, d, phsmr.endpnt-phsmr.begpnt, phsmr.endi-phsmr.begi+1, 1000.*(phsmr.endi-phsmr.begi+1)/(phsmr.endpnt-phsmr.begpnt) ) )
                for ( j in phsmr.begi:phsmr.endi ) {
                  print( sampledata[j,] )
                }
              }
              phsmr.name <- as.character(sampledata$chr[gdbegpnt[i]])
              phsmr.begi <- gdbegpnt[i]
              phsmr.begpnt <- sampledata$start[gdbegpnt[i]]
              phsmr.endi <- gdendpnt[i]
              phsmr.endpnt <- sampledata$start[gdendpnt[i]]
            }
            if (checkflag) {
              cat( sprintf("Poor hsm region %s from %d to %d, %2d-distance:\t%6d\n",sampledata$chr[gdbegpnt[i]], sampledata$start[gdbegpnt[i]], sampledata$start[gdendpnt[i]], d, sampledata$start[gdendpnt[i]] - sampledata$start[gdbegpnt[i]] ) )
              for ( j in gdbegpnt[i]:gdendpnt[i] ) {
                print( sampledata[j,] )
              }
            }
          }
        }
      if ( d > 1 ) {
        hsmr.record <- data.frame( Rtype = "Rich", Sample = slabel, SeqID = dhsmr.name, From = dhsmr.begpnt, To = dhsmr.endpnt, Rlgth = dhsmr.endpnt-dhsmr.begpnt+1, NbrSites = dhsmr.endi-dhsmr.begi+1, Sdnsty = 1000.*(dhsmr.endi-dhsmr.begi+1)/(dhsmr.endpnt-dhsmr.begpnt), Fsite = dhsmr.begi, Tsite = dhsmr.endi, Dvalue = d )
        DFhsmr <- rbind(DFhsmr,hsmr.record)
        if (checkflag) {
          cat( sprintf("Rich \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\t:number of sites, %4d, per-kb-density: %6.2f\n", slabel, dhsmr.name, dhsmr.begpnt, dhsmr.begi, dhsmr.endpnt, dhsmr.endi, d, dhsmr.endpnt-dhsmr.begpnt, dhsmr.endi-dhsmr.begi+1, 1000.*(dhsmr.endi-dhsmr.begi+1)/(dhsmr.endpnt-dhsmr.begpnt) ) )
          for ( j in dhsmr.begi:dhsmr.endi ) {
            print( sampledata[j,] )
          }
        }
      }
      if ( d == 1 ) {
        hsmr.record <- data.frame( Rtype = "Poor", Sample = slabel, SeqID = phsmr.name, From = phsmr.begpnt, To = phsmr.endpnt, Rlgth = phsmr.endpnt-phsmr.begpnt+1, NbrSites = phsmr.endi-phsmr.begi+1, Sdnsty = 1000.*(phsmr.endi-phsmr.begi+1)/(phsmr.endpnt-phsmr.begpnt), Fsite = phsmr.begi, Tsite = phsmr.endi, Dvalue = d )
        DFhsmr <- rbind(DFhsmr,hsmr.record)
        if (checkflag) {
          cat( sprintf("Poor \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\t:number of sites, %4d, per-kb-density: %6.2f\n", slabel, phsmr.name, phsmr.begpnt, phsmr.begi, phsmr.endpnt, phsmr.endi, d, phsmr.endpnt-phsmr.begpnt, phsmr.endi-phsmr.begi+1, 1000.*(phsmr.endi-phsmr.begi+1)/(phsmr.endpnt-phsmr.begpnt) ) )
          for ( j in phsmr.begi:phsmr.endi ) {
            print( sampledata[j,] )
          }
        }
      }
    
      if (d == 1) {
        cat( sprintf("\nOrdered lists of methylation-poor regions:\n" ) )
      } else {
        cat( sprintf("\nOrdered lists of methylation-rich regions:\n" ) )
      }
      if (checkflag) {
        print( DFhsmr[order(DFhsmr$Sdnsty),] )
      }
    
      if (d > 1) {
        DFhsmrR <- subset(DFhsmr, Rtype == "Rich" & Sample == slabel, select = c(Rtype,Sample,SeqID,From,To,Rlgth,NbrSites,Sdnsty))
        DFhsmrR <- arrange(DFhsmrR, -Sdnsty)
        slabel.DFhsmrR.df <- as.data.frame(DFhsmrR)
        write.table(slabel.DFhsmrR.df, sprintf("mrr-%s_%s.txt",outfstem,slabel), quote=F, sep="\t", row.names=F, col.names=F)
        cat( sprintf("\n Top %3d (out of %5d) methylation-rich regions in \"%s\":\n", nr2d, dim(DFhsmrR)[1], slabel) )
        print( head(DFhsmrR, n = nr2d) )
      }
    
      if (d == 1) {
        DFhsmrP <- subset(DFhsmr, Rtype == "Poor" & Sample == slabel, select = c(Rtype,Sample,SeqID,From,To,Rlgth,NbrSites,Sdnsty))
        DFhsmrP <- arrange(DFhsmrP, +Sdnsty)
        slabel.DFhsmrP.df <- as.data.frame(DFhsmrP)
        write.table(slabel.DFhsmrP.df, sprintf("mpr-%s_%s.txt",outfstem,slabel), quote=F, sep="\t", row.names=F, col.names=F)
        cat( sprintf("\n Top %3d (out of %5d) methylation-poor regions in \"%s\":\n", nr2d, dim(DFhsmrP)[1], slabel) )
        print( head(DFhsmrP, n = nr2d) )
      }
      cat( sprintf("\n\n") )
    
      if (doplots) {
        pdffile <- sprintf("%dds-%s.pdf",d,outfstem)
        mytitle <- sprintf( "\n%s\n%2d-distance distributions", slabel, d)
        ggsave(pdffile, do.call(marrangeGrob, list(grobs=plotme, nrow=1, ncol=2, top = mytitle)), width = 7, height = 7)
        graphics.off()
      }
    }
    
    DFhsmrPRT <- data.frame( DFhsmr$SeqID, DFhsmr$From, DFhsmr$To, DFhsmr$Rtype, DFhsmr$Sample, DFhsmr$Rlgth, DFhsmr$NbrSites, DFhsmr$Sdnsty)
    names(DFhsmrPRT) <- c("SeqID", "From", "To", "Rtype", "Sample", "Rlgth", "NbrSites", "Sdnsty")
    
    tabfile <- sprintf("mdr-%s.tab",outfstem)
    bedfile <- sprintf("mdr-%s.bed",outfstem)
    sink(tabfile)
    write.table( format(DFhsmrPRT,scientific=FALSE,digits=1L), row.names = FALSE, quote = FALSE, sep="\t" )
    sink()
    cmd <- sprintf("sed -e 's/ //g' %s > %s",tabfile,bedfile)
    system(cmd)
    
    cat( sprintf( "\n\n" ) )
}


# function ddstats - calculates and plots statistics for d-distances
# 
#
ddstats <- function(v,xlabel,mycolor) {
  checkflag <- 0
  if (length(v) < 2) {
    return()
  }
  dstats <- stat.pen(v,desc=T)
  if (checkflag) {
    print( dstats )
  }
  dstncdf <- as.data.frame(v)
  pme <- subset(dstncdf,v<=600)
  plot1 <-  ggplot(pme, aes(x=v)) + scale_x_continuous(xlabel,limits=c(0,600),breaks=seq(0,600,60)) + geom_histogram(aes(y=..density..),binwidth=20,closed="right",color="black",fill="white") + scale_y_continuous(limits=c(0,0.003)) + geom_density(alpha=.2,color=mycolor)
  pme <- subset(dstncdf,v<=60)
  plot2 <-  ggplot(pme, aes(x=v)) + scale_x_continuous(xlabel,limits=c(0,60),breaks=seq(0,60,4)) + geom_histogram(aes(y=..density..),binwidth=2,closed="right",color="black",fill="white") + scale_y_continuous(limits=c(0,0.03)) + geom_density(alpha=.2,color=mycolor)
  rlist <- list("dstats" = dstats, "plot1" = plot1, "plot2" = plot2)

  return(rlist)
}
