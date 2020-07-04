#' det_mrpr()
#'   This function determines methylation-rich and -poor regions based on clustering
#'   of methylation (hsm sites). First, the ddstats() function is used to establish
#'   the distribution of distances between d-nearest neighbors.  Short distances
#'   indcate clustering of hsm sites, and long distances indicate regions relatively
#'   devoid of hsm sites.
#'
#' @param mrobj A methylRaw object.
#' @param sampleL Label used in output to identify the sample.
#' @param ddset Vector of distances; d=1 is required to detect methylation-poor
#'   regions; default: c(1,5)
#' @param outfile If specified, then output is saved in the specified file name.
#' @param nr2d integer; number of regions to display
#' @param doplots Logical; if true, then show plots of distances
#'
#' @return List of two data frames, recording hsm rich and poor regions.
#'
#' @importFrom methylKit getData getCoverageStats getMethylationStats
#' @importFrom pastecs stat.pen
#' @import ggplot2 R.devices
#' @import gridExtra
#'
#' @examples
#'   mydatf <- system.file("extdata","Am.dat",package="BWASPR")
#'   myparf <- system.file("extdata","Am.par",package="BWASPR")
#'   myfiles <- setup_BWASPR(datafile=mydatf,parfile=myparf)
#'   AmHE <- mcalls2mkobj(myfiles$datafiles,species="Am",study="HE",
#'                        type="CpGhsm", mincov=1,assembly="Amel-4.5")
#'   det_mrpr(AmHE[[1]],"Am_HE_fr",ddset=c(1,5),outfile="dst-Am_HE_fr.txt",nr2d=10L,doplots=TRUE)
#'
#' @export

det_mrpr <- function(mrobj,sampleL,ddset=c(1,5),outfile="",nr2d=10L,doplots=TRUE) {
    message("... determining methylation site rich and poor regions for " ,sampleL," ...")
# ... set checkflag=1 to print out loads of detail to stdout ...:
    checkflag <- 0
    sampledata <- getData(mrobj)
    sampleID <- mrobj@sample.id
    if (outfile != "") {
        sink(outfile)
    }

    DFhsmr <- data.frame( Rtype = as.character(), SeqId = as.character(), From = as.numeric(), To = as.numeric(), Rlgth = as.numeric(), NbrSites = as.numeric(), Sdnsty = as.numeric(), Fsite = as.numeric(), Tsite = as.numeric(), Dvalue = as.numeric() )
    
    for ( d in ddset ) {
      cat( sprintf("\n") )
      plotme <- list()
    
      xlabel <- sprintf("%s (%2d-distances)",sampleL,d)
      sampledata <- getData(mrobj)
      
      cat( sprintf("Analysis of %2d-distances for sample \"%s\":\n\n",d,sampleL) )
      
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
            hsmr.record <- data.frame( Rtype = "Rich", Sample = sampleID, SeqID = dhsmr.name, From = dhsmr.begpnt, To = dhsmr.endpnt, Rlgth = dhsmr.endpnt-dhsmr.begpnt+1, NbrSites = dhsmr.endi-dhsmr.begi+1, Sdnsty = 1000.*(dhsmr.endi-dhsmr.begi+1)/(dhsmr.endpnt-dhsmr.begpnt+1), Fsite = dhsmr.begi, Tsite = dhsmr.endi, Dvalue = d )
            DFhsmr <- rbind(DFhsmr,hsmr.record)
            if (checkflag) {
              cat( sprintf("Rich \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\t:number of sites, %4d, per-kb-density: %6.2f\n", sampleL, dhsmr.name, dhsmr.begpnt, dhsmr.begi, dhsmr.endpnt, dhsmr.endi, d, dhsmr.endpnt-dhsmr.begpnt+1, dhsmr.endi-dhsmr.begi+1, 1000.*(dhsmr.endi-dhsmr.begi+1)/(dhsmr.endpnt-dhsmr.begpnt+1) ) )
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
            hsmr.record <- data.frame( Rtype = "Poor", Sample = sampleID, SeqID = phsmr.name, From = phsmr.begpnt, To = phsmr.endpnt, Rlgth = phsmr.endpnt-phsmr.begpnt+1, NbrSites = phsmr.endi-phsmr.begi+1, Sdnsty = 1000.*(phsmr.endi-phsmr.begi+1)/(phsmr.endpnt-phsmr.begpnt+1), Fsite = phsmr.begi, Tsite = phsmr.endi, Dvalue = d )
            DFhsmr <- rbind(DFhsmr,hsmr.record)
            if (checkflag) {
              cat( sprintf("Poor \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\t:number of sites, %4d, per-kb-density: %6.2f\n", sampleL, phsmr.name, phsmr.begpnt, phsmr.begi, phsmr.endpnt, phsmr.endi, d, phsmr.endpnt-phsmr.begpnt+1, phsmr.endi-phsmr.begi+1, 1000.*(phsmr.endi-phsmr.begi+1)/(phsmr.endpnt-phsmr.begpnt+1) ) )
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
            cat( sprintf("Poor hsm region %s from %d to %d, %2d-distance:\t%6d\n",sampledata$chr[gdbegpnt[i]], sampledata$start[gdbegpnt[i]], sampledata$start[gdendpnt[i]], d, sampledata$start[gdendpnt[i]]-sampledata$start[gdbegpnt[i]]+1) )
            for ( j in gdbegpnt[i]:gdendpnt[i] ) {
              print( sampledata[j,] )
            }
          }
        }
      }
      if ( d > 1 ) {
        hsmr.record <- data.frame( Rtype = "Rich", Sample = sampleID, SeqID = dhsmr.name, From = dhsmr.begpnt, To = dhsmr.endpnt, Rlgth = dhsmr.endpnt-dhsmr.begpnt+1, NbrSites = dhsmr.endi-dhsmr.begi+1, Sdnsty = 1000.*(dhsmr.endi-dhsmr.begi+1)/(dhsmr.endpnt-dhsmr.begpnt), Fsite = dhsmr.begi, Tsite = dhsmr.endi, Dvalue = d )
        DFhsmr <- rbind(DFhsmr,hsmr.record)
        if (checkflag) {
          cat( sprintf("Rich \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\t:number of sites, %4d, per-kb-density: %6.2f\n", sampleL, dhsmr.name, dhsmr.begpnt, dhsmr.begi, dhsmr.endpnt, dhsmr.endi, d, dhsmr.endpnt-dhsmr.begpnt+1, dhsmr.endi-dhsmr.begi+1, 1000.*(dhsmr.endi-dhsmr.begi+1)/(dhsmr.endpnt-dhsmr.begpnt+1) ) )
          for ( j in dhsmr.begi:dhsmr.endi ) {
            print( sampledata[j,] )
          }
        }
      }
      if ( d == 1 ) {
        hsmr.record <- data.frame( Rtype = "Poor", Sample = sampleID, SeqID = phsmr.name, From = phsmr.begpnt, To = phsmr.endpnt, Rlgth = phsmr.endpnt-phsmr.begpnt+1, NbrSites = phsmr.endi-phsmr.begi+1, Sdnsty = 1000.*(phsmr.endi-phsmr.begi+1)/(phsmr.endpnt-phsmr.begpnt), Fsite = phsmr.begi, Tsite = phsmr.endi, Dvalue = d )
        DFhsmr <- rbind(DFhsmr,hsmr.record)
        if (checkflag) {
          cat( sprintf("Poor \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\t:number of sites, %4d, per-kb-density: %6.2f\n", sampleL, phsmr.name, phsmr.begpnt, phsmr.begi, phsmr.endpnt, phsmr.endi, d, phsmr.endpnt-phsmr.begpnt+1, phsmr.endi-phsmr.begi+1, 1000.*(phsmr.endi-phsmr.begi+1)/(phsmr.endpnt-phsmr.begpnt+1) ) )
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
        DFhsmrR <- subset(DFhsmr, Rtype == "Rich" & Sample == sampleID, select = c(Rtype,Sample,SeqID,From,To,Rlgth,NbrSites,Sdnsty))
        DFhsmrR <- arrange(DFhsmrR, -Sdnsty)
        sampleL.DFhsmrR.df <- as.data.frame(DFhsmrR)
        write.table(sampleL.DFhsmrR.df, sprintf("mrr-%s.txt",sampleL), quote=F, sep="\t", row.names=F, col.names=F)
        cat( sprintf("\n Top %3d (out of %5d) methylation-rich regions in \"%s\":\n", nr2d, dim(DFhsmrR)[1], sampleL) )
        print( head(DFhsmrR, n = nr2d) )
      }
    
      if (d == 1) {
        DFhsmrP <- subset(DFhsmr, Rtype == "Poor" & Sample == sampleID, select = c(Rtype,Sample,SeqID,From,To,Rlgth,NbrSites,Sdnsty))
        DFhsmrP <- arrange(DFhsmrP, +Sdnsty)
        sampleL.DFhsmrP.df <- as.data.frame(DFhsmrP)
        write.table(sampleL.DFhsmrP.df, sprintf("mpr-%s.txt",sampleL), quote=F, sep="\t", row.names=F, col.names=F)
        cat( sprintf("\n Top %3d (out of %5d) methylation-poor regions in \"%s\":\n", nr2d, dim(DFhsmrP)[1], sampleL) )
        print( head(DFhsmrP, n = nr2d) )
      }
      cat( sprintf("\n\n") )
    
      if (doplots) {
        pdffile <- sprintf("%dds-%s.pdf",d,sampleL)
        mytitle <- sprintf( "\n%s\n%2d-distance distributions", sampleL, d)
        suppressGraphics(ggsave(pdffile, do.call(marrangeGrob, list(grobs=plotme, nrow=1, ncol=2, top = mytitle)), width=7, height=7, units="in"))
      }
    }
    
    DFhsmrPRT <- data.frame( DFhsmr$SeqID, DFhsmr$From, DFhsmr$To, DFhsmr$Rtype, DFhsmr$Sample, DFhsmr$Rlgth, DFhsmr$NbrSites, DFhsmr$Sdnsty)
    names(DFhsmrPRT) <- c("SeqID", "From", "To", "Rtype", "Sample", "Rlgth", "NbrSites", "Sdnsty")
    
    tabfile <- sprintf("mdr-%s.tab",sampleL)
    bedfile <- sprintf("mdr-%s.bed",sampleL)
    sink(tabfile)
    write.table( format(DFhsmrPRT,scientific=FALSE,digits=1L), row.names = FALSE, quote = FALSE, sep="\t" )
    sink()
    cmd <- sprintf("sed -e 's/ //g' %s > %s",tabfile,bedfile)
    system(cmd)

    cat( sprintf( "\n\n" ) )
    sink()
    message('... det_mrpr() finished ...')
    return(list('hsmrR' = DFhsmrR,
                'hsmrP' = DFhsmrP))
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
  plot1 <-  ggplot(pme, aes(x=v)) + scale_x_continuous(xlabel,limits=c(0,600),breaks=seq(0,600,60)) + geom_histogram(aes(y=..density..),binwidth=20,closed="right",color="black",fill="white") + scale_y_continuous(limits=c(0,0.005)) + geom_density(alpha=.2,color=mycolor)
  pme <- subset(dstncdf,v<=60)
  plot2 <-  ggplot(pme, aes(x=v)) + scale_x_continuous(xlabel,limits=c(0,60),breaks=seq(0,60,4)) + geom_histogram(aes(y=..density..),binwidth=2,closed="right",color="black",fill="white") + scale_y_continuous(limits=c(0,0.05)) + geom_density(alpha=.2,color=mycolor)
  rlist <- list("dstats" = dstats, "plot1" = plot1, "plot2" = plot2)

  return(rlist)
}
