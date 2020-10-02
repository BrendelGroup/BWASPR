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
#' @param nbrxtrms integer; number of extreme regions to evaluate
#' @param outfile If specified, then output is saved in the specified file name.
#' @param doplots Logical; if true, then show plots of distances
#'
#' @return List of two data frames, recording hsm rich and poor regions.
#'
#' @importFrom methylKit getData getCoverageStats getMethylationStats
#' @importFrom dplyr arrange
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
#'   det_mrpr(AmHE[[1]],"Am_HE_fr",ddset=c(1,5),nbrxtrms=100L,
#'            outfile="dst-Am_HE_fr.txt",doplots=TRUE)
#'
#' @export

det_mrpr <- function(mrobj,sampleL,ddset=c(1,5),nbrxtrms=100L,outfile="",doplots=TRUE) {
    message("... determining methylation site rich and poor regions for " ,sampleL," ...")
# ... set checkflag=1 to print out loads of detail to stdout ...:
    checkflag <- 0
    sampledata <- getData(mrobj)
    sampleID <- mrobj@sample.id
    if (outfile != "") {
        sink(outfile)
    }

    DFhsmr <- data.frame( Rtype = as.character(), SeqId = as.character(), From = as.integer(),
                          To = as.integer(), Rlgth = as.integer(), NbrSites = as.integer(),
                          Sdnsty = as.numeric(), Fsite = as.integer(), Tsite = as.integer(),
                          Dvalue = as.numeric()
                        )
    
    for ( d in ddset ) {
      cat( sprintf("\n") )
      plotme <- list()
      xlabel <- sprintf("%s (%2d-distances)",sampleL,d)
      
      cat( sprintf("Analysis of %2d-distances for sample \"%s\" (%d sites):\n\n",
           d,sampleL,nrow(sampledata)) )
      
      # gdvector = vector of d-distances over the entire data set:
      #
      gdvector <- matrix(NA,1,nrow(sampledata))
      gdbegidx <- matrix(NA,1,nrow(sampledata))
      gdendidx <- matrix(NA,1,nrow(sampledata))
      
      seqid <- sampledata$chr[1]
      gndstnc= 0
      nsites= 0
      
      for ( i in 1:nrow(sampledata) ) {
        if ( !identical(sampledata$chr[i],seqid) ) {
          if (checkflag) {
            cat( sprintf("Check1: %s\thas %4d sites and %4d %2d-distances\n",
                         seqid, nsites, nsites-d, d) )
          }
          seqid <- sampledata$chr[i]
          nsites= 1
          if ( i == nrow(sampledata)  &&  checkflag ) {
            cat( sprintf("Check2: %s\thas %4d sites and %4d %2d-distances\n",
                         sampledata$chr[i], 1, 0, d) )
          }
        } else if ( i == nrow(sampledata) ) {
          nsites <- nsites + 1
          if ( nsites >= d ) {
            dstnc = sampledata$start[i] - sampledata$start[i-d]
            if (checkflag) {
              cat( sprintf("\t%s\t%9d %9d %d-dstnc:\t%6d\n",
                   sampledata$chr[i], sampledata$start[i-d],
                   sampledata$start[i], d, dstnc) ) 
            }
            gndstnc <- gndstnc + 1
            gdvector[gndstnc] <- dstnc
            gdbegidx[gndstnc] <- i-d
            gdendidx[gndstnc] <- i
          }
          if (checkflag) {
            cat( sprintf("Check3: %s\thas %4d sites and %4d %2d-distances\n",
                 seqid, nsites, max(0,nsites-d), d) )
          }
        } else {
          nsites <- nsites + 1
          if ( nsites >  d ) {
            dstnc = sampledata$start[i] - sampledata$start[i-d]
            if (checkflag) {
              cat( sprintf("\t%s\t%9d %9d %d-dstnc:\t%6d\n",
                   sampledata$chr[i], sampledata$start[i-d],
                   sampledata$start[i], d, dstnc) ) 
            }
            gndstnc <- gndstnc + 1
            gdvector[gndstnc] <- dstnc
            gdbegidx[gndstnc] <- i-d
            gdendidx[gndstnc] <- i
          }
        }
      }
      cat( sprintf("  Total number of %2d-distances:\t%6d\n", d, gndstnc) )

      # ... cutting gd* to included only the values calculated:
      #
      gdvector <- gdvector[1:gndstnc]
      gdbegidx <- gdbegidx[1:gndstnc]
      gdendidx <- gdendidx[1:gndstnc]

      # ... finding the top and bottom values:
      #
      sgdvector <- sort(gdvector)
      if (length(sgdvector) >= nbrxtrms) {
        ldcutoff <- sgdvector[length(sgdvector)-nbrxtrms+1]
        hdcutoff <- sgdvector[nbrxtrms]
      }
      else { # ... this would be very odd, indeed, but to be safe:
        ldcutoff <- sgdvector[length(sgdvector)]
        hdcutoff <- sgdvector[1]
      }

      # ... getting distance statistics:
      #
      ddlist <- ddstats(gdvector,xlabel,"black")
      dstats <- ddlist$dstats
      plotme[[1]] <- ddlist$plot1
      plotme[[2]] <- ddlist$plot2
      cat( sprintf("  Median: %7.1f Mean: %7.1f Std: %7.1f\n",
           dstats["median"], dstats["mean"], dstats["std.dev"]) )
      cat( sprintf("  Quantiles:\n") )
      q <- quantile(gdvector,probs=seq(0,1,0.05))
      print( q )

      cat( sprintf("\n") )
      cat( sprintf("  Low  density region distance cutoff (lowest  100): %7d\n", ldcutoff) )
      cat( sprintf("  High density region distance cutoff (highest 100): %7d\n", hdcutoff) )
      cat( sprintf("\n") )

      dhsmr.name <- ""
      phsmr.name <- ""

      for ( i in 1:length(gdvector) ) {
        if ( d > 1  &&  gdvector[i] <= hdcutoff ) {
          if (dhsmr.name == "") {
            dhsmr.name <- as.character(sampledata$chr[gdbegidx[i]])
            dhsmr.begi <- gdbegidx[i]
            dhsmr.begpnt <- sampledata$start[gdbegidx[i]]
            dhsmr.endi <- gdendidx[i]
            dhsmr.endpnt <- sampledata$start[gdendidx[i]]
          } else if (as.character(sampledata$chr[gdbegidx[i]]) == dhsmr.name  &&
                     sampledata$start[gdbegidx[i]] <= dhsmr.endpnt) {
            dhsmr.endi <- gdendidx[i]
            dhsmr.endpnt <- sampledata$start[gdendidx[i]]
          } else {
            # ... overlapping regions are merged to indicate the extent of the
            #     methylation-rich region:
            #
            Rlgth = dhsmr.endpnt-dhsmr.begpnt+1
            NbrSites = dhsmr.endi-dhsmr.begi+1
            Sdnsty = 1000.*NbrSites/Rlgth
            hsmr.record <- data.frame( Rtype = "Rich", Sample = sampleID, SeqID = dhsmr.name,
              From = dhsmr.begpnt, To = dhsmr.endpnt, Rlgth = Rlgth, NbrSites = NbrSites,
              Sdnsty = Sdnsty, Fsite = dhsmr.begi, Tsite = dhsmr.endi, Dvalue = d )
            DFhsmr <- rbind(DFhsmr,hsmr.record)
            if (checkflag) {
              cat( sprintf("Rich \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\tnumber of sites: %4d, per-kb-density: %6.2f\n",
                sampleL, dhsmr.name, dhsmr.begpnt, dhsmr.begi, dhsmr.endpnt, dhsmr.endi, d,
                Rlgth, NbrSites, 1000.*NbrSites/Rlgth ) )
              for ( j in dhsmr.begi:dhsmr.endi ) {
                print( sampledata[j,] )
              }
            }
            dhsmr.name <- as.character(sampledata$chr[gdbegidx[i]])
            dhsmr.begi <- gdbegidx[i]
            dhsmr.begpnt <- sampledata$start[gdbegidx[i]]
            dhsmr.endi <- gdendidx[i]
            dhsmr.endpnt <- sampledata$start[gdendidx[i]]
          }
          if (checkflag) {
            cat( sprintf("Rich hsm region %s from %d to %d, %2d-distance:\t%6d\n",
              sampledata$chr[gdbegidx[i]], sampledata$start[gdbegidx[i]],
              sampledata$start[gdendidx[i]], d,
              sampledata$start[gdendidx[i]] - sampledata$start[gdbegidx[i]] ) )
            for ( j in gdbegidx[i]:gdendidx[i] ) {
              print( sampledata[j,] )
            }
          }
        } else if ( d == 1  &&  gdvector[i] >= ldcutoff ) {
          if (phsmr.name == "") {
            phsmr.name <- as.character(sampledata$chr[gdbegidx[i]])
            phsmr.begi <- gdbegidx[i]
            phsmr.begpnt <- sampledata$start[gdbegidx[i]]
            phsmr.endi <- gdendidx[i]
            phsmr.endpnt <- sampledata$start[gdendidx[i]]
          } else if (as.character(sampledata$chr[gdbegidx[i]]) == phsmr.name  &&
                     sampledata$start[gdbegidx[i]] <= phsmr.endpnt) {
            phsmr.endi <- gdendidx[i]
            phsmr.endpnt <- sampledata$start[gdendidx[i]]
          } else {
            Rlgth = phsmr.endpnt-phsmr.begpnt+1
            NbrSites = phsmr.endi-phsmr.begi+1
            Sdnsty = 1000.*NbrSites/Rlgth
            hsmr.record <- data.frame( Rtype = "Poor", Sample = sampleID, SeqID = phsmr.name,
              From = phsmr.begpnt, To = phsmr.endpnt, Rlgth = Rlgth, NbrSites = NbrSites,
              Sdnsty = Sdnsty, Fsite = phsmr.begi, Tsite = phsmr.endi, Dvalue = d )
            DFhsmr <- rbind(DFhsmr,hsmr.record)
            if (checkflag) {
              cat( sprintf("Poor \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\tnumber of sites: %4d, per-kb-density: %6.2f\n",
                sampleL, phsmr.name, phsmr.begpnt, phsmr.begi, phsmr.endpnt, phsmr.endi, d,
                Rlgth, NbrSites, 1000.*NbrSites/Rlgth ) )
              for ( j in phsmr.begi:phsmr.endi ) {
                print( sampledata[j,] )
              }
            }
            phsmr.name <- as.character(sampledata$chr[gdbegidx[i]])
            phsmr.begi <- gdbegidx[i]
            phsmr.begpnt <- sampledata$start[gdbegidx[i]]
            phsmr.endi <- gdendidx[i]
            phsmr.endpnt <- sampledata$start[gdendidx[i]]
          }
          if (checkflag) {
            cat( sprintf("Poor hsm region %s from %d to %d, %2d-distance:\t%6d\n",
              sampledata$chr[gdbegidx[i]], sampledata$start[gdbegidx[i]],
              sampledata$start[gdendidx[i]], d,
              sampledata$start[gdendidx[i]]-sampledata$start[gdbegidx[i]]+1) )
            for ( j in gdbegidx[i]:gdendidx[i] ) {
              print( sampledata[j,] )
            }
          }
        }
      }
      if ( d > 1 ) {
        Rlgth = dhsmr.endpnt-dhsmr.begpnt+1
        NbrSites = dhsmr.endi-dhsmr.begi+1
        Sdnsty = 1000.*NbrSites/Rlgth
        hsmr.record <- data.frame( Rtype = "Rich", Sample = sampleID, SeqID = dhsmr.name,
          From = dhsmr.begpnt, To = dhsmr.endpnt, Rlgth = Rlgth, NbrSites = NbrSites,
          Sdnsty = Sdnsty, Fsite = dhsmr.begi, Tsite = dhsmr.endi, Dvalue = d )
        DFhsmr <- rbind(DFhsmr,hsmr.record)
        if (checkflag) {
          cat( sprintf("Rich \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\tnumber of sites: %4d, per-kb-density: %6.2f\n",
            sampleL, dhsmr.name, dhsmr.begpnt, dhsmr.begi, dhsmr.endpnt, dhsmr.endi, d,
            Rlgth, NbrSites, 1000.*NbrSites/Rlgth ) )
          for ( j in dhsmr.begi:dhsmr.endi ) {
            print( sampledata[j,] )
          }
        }
      }
      if ( d == 1 ) {
        Rlgth = phsmr.endpnt-phsmr.begpnt+1
        NbrSites = phsmr.endi-phsmr.begi+1
        Sdnsty = 1000.*NbrSites/Rlgth
        hsmr.record <- data.frame( Rtype = "Poor", Sample = sampleID, SeqID = phsmr.name,
          From = phsmr.begpnt, To = phsmr.endpnt, Rlgth = Rlgth, NbrSites = NbrSites,
          Sdnsty = Sdnsty, Fsite = phsmr.begi, Tsite = phsmr.endi, Dvalue = d )
        DFhsmr <- rbind(DFhsmr,hsmr.record)
        if (checkflag) {
          cat( sprintf("Poor \"%s\"\thsm region %s\tfrom %8d (site %6d) to %8d (site %6d), (adjusted) %2d-distance:\t%6d\tnumber of sites: %4d, per-kb-density: %6.2f\n",
            sampleL, phsmr.name, phsmr.begpnt, phsmr.begi, phsmr.endpnt, phsmr.endi, d,
            Rlgth, NbrSites, 1000.*NbrSites/Rlgth ) )
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
        DFhsmrR <- subset(DFhsmr, Rtype == "Rich" & Sample == sampleID,
                          select = c(Rtype,Sample,SeqID,From,To,Rlgth,NbrSites,Sdnsty))
        DFhsmrR <- arrange(DFhsmrR, -Sdnsty)
        sampleL.DFhsmrR.df <- as.data.frame(DFhsmrR)
        write.table(sampleL.DFhsmrR.df, sprintf("mrr-%s.txt",sampleL),
                    quote=F, sep="\t", row.names=F, col.names=F)
        cat( sprintf("\n Top %3d (out of %5d) methylation-rich regions in \"%s\":\n",
                     nbrxtrms, dim(DFhsmrR)[1], sampleL) )
        print( head(DFhsmrR, n = nbrxtrms) )
      }
    
      if (d == 1) {
        DFhsmrP <- subset(DFhsmr, Rtype == "Poor" & Sample == sampleID,
                          select = c(Rtype,Sample,SeqID,From,To,Rlgth,NbrSites,Sdnsty))
        DFhsmrP <- arrange(DFhsmrP, +Sdnsty)
        sampleL.DFhsmrP.df <- as.data.frame(DFhsmrP)
        write.table(sampleL.DFhsmrP.df, sprintf("mpr-%s.txt",sampleL),
                    quote=F, sep="\t", row.names=F, col.names=F)
        cat( sprintf("\n Top %3d (out of %5d) methylation-poor regions in \"%s\":\n",
                     nbrxtrms, dim(DFhsmrP)[1], sampleL) )
        print( head(DFhsmrP, n = nbrxtrms) )
      }
      cat( sprintf("\n\n") )
    
      if (doplots) {
        pdffile <- sprintf("%dds-%s.pdf",d,sampleL)
        mytitle <- sprintf( "\n%s\n%2d-distance distributions", sampleL, d)
        suppressGraphics(ggsave(pdffile, do.call(marrangeGrob,
          list(grobs=plotme, nrow=1, ncol=2, top = mytitle)), width=7, height=7, units="in"))
      }
    }
    
    DFhsmrPRT <- data.frame( DFhsmr$SeqID, DFhsmr$From, DFhsmr$To, DFhsmr$Rtype,
                             DFhsmr$Sample, DFhsmr$Rlgth, DFhsmr$NbrSites, DFhsmr$Sdnsty)
    names(DFhsmrPRT) <- c("SeqID", "From", "To", "Rtype", "Sample", "Rlgth", "NbrSites", "Sdnsty")
    
    tabfile <- sprintf("mdr-%s.tab",sampleL)
    bedfile <- sprintf("mdr-%s.bed",sampleL)
    sink(tabfile)
    write.table( format(DFhsmrPRT,scientific=FALSE,digits=1L),
                 row.names = FALSE, quote = FALSE, sep="\t" )
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
  plot1 <-  ggplot(pme, aes(x=v)) + scale_x_continuous(xlabel,limits=c(0,600),
    breaks=seq(0,600,60)) + geom_histogram(aes(y=..density..),binwidth=20,closed="right",
    color="black",fill="white") + scale_y_continuous(limits=c(0,0.005)) +
    geom_density(alpha=.2,color=mycolor)
  pme <- subset(dstncdf,v<=60)
  plot2 <-  ggplot(pme, aes(x=v)) + scale_x_continuous(xlabel,limits=c(0,60),
    breaks=seq(0,60,4)) + geom_histogram(aes(y=..density..),binwidth=2,closed="right",
    color="black",fill="white") + scale_y_continuous(limits=c(0,0.05)) +
    geom_density(alpha=.2,color=mycolor)
  rlist <- list("dstats" = dstats, "plot1" = plot1, "plot2" = plot2)

  return(rlist)
}
