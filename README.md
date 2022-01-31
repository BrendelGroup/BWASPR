# BWASPR
### R functions and scripts to process output from the BWASP workflow

This repository contains R functions and scripts we use to analyze the
\*.mcalls output files from the [BWASP](https://github.com/brendelgroup/BWASP) workflow.


## Installation
Please find detailed installation instructions and options in the
[INSTALL](./INSTALL.md) document.


## Reference

__Claire Morandin and Volker P. Brendel (2021)__
 _Tools and applications for integrative analysis of DNA methylation in social insects._
Molecular Ecology Resources, 00, 1-19. [https://doi.org/10.1111/1755-0998.13566](https://doi.org/10.1111/1755-0998.13566).

Original pre-print: [at BioRxiv](https://www.biorxiv.org/content/10.1101/2021.08.19.457008v3).


## Contact

Please direct all comments and suggestions to
[Volker Brendel](<mailto:vbrendel@indiana.edu>)
at [Indiana University](http://brendelgroup.org/).


## What BWASPR does
Required input to the _BWASPR_ workflow consists of the \*.mcalls files (tab
delimited data for the named columns)

```{bash}
SeqID.Pos SequenceID  Position  Strand Coverage  Prcnt_Meth  Prcnt_Unmeth
```
and two files specifying the data labels and \*.mcalls file locations and
certain parameters, respectively.
Let's look at the example files in [inst/extdata](./inst/extdata):

```
AmHE.dat
================================================================================
# Samples from Herb et al. (2012) Nature Neuroscience:
#
Am      HE      forager 0       CpGhsm  ../inst/extdata/Amel-forager.CpGhsm.mcalls
Am      HE      forager 0       CpGscd  ../inst/extdata/Amel-forager.CpGscd.mcalls
Am      HE      nurse   0       CpGhsm  ../inst/extdata/Amel-nurse.CpGhsm.mcalls
Am      HE      nurse   0       CpGscd  ../inst/extdata/Amel-nurse.CpGscd.mcalls
```

```
AmHE.par
================================================================================
SPECIESNAME     Apis mellifera
ASSEMBLYVERSION Amel_4.5
GENOMESIZE      250270657
TOTALNBRPMSITES 20307353
SPECIESGFF3DIR  ../inst/extdata/AmGFF3DIR
GENELISTGFF3    Amel.gene.gff3
EXONLISTGFF3    Amel.exon.gff3
PCGEXNLISTGFF3  Amel.pcg-exon.gff3
PROMOTRLISTGFF3 Amel.promoter.gff3
CDSLISTGFF3     Amel.pcg-CDS.gff3
UTRFLAGSET      1
5UTRLISTGFF3    Amel.pcg-5pUTR.gff3
3UTRLISTGFF3    Amel.pcg-3pUTR.gff3
```

The first file has columns for _species_ (here _Am_); _study_ (here _HE_);
_sample_ (here _forager_ and _nurse_"); replicate number (here 0, indicating
single samples or, as in the case of this study, aggregates over replicates);
and file locations (here for the _CpGhsm_ and _CpGscd_ \*.mcalls files);
note that the file locations in this example are relative links, assuming you
will run the example discussed in the [demo](./demo) directory.
The second file specifies the species name, genome assembly version, genome
size (in base pairs), total number of potential methylation sites (CpGs), and
file names for GFF3 annotation of various genomic features (_UTRFLAGSET_ is set
to 1 to use UTR annotation in the GFF3 file).

A typical *BWASPR* workflow will read the specified \*.mcalls files and
generate various output tables and plots, labeled in various ways with
_species_\_ _study_\_ _sample_\_ _replicate_ labels.
The [demo/Rscript.BWASPR](./demo/Rscript.BWASPR) file shows a template
workflow.
Initial customization is done at the top of the file and mostly from
inclusion of a configuration file such as
[demo/sample.conf](./demo/sample.conf).
The following table summarizes the successive workflow steps.
You may want to open the [demo/Rscript.BWASPR](./demo/Rscript.BWASPR) and
[demo/sample.conf](./demo/sample.conf) in separate windows as a reference
while viewing the table.
Details on running the workflow with the demo data are given in
[demo/README](./demo/README.md).


## RUNflag to expected output correspondence

| RUNflag    | input                 | (select) parameters        | function                      | theme                                           | output files                                                                                                         |
|------------|-----------------------|----------------------------|-------------------------------|-------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|
| RUNcms     | studymk               | covlist, locount, hicount  | cmStats()                     | sample coverage and methylation statistics      | cms-\*.txt<br/>cms-\*.pdf                                                                                            |
| RUNpwc     | studymk<br/>studymc   | -                          | cmpSites()                    | pairwise sample comparisons                     | pwc-\*.vs.\*.txt                                                                                                     |
| RUNcrl     | studymk               | destrand                   | cmpSamples()                  | correlations between aggregate samples          | crl-\*.txt<br/>crl-\*.pdf                                                                                            |
|            |                       |                            |                               |                                                 |                                                                                                                      |
| RUNrepcms  | replicate \*.mcalls   | repcovlist,<br/>replocount, rephicount | cmStats()                     | replicate coverage and methylation statistics   | repcms-\*.txt<br/>repcms-\*.pdf                                                                                      |
| RUNrepcrl  | replicate \*.mcalls   | destrand                 | cmpSamples()                    | correlations between replicates                 | repcrl-\*.txt<br/>repcrl-\*.pdf                                                                                        |
|            |                       |                            |                               |                                                 |                                                                                                                      |
| RUNmmp     | studymk               | -                          | map_methylome()               | methylation to annotation maps                  | mmp-\*.txt                                                                                                           |
| RUNacs     | studymk               | destrand                   | annotate_methylome()          | annotation of common sites                      | acs-\*.txt                                                                                                           |
| RUNrnk     | studymk               | genome_ann$*region*        | rank_rbm()                    | ranked genes and promoters                      | ranked-\*.txt<br/>sites-in-\*.txt<br/>rnk-sig-\*.pdf<br/>sip-\*.txt<br/>rnk-sip-\*.txt<br/>rnk-sip-\*.pdf                |
| RUNmrpr    | studymk               | ddset<br/>nr2d<br/>doplots | det_mrpr()                    | methylation-rich and -poor regions              | dst-\*.txt<br/>\*ds-\*.pdf<br/>mdr-\*.tab<br/>mdr-\*.bed<br/>mpr-\*.txt<br/>mrr-\*.txt<br/>rmp-\*.txt<br/>gwr-\*.txt |
|            |                       |                            |                               |                                                 |                                                                                                                      |
| RUNdmt     | studymc               | wsize, stepsize            | det_dmt()                     | differentially methylated tiles and genes       | dmt-\*.txt<br/>dmg-\*.txt                                                                                            |
| RUNdmsg    | sample \*.mcalls<br/> | highcoverage<br/>destrand  | det_dmsg()                    | differentially methylated sites and genes       | dms-\*.txt<br/>dmg-\*.txt                                                                                            |
| RUNdmgdtls | studyhc               | destrand                   | show_dmsg()                   | details for differentially methylated genes     | dmg-\*.vs.\*_details.txt<br/>dmg-\*.vs.*_heatmaps.pdf                                                              |
| RUNogl     | studyhc               | -                          | explore_dmsg()<br/>rank_dmg() | ranked lists of differentially methylated genes | ogl-\*.txt<br/>rnk-dmg-\*.vs.\*.txt<br/>rnk-dmg-\*.vs.\*.pdf<br/>wrt-*.txt                                           |
|            |                       |                            |                               |                                                 |                                                                                                                      |
| RUNsave    | workflow output       | -                          | save.image()                  | save image of workflow output                   | \*.RData                                                                                                             |
