# BWASPR
### R functions and scripts to process output from the BWASP workflow

This repository contains R functions and scripts we use to analyze the
\*.mcalls output files from the
[BWASP](https://github.com/brendelgroup/BWASP) workflow.
Required input consists of the \*.mcalls files (tab delimited data for the named
columns

```{bash}
SeqID.Pos SequenceID  Position  Strand Coverage  Prcnt_Meth  Prcnt_Unmeth
```
and two files specifying the data labels and \*.mcall file locations and certain
parameters, respectively. Let's look at the example files:

```
Am.dat
================================================================================
# Samples from Herb et al. (2012) Nature Neuroscience:
#
Am      HE      forager 0       CpGhsm  ./BWASPR/extdata/Amel-forager.CpGhsm.mcalls
Am      HE      forager 0       CpGscd  ./BWASPR/extdata/Amel-forager.CpGscd.mcalls
Am      HE      nurse   0       CpGhsm  ./BWASPR/extdata/Amel-nurse.CpGhsm.mcalls
Am      HE      nurse   0       CpGscd  ./BWASPR/extdata/Amel-nurse.CpGscd.mcalls
```

```
Am.par
================================================================================
SPECIESNAME     Apis mellifera
TOTALNBRPMSITES 20307353
ASSEMBLYVERSION Amel_4.5
GENOMESIZE      250270657
SPECIESGFF3DIR  ./BWASPR/extdata/AmGFF3DIR
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
_sample_ (here _forager_ and _nurse_"); replicate number (here _0_, indicating
single samples or, as in the case of this study, aggregates over replicates);
and file locations (here for the _CpGhsm_ and _CpGscd_ \*.mcalls files).
The second file species the genome assembly version, genome size (in base
pairs), total number of potential methylation sites (CpGs), and file names
for GFF3 annotation of various genomic features.

A typical *BWASPR* workflow will read the specified \*.mcalls files and
generate various output tables and plots, labeled in various ways with
_species_\__study_\__sample_\__replicate_ labels.
The [demo/Rscript.BWASPR](./demo/Rscript.BWASPR) file shows a template
workflow.
Initial customization is done at the top of the file and most from
inclusion of a configuration file such as
[demo/sample.conf](./demo/sample.conf).
The following table summarizes the successive workflow steps.
You may want to open the [demo/Rscript.BWASPR](./demo/Rscript.BWASPR) and
[demo/sample.conf](./demo/sample.conf) in separate windows as a reference
while viewing the table.

| RUNflag    | input                 | parameters                 | function                      | theme                                           | output files                                                           |
|------------|-----------------------|----------------------------|-------------------------------|-------------------------------------------------|------------------------------------------------------------------------|
| RUNcms     | studymk               | covlist                    | cmStats()                     | coverage and methylation statistics             | cms-\*.txt<br/>cms-\*.pdf                                              |
| RUNrnk     | studymk               | genome_ann$*region*        | rank_rbm()                    | ranked genes and promoters                      | rnk-sig-\*.txt<br/>rnk-sip-\*.txt                                      |
| RUNprr     | studymk               | ddset<br/>nr2d<br/>doplots | det_mprr()                    | methylation poor and rich regions               | dst-\*.txt<br/>1ds-\*.pdf<br/>5ds-\*.pdf<br/>mdr-\*.txt<br/>mdr-\*.bed |
| RUNrepcms  | replicate \*.mcalls   | covlist                    | cmStats()                     | coverage and methylation statistics             | cms-\*.txt<br/>cms-\*.pdf                                              |
| RUNrepcrl  | replicate \*.mcalls   | destrand                   | cmpSamples()                  | correlations between replicates                 | crl-\*.txt                                                             |
| RUNpwc     | studymk<br/>studymc   | -                          | cmpSites()                    | pairwise sample comparisons                     | pwc-\*.txt                                                             |
| RUNagscrl  | studymk               | destrand                   | cmpSamples()                  | correlations between aggregate samples          | crl-\*.txt                                                             |
| RUNacs     | studymk               | destrand                   | annotate_methylome()          | annotation of common sites                      | acs-\*.txt                                                             |
| RUNmmp     | studymk               | -                          | map_methylome()               | methylation to annotation maps                  | mmp-\*.txt                                                             |
| RUNdmsg    | sample \*.mcalls<br/> | highcoverage               | det_dmsg()                    | differentially methylated sites and genes       | dms-\*.txt<br/>dmg-\*.txt                                              |
| RUNdmgdtls | studyhc               | destrand                   | show_dmsg()                   | details for differentially methylated genes     | dmg-\*details\*.txt                                                    |
| RUNogl     | studyhc               | -                          | explore_dmsg()<br/>rank_dmg() | ranked lists of differentially methylated genes | ogl-\*.txt<br/>rnk-dmg-\*.txt                                          |
| RUNsave    | workflow output       | -                          | save.image()                  | save image of workflow output                   | \*.RData                                                               |
