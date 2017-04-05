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

A typical *BWASPR* workflow will read the specified \*.mcalls files and
generate various output tables and plots, labeled in various ways with
_species_\__study_\__sample_\__replicate_ labels.
