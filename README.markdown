# What is the DSPR?

The *Drosophila* Synthetic Population Resource (DSPR) consists of a new panel of over 1700 recombinant inbred lines (RILs) of *Drosophila melanogaster*, derived from two highly recombined synthetic populations, each created by intercrossing a different set of 8 inbred founder lines (with one founder line common to both populations). Complete genome sequence data for the founder lines are available, and in addition, there is a high resolution genetic map for each RIL. The DSPR has been developed as a community resource for high-resolution QTL mapping and is intended to be used widely by the *Drosophila* community.

# Advantages of the DSPR

## Benefits of RILs

Because the DSPR consists of recombinant inbred lines, genotyping need only be carried out once, phenotypes that cannot easily be scored on single individuals can be assayed (i.e., those of low heritability and/or with appreciable measurement error), and RILs facilitate the examination of genotype-by-sex, genotype-by-environment, and so on, interactions.

## Wide range of natural variation

The DSPR was created from a worldwide sample of 15 founder lines (South Africa, 2 lines; Asia, 4 lines; Australasia, 1 line; Europe, 2 lines; North America, 4 lines; South America, 2 lines). This panel therefore captures a large proportion of natural trait variation and provides a robust picture of trait genetic architecture.

## Extremely high mapping resolution

The DSPR is derived from highly recombined synthetic populations. This design allows for QTLs to be mapped to very narrow regions of the genome (1-2cM), encompassing just a handful of genes.

## Joint estimate of QTL effect and frequency

Each set of RILs is derived from eight inbred founders. This multi-founder composition permits an estimate of BOTH the effect and frequency of mapped QTL, allowing for an estimate of the contribution of each QTL to natural trait variation.

## Identification of causative sites

Researchers using the DSPR will be able to identify putative causative sites. For each mapped QTL, a QTL allele can be assigned to each founder. The whole genome re-sequencing data allows the small proportion of polymorphisms in phase with the QTL allelic configuration to be identified. This small set is highly likely to include the causative site.

# Installing the DSPRqtl Package

## Stable Release

Directions for installation of the stable release are in the [DSPRqtl tutorial (PDF)](http://wfitch.bio.uci.edu/DatFILES/DSPRqtl-intro.pdf).

## Development Version

You can install the most recent version directly from github using
`install_github()` from the [devtools
package](https://github.com/hadley/devtools).

```R
require("devtools")
install_github("DSPRqtl", "egking", branch = "develop")
```

Note that if you are using Windows, you will first need to install
Rtools. Start at <http://www.murdoch-sutherland.com/Rtools/> and
follow the links to CRAN.
