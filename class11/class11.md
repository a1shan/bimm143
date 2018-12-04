Structural Bioinfo Class11
================

PDB database composition statistics
-----------------------------------

Download PDB CSV Statistiscs from <http://www.rcsb.org/stats/summary>

``` r
pdbstats <- read.csv("Data Export Summary.csv", row.names = 1)
pdbstats
```

    ##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## X-Ray                 122263          1960               6333    10 130566
    ## NMR                    10898          1263                253     8  12422
    ## Electron Microscopy     1822            31                657     0   2510
    ## Other                    244             4                  6    13    267
    ## Multi Method             119             5                  2     1    127

Edit table format

``` r
library(knitr)
kable(pdbstats)
```

|                     |  Proteins|  Nucleic.Acids|  Protein.NA.Complex|  Other|   Total|
|---------------------|---------:|--------------:|-------------------:|------:|-------:|
| X-Ray               |    122263|           1960|                6333|     10|  130566|
| NMR                 |     10898|           1263|                 253|      8|   12422|
| Electron Microscopy |      1822|             31|                 657|      0|    2510|
| Other               |       244|              4|                   6|     13|     267|
| Multi Method        |       119|              5|                   2|      1|     127|

Determine percentage of structures by X-ray & Electron Microscopy

``` r
#Determine TOTAL structures
nstruc <- sum(pdbstats$Total)

#Determine % structures by experimental method
EMpro <- pdbstats$Total / nstruc
per <- round(EMpro, 4)*100
```

There are 89.49% X-ray structures and 1.72% Electron Microscopy structures in the PDB database as of 2018-11-06

``` r
#add percent to table
nstats <- pdbstats
nstats$Percent <- per
kable(nstats)
```

|                     |  Proteins|  Nucleic.Acids|  Protein.NA.Complex|  Other|   Total|  Percent|
|---------------------|---------:|--------------:|-------------------:|------:|-------:|--------:|
| X-Ray               |    122263|           1960|                6333|     10|  130566|    89.49|
| NMR                 |     10898|           1263|                 253|      8|   12422|     8.51|
| Electron Microscopy |      1822|             31|                 657|      0|    2510|     1.72|
| Other               |       244|              4|                   6|     13|     267|     0.18|
| Multi Method        |       119|              5|                   2|      1|     127|     0.09|

Proportion of structures are protein?

``` r
proteinpro <- round(sum(nstats$Proteins) / nstruc *100, 2)
```

92.77% of structures that are protein as of 2018-11-06
