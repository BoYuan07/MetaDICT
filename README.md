# MetaDICT: Microbiome Data Integration via Shared Dictionary Learning

## Overview:
MetaDICT is a method for the integration of microbiome datasets. This method is designed to remove batch effects and preserve biological variation while integrating heterogeneous datasets. MetaDICT consists of two stages: the first stage provides an initial estimation of batch effect via covariate balancing, while batch effects is defined as heterogeneous capturing efficiency in sequencing measurement; the second stage refines the estimation by shared dictionary learning. Compared with other existing methods, MetaDICT can better avoid overcorrection when unobserved confounding variables are present. The integrated data set can be applied to downstream analysis such as PCoA, taxa/sample community detection and differential abundance test.

## Installation:
The package can be downloaded from github: 

```
devtools::install_github("BoYuan07/MetaDICT, build_vignettes = TRUE")
```

Load the package:

```
library(MetaDICT)
```

# Reference

Bo Yuan, Shulei Wang,
<b>Microbiome Data Integration via Shared Dictionary Learning</b>
(2024).
[<a href=https://www.biorxiv.org/content/10.1101/2024.10.04.616752v1>link</a>]
