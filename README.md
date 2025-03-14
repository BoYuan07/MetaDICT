# MetaDICT: Microbiome Data Integration via Shared Dictionary Learning

## Contents

[Overview](#Overview)

[Installation](#Installation)

[Demo](#Demo)


## Overview:
MetaDICT is a method for the integration of microbiome datasets. This method is designed to remove batch effects and preserve biological variation while integrating heterogeneous datasets. MetaDICT consists of two stages: the first stage provides an initial estimation of batch effect via covariate balancing, while batch effects are defined as heterogeneous capturing efficiency in sequencing measurement; the second stage refines the estimation by shared dictionary learning. Compared with other existing methods, MetaDICT can better avoid overcorrection when unobserved confounding variables are present. The integrated data set can be applied to downstream analysis such as PCoA, taxa/sample community detection, and differential abundance test.

## Installation:

Our R package has been tested on R version 4.4.2. Users should install the following packages before installing MetaDICT:

```
install.packages("c('stats', 'RANN', 'igraph', 'vegan', 'edgeR', 'ecodist', 'ggplot2', 'viridis')")
```

The package can be downloaded from GitHub: 

```
devtools::install_github("BoYuan07/MetaDICT", build_vignettes = TRUE)
```

Load the package:

```
library(MetaDICT)
```

The package should take approximately 223 seconds to install with vignette.



## Demo:

A demo dataset is included in the package. The tutorial can be found in the vignette:

```
vignette("MetaDICT")
```
The estimated running time of MetaDICT on the demo dataset is approximately 65 seconds.

# Reference

Bo Yuan, Shulei Wang,
<b>Microbiome Data Integration via Shared Dictionary Learning</b>
(2024).
[<a href=https://www.biorxiv.org/content/10.1101/2024.10.04.616752v1>link</a>]
