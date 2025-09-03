# MetaDICT: Microbiome Data Integration via Shared Dictionary Learning

## Contents

[Overview](#Overview)

[Installation](#Installation)

[Demo](#Demo)


## Overview:
MetaDICT is a method for the integration of microbiome datasets. This method is designed to remove batch effects and preserve biological variation while integrating heterogeneous datasets. MetaDICT consists of two stages: the first stage provides an initial estimation of batch effect via covariate balancing, while batch effects are defined as heterogeneous capturing efficiency in sequencing measurement; the second stage refines the estimation by shared dictionary learning. Compared with other existing methods, MetaDICT can better avoid overcorrection when unobserved confounding variables are present. The integrated data set can be applied to downstream analysis such as PCoA, taxa/sample community detection, and differential abundance test.

## Installation:

Our R package has been tested on R version 4.4.2. The package can be downloaded from GitHub: 

```
devtools::install_github("BoYuan07/MetaDICT", build_vignettes = TRUE)
```

Load the package:

```
library(MetaDICT)
```




## Demo:

A demo dataset is included in the package. The tutorial can be found in the vignette:

```
vignette("MetaDICT")
```

# Reference

Bo Yuan, Shulei Wang,
<b>Microbiome Data Integration via Shared Dictionary Learning</b>
(2025), Nature Communications, 16(1), 8147.
[<a href=https://www.nature.com/articles/s41467-025-63425-y>link</a>]
