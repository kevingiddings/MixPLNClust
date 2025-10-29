# MixPLNClust
R package for clustering multivariate discrete count data with incomplete values. Uses the EM algorithm to cluster according to mixtures of multivaraite Poisson log-normal distributitions and estimate latent parameters

## Installation

Requires EdgeR, a Biocunductor pacakge. This package can be installed with the following code:

```{r}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
```

```{r}
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("kevingiddings/MixPLNClust")
```

## Vignette

There is a vignette avaiable that discusses a simple use case for the package:

```{r}
vignette("MixPLNClust")
```

## Citations

Chen Y, Chen L, Lun ATL, Baldoni P, Smyth GK (2025). “edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets.” Nucleic Acids Research, 53(2), gkaf018. doi:10.1093/nar/gkaf018. 
