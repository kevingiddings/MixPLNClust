# MixPLNClust
R package for clustering multivariate discrete count data with incomplete values. Uses the EM algorithm to cluster according to mixtures of multivaraite Poisson log-normal distributitions and estimate latent parameters

## Installation

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
