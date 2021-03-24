# popest_clstr

Codes to run the hierarchical clustering model to estimate number of groups and group-size. 
The code uses two dataset, the spatial location of the capture and the maximum number of individuals captured at the location. 

I used the WeightedCluster(http://mephisto.unige.ch/weightedcluster/), factoextra (https://rpkgs.datanovia.com/factoextra/index.html) and Nbclust (https://cran.r-project.org/web/packages/NbClust/index.html) package in R for the analysis.

```{r}
library(WeightedCluster)
library(factoextra)
library(NbClust)
```
