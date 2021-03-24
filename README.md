# popest_clstr

Codes to run the hierarchical clustering model to estimate number of groups and group-size. 
The code uses two dataset, the spatial location of the capture and the maximum number of individuals captured at the location. 

I used the WeightedCluster (http://mephisto.unige.ch/weightedcluster/), factoextra (https://rpkgs.datanovia.com/factoextra/index.html) and Nbclust (https://cran.r-project.org/web/packages/NbClust/index.html) package in R for the analysis.

```{r}
library(WeightedCluster)
library(factoextra)
library(NbClust)
```
The dataset requird need to have the following format,
```{r}
Latitude Longtitude Number of indiviudals 
                        photo-captured
20.27764  79.38614  4
20.29808  79.37403  3
20.252    79.38497  1
20.32906  79.29086  2
20.33533  79.34514  2
20.36006  79.31292  3
```
Running clustering with the data
```{r}
x<- cbind(Longitude, Latitude)
dst <-dist(x)
clstr <- hclust(dst, method = "complete", members = Number of individuals)
```

