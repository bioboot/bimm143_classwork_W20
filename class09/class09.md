Class 9: Machine learning pt1
================
Barry Grant
2/4/2020

## K-means clustering

The main k-means function in R is called `kmeans()`. Lets play with it
here.

``` r
# Generate some example data for clustering
tmp <- c(rnorm(30,-3), rnorm(30,3))
x <- cbind(x=tmp, y=rev(tmp))

plot(x)
```

![](class09_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Use the kmeans() function setting k to 2 and nstart=20

``` r
km <- kmeans(x, centers = 2, nstart = 20)
km
```

    ## K-means clustering with 2 clusters of sizes 30, 30
    ## 
    ## Cluster means:
    ##           x         y
    ## 1 -2.872733  3.317305
    ## 2  3.317305 -2.872733
    ## 
    ## Clustering vector:
    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    ## [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ## 
    ## Within cluster sum of squares by cluster:
    ## [1] 56.08682 56.08682
    ##  (between_SS / total_SS =  91.1 %)
    ## 
    ## Available components:
    ## 
    ## [1] "cluster"      "centers"      "totss"        "withinss"     "tot.withinss"
    ## [6] "betweenss"    "size"         "iter"         "ifault"

Inspect/print the results

Q. How many points are in each cluster? Q. What ‘component’ of your
result object details - cluster size? - cluster assignment/membership? -
cluster center?

``` r
# Cluster size
km$size
```

    ## [1] 30 30

``` r
km$cluster
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    ## [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
length(km$cluster)
```

    ## [1] 60

``` r
table(km$cluster)
```

    ## 
    ##  1  2 
    ## 30 30

Plot x colored by the kmeans cluster assignment and add cluster centers
as blue points

``` r
plot(x, col=km$cluster)
points(km$centers, col="blue", pch=15, cex=3)
```

![](class09_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Hierarchical clustering in R

The main Hierarchical clustering function in R is called `hclust()` An
important point here is that you have to calculate the distance matrix
deom your input data before calling `hclust()`.

For this we will use the `dist()` function first.

``` r
# We will use our x again from above...
d <- dist(x)
hc <- hclust(d)
hc
```

    ## 
    ## Call:
    ## hclust(d = d)
    ## 
    ## Cluster method   : complete 
    ## Distance         : euclidean 
    ## Number of objects: 60

Folks often view the results of Hierarchical clustering graphicaly. Lets
try passing this to the `plot()` function.

``` r
plot(hc)
abline(h=6, col="red", lty=2)
abline(h=4, col="blue")
```

![](class09_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

To get cluster membership vector I need to “cut” the tree at a certain
height to yield my seperate cluster
    branches.

``` r
cutree(hc, h=6)
```

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2
    ## [39] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

``` r
gp4 <- cutree(hc, k=6)
table(gp4)
```

    ## gp4
    ##  1  2  3  4  5  6 
    ##  8  7 15 15  7  8

Let’s try with some more real like data…

``` r
# Step 1. Generate some example data for clustering
x <- rbind(
  matrix(rnorm(100, mean=0, sd = 0.3), ncol = 2),   # c1
  matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2), # c2
  matrix(c(rnorm(50, mean = 1, sd = 0.3),           # c3
           rnorm(50, mean = 0, sd = 0.3)), ncol = 2))
colnames(x) <- c("x", "y")

# Step 2. Plot the data without clustering
plot(x)
```

![](class09_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
# Step 3. Generate colors for known clusters 
#         (just so we can compare to hclust results)
col <- as.factor( rep(c("c1","c2","c3"), each=50) ) 

plot(x, col=col)
```

![](class09_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

Q. Use the `dist()`, `hclust()`, `plot()` and `cutree()` functions to
return 2 and 3 clusters

Q. How does this compare to your known ‘col’ groups?

``` r
hc <- hclust( dist(x) )
plot(hc)
abline(h=1.7, col="red")
```

![](class09_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

To get cluster membership vector use `cutree()` and then use `table()`
to tabulate up how many members in each cluster we have.

``` r
grps <- cutree(hc, k=3)
table(grps)
```

    ## grps
    ##  1  2  3 
    ## 59 33 58

Make a plot with our cluster results

``` r
plot(x, col=grps)
```

![](class09_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## PCA of UK food data

``` r
x <- read.csv("UK_foods.csv", row.names = 1)
x
```

    ##                     England Wales Scotland N.Ireland
    ## Cheese                  105   103      103        66
    ## Carcass_meat            245   227      242       267
    ## Other_meat              685   803      750       586
    ## Fish                    147   160      122        93
    ## Fats_and_oils           193   235      184       209
    ## Sugars                  156   175      147       139
    ## Fresh_potatoes          720   874      566      1033
    ## Fresh_Veg               253   265      171       143
    ## Other_Veg               488   570      418       355
    ## Processed_potatoes      198   203      220       187
    ## Processed_Veg           360   365      337       334
    ## Fresh_fruit            1102  1137      957       674
    ## Cereals                1472  1582     1462      1494
    ## Beverages                57    73       53        47
    ## Soft_drinks            1374  1256     1572      1506
    ## Alcoholic_drinks        375   475      458       135
    ## Confectionery            54    64       62        41

Lets make some plots to explore our data a bit more

``` r
barplot(as.matrix(x), beside=F, col=rainbow(nrow(x)))
```

![](class09_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
barplot(as.matrix(x), beside=T, col=rainbow(nrow(x)))
```

![](class09_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

A “pairs” plot can help when we have small datasets like this one but
often we are dealing with data that is too large for these appraches :-(

``` r
pairs(x, col=rainbow(10), pch=16)
```

![](class09_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Principal Component Analysis (PCA) with the `prcomp()` function.

``` r
pca <- prcomp( t(x) )
#pca
```

What is in my result object `pca`? I can check the attributes…

``` r
attributes(pca)
```

    ## $names
    ## [1] "sdev"     "rotation" "center"   "scale"    "x"       
    ## 
    ## $class
    ## [1] "prcomp"

``` r
plot(pca$x[,1], pca$x[,2], xlab="PC1", ylab="PC2")
text(pca$x[,1], pca$x[,2], colnames(x), 
     col=c("black", "red","blue","darkgreen"))
```

![](class09_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
summary(pca)
```

    ## Importance of components:
    ##                             PC1      PC2      PC3       PC4
    ## Standard deviation     324.1502 212.7478 73.87622 4.189e-14
    ## Proportion of Variance   0.6744   0.2905  0.03503 0.000e+00
    ## Cumulative Proportion    0.6744   0.9650  1.00000 1.000e+00
