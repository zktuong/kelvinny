[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/gh/zktuong/kelvinny/branch/master/graph/badge.svg)](https://codecov.io/gh/zktuong/kelvinny)
[![Build Status](https://travis-ci.com/zktuong/kelvinny.svg?branch=master)](https://travis-ci.com/zktuong/kelvinny)
# kelvinny
Kelvin's wrapper scripts for R plotting functions and other stuff.
v2.0 - updated the functionality of the logistic regression functions, including a bootstrap wrapper.


## Installation instructions
You can install the package via ```devtools::install_github()``` function in R
```R
library(devtools)
devtools::install_github('zktuong/kelvinny', dependencies = TRUE)

# one function requires SummarizedExperiment from bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
```
## Usage instructions
```R
library(kelvinny)
```
The package contains a couple of wrapper functions for plotting in R, mostly revolving around the use of ggplot2, pheatmap and viridis etc. I will update the package as i start writing more. use ?functionname to find out more options that each function can take.
The plotHeat function has a small tutorial on how to generate some basic heatmaps. 

### plotHeat
A shortcut to plotting heatmaps, provided the table you want to plot can be uploaded into R.
```R
data <- read.delim("table.txt") # a typical way to load a tab-delimited text file
# or in this specific example, i'm using a data table that is already preloaded in R
data(iris)
head(iris)
data <- iris[,1:4] # i'm removing the last column because it's not numeric
plotHeat(data, col = "RdWhBlu")
```
![heatmap](exampleImages/heat.png)
```R
plotHeat(data, col="viridis")
```
![heatmap](exampleImages/heat2.png)
```R
plotHeat(data, col=c("green","black","red"))
```
![heatmap](exampleImages/heat3.png)

plotHeat will take additional arguments from pheatmap. So, using example for pheatmap: 
```R
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")
# Generate annotations for rows and columns
annotation_col = data.frame(
  CellType = factor(rep(c("CT1", "CT2"), 5)),
  Time = 1:5
)
rownames(annotation_col) = paste("Test", 1:10, sep = "")
annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")
# plot
plotHeat(test, annotation_col=annotation_col, annotation_row=annotation_row)
```
![heatmap](exampleImages/heat4.png)
```R
plotHeat(test, annotation_col=annotation_col, annotation_row=annotation_row, col="viridis")
```
![heatmap](exampleImages/heat5.png)

### gg_color_hue
Generates the standard colors used by ggplot.
```R
gg_color_hue(20)
#  [1] "#F8766D" "#EA8331" "#D89000" "#C09B00" "#A3A500" "#7CAE00" "#39B600" "#00BB4E" "#00BF7D" "#00C1A3" "#00BFC4" "#00BAE0" "#00B0F6"
# [14] "#35A2FF" "#9590FF" "#C77CFF" "#E76BF3" "#FA62DB" "#FF62BC" "#FF6A98"
```

### dirCreate/createDir
same in function as {base} dir.create, but always recursive
```R
dirCreate("/path/to/path/to/file")
createDir("/path/to/path/to/file/2")
```

### pbcopy/pbpaste
pbcopy lets you copy any object from R to paste outside as a dataframe/vector as you wish.
```R
pbcopy(data)
```
pbpaste does the reverse: converts what you copy outside and paste as a vector in R.
```R
pasted_data <- pbpaste()
```

### RFclassifier/RFpredictor
Uses RandomForest algorithm to classify data, for example seurat single-cell data
```R
# library(kelvinny)
## running RF classifier neat
classifier.class <- RFclassifier(train.seurat, training.classes = train.seurat@ident)
prediction.class <- RFpredictor(classifier.class, as.matrix(test.seurat@data))
RF_class <- prediction.class$prediction
names(RF_class) <- colnames(test.seurat@data) 
test.seurat <- AddMetaData(test.seurat, metadata = RF_class, "RF_class")

## running RF classifier but extract the probability
classifier <- RFclassifier(train.seurat, training.classes = train.seurat@ident, importance = "impurity", probability = TRUE)
prediction <- RFpredictor(classifier, as.matrix(test.seurat@data))

## plot this as a heatmap:
prediction_mat <- prediction$predictions
test.seurat.metadata <- cbind(test.seurat@meta.data, prediction$predictions)
test.seurat.metadata <- test.seurat.metadata[order(test.seurat.metadata$combined.clusters), ]

library(viridis)
plotmat <- t(test.seurat.metadata[,9:ncol(test.seurat.metadata)])
colAnno <- test.seurat.metadata[,c(5,7,8)]

library(paletteer)
palette <- paletteer_d(package = 'rcartocolor', palette = 'Vivid', 12)

RF_class <- gg_color_hue(15)
names(RF_class) <- levels(test.seurat.metadata$RF_class)[-4]
combined.clusters <- palette
names(combined.clusters) <- levels(test.seurat.metadata$combined.clusters)
protocol <- viridis::viridis(2)
names(protocol) <- unique(test.seurat.metadata$protocol)

anno_color <- list(RF_class = RF_class, combined.clusters = combined.clusters, protocol = protocol)

plotHeat(plotmat, color = inferno(50), scale = "none", annotation_col = colAnno, annotation_colors = anno_color, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE)
```
![heatmap](exampleImages/RFheat.png)

### train_model_glmnet/test_model_glmnet
Uses glmnet algorithm to predict data.
```R
model <- train_model_glmnet(data, variable_colname = "disease", alpha = 0.5, cutOff = 0.5, nfolds = ncol(data)) # LOOCV Elastic net regularization
pred <- test_model_glmnet(model = model, new_data = newdat, type = "link")
```

### trainScSimilarity/predScSimilarity
A cleaner version of above. Uses glmnet algorithm to predict similarity of cells to reference/training data.
Can take expression matrix as well as Seurat or SummarizedExperiment objects.
```R
### SummarizedExperiment object
model <- trainScSimilarity(train.sce, colData(train.sce)$CellType, test.sce, nfolds = dim(train.sce)[2])
pred <- predScSimilarity(model, test.sce)
### Seurat v2 object
model <- trainScSimilarity(train.seurat, seurat@ident, test.sce, nfolds = dim(train.seurat@data)[2])
pred <- predScSimilarity(model, test.seurat)

# bootstrap
# using Seurat V3 object as an example
pred <- similarity_bootstrap(train.seurat, Idents(train.seurat), test.seurat, nboots = 50, simplify = TRUE, verbose = FALSE) # simplify will toggle whether to return the mean and SD of the prediction, or all as nested lists, and verbose will toggle hide/unhide of the messages
```

### codon
Converts codon nucleotide sequences to amino acids
```R
codon("ttt")
# [1] "Phe"
codon("ttt", "ttg")
# [1] "Phe>Leu"
```

### parse_gmt
reads a .gmx file and automatically convert to a table, like a .gmt file
```R
parse_gmt("file.gmt")
```

### mtx_to_h5/mtx_to_h5totxt
writes a matrix to .h5 format
A complementary python module is also available to convert the .h5 file to txt
```R
mtx_to_h5(counts, "counts.h5")
```
and then in bash
```bash
pip install kelvinnypy 
# or 
conda install -c kt16 kelvinnypy

# after installation
h5totxt counts.h5
```
or if you want to convert to .h5 and .txt in a single function in R
```R
mtx_to_h5totxt(counts, "counts.h5")
```
### As a comparison
```R
norm_counts <- matrix(rexp(2e8, rate=.1), ncol=20000, dimnames = list(paste0("gene", 1:10000), paste0("cell", 1:20000)))
dim(norm_counts)
# [1] 10000 20000
```
### write to .txt file
```R
start <- Sys.time()
write.table(norm_counts, "./norm_counts.txt", quote = FALSE, sep = "\t", row.names = TRUE)
end <- Sys.time()
end - start
Time difference of 4.033704 mins
```
### write to .txt.gz file
```R
start <- Sys.time()
write.table(norm_counts, gzfile("./norm_counts.txt.gz"), quote = FALSE, sep = "\t", row.names = TRUE)
end <- Sys.time()
end - start
Time difference of 9.407493 mins
```
### write to .h5 file
```R
# library(kelvinny)
# write to .h5 
start <- Sys.time()
mtx_to_h5(norm_counts, "norm_counts.h5")
end <- Sys.time()
end - start
Time difference of 2.933439 mins
```
this file is compressed too. Can be used for transferring rather than transfer the .txt file, or gzipped .txt.gz file. 
### convert .h5 file
```bash
time h5totxt norm_counts.h5
real    2m46.901s
user    2m40.488s
sys     0m6.380s
```

### write to .h5 and convert to .txt 
```R
start <- Sys.time()
mtx_to_h5totxt(norm_counts, "norm_counts.h5")
end <- Sys.time()
end - start
Time difference of 5.651026 mins
```

I think ... you save a bit of time from the gzipping, transfer, and unzipping of the large .txt if you need to move them. But you could always move the .RDS to begin with... ¯\_(ツ)_/¯
