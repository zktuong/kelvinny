# kelvinny
Kelvin's wrapper scripts for R ggplot2 plotting functions

## Installation instructions
You can install the package via ```devtools::install_github()``` function in R
```R
library(devtools)
devtools::install_github('zktuong/kelvinny')
```
## Usage instructions
```R
library(kelvinny)
```
The package contains a couple of wrapper functions for plotting in R, mostly revolving around the use of ggplot2. I will update the package as i start writing more. use ?functionname to find out more options that each function can take.

### plotHeat
A shortcut to plotting heatmaps, provided the table you want to plot can be uploaded into R.
```R
data <- read.delim("table.txt") # a typical way to load a tab-delimited text file
# or in this specific example, i'm using a data table that is already preloaded in R
data(iris)
head(iris)
data <- iris[,1:4] # i'm removing the last column because it's not numeric
plotHeat(data)
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

### simple.tSNE
A simple function to overlay colours onto tSNE plots. It is currently written to accept a gene name and plot which cells express at least 1 transcript and separate the colours depending on the sample type.
```R
simple.tSNE('Epcam')
```
![basic simple.tSNE plot](exampleImages/simpleEpcam.png)
### ascendtSNE
An altered version of simple.tSNE to work with data from ascend. 
```R
ascendtSNE.gene('Epcam')
```
![basic ascendtSNE plot](exampleImages/ascendEpcam.png)
```R
ascendtSNE.gene('Epcam', heat=TRUE)
```
![ascendtSNE plot as a heatmap](exampleImages/ascendEpcamHeat.png)

```R
ascendtSNE.info('batch') # only accepts 'batch' or 'cluster' at the moment
```
![ascendtSNE plot with batch information](exampleImages/batch.png)

### colorMyTrajectory
Overlays what the simple.tSNE plot does in terms of coloring, onto the Monocle pseudotime trajectory generated with their CellDataSet (called HSMM in the tutorial).
```R
colorMyTrajectory('Epcam') 
```
![colorTrajectory](exampleImages/colorTrajectory.png)
