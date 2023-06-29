
# Create a new RStudio project

Open RStudio and create a new project, for more info see [Using-Projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects)

* File > New Project > New Directory > New Project (name the new directory, Ex. mRNA_Seq_Workshop).

Learn more about [renv](https://rstudio.github.io/renv/articles/renv.html)

## Install the needed R packages

Set some options and make sure the packages edgeR, gplots, RColorBrewer, topGO, KEGGREST, Rgraphviz and org.Mm.eg.db are installed (if not install it), and then load.

In the R console run the following commands *one at a time*:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!any(rownames(installed.packages()) == "edgeR")){
  BiocManager::install("edgeR")
}
library(edgeR)

if (!any(rownames(installed.packages()) == "topGO")){
  BiocManager::install("topGO")
}
library(topGO)

if (!any(rownames(installed.packages()) == "clusterProfiler")){
  BiocManager::install("clusterProfiler")
}
library(clusterProfiler)

if (!any(rownames(installed.packages()) == "dplyr")){
  BiocManager::install("dplyr")
}
library(dplyr)

if (!any(rownames(installed.packages()) == "Rgraphviz")){
  BiocManager::install("Rgraphviz")
}
library(Rgraphviz)

if (!any(rownames(installed.packages()) == "org.Mm.eg.db")){
  BiocManager::install("org.Mm.eg.db")
}
library(org.Mm.eg.db)

if (!any(rownames(installed.packages()) == "gplots")){
    BiocManager::install("gplots")
}
library(gplots)

if (!any(rownames(installed.packages()) == "RColorBrewer")){
    BiocManager::install("RColorBrewer")
}
library(RColorBrewer)

if (!any(rownames(installed.packages()) == "ggplot2")){
    BiocManager::install("ggplot2")
}
library(ggplot2)


if (!any(rownames(installed.packages()) == "devtools")){
  BiocManager::install("devtools")
}
library(devtools)

devtools::install_github("javadnoorb/pathview")
```

Note about pathview: As of June 2022, the version of pathview on Bioconductor is (presumably temporarily) broken due to KEGG's move from http to https. Therefore, the instructions above install a patched version of pathview from Github.  

## Download the template Markdown workshop document and open it

In the R console run the following command

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-June-RNA-Seq-Analysis/master/data_analysis/DE_Analysis_mm.Rmd", "DE_Analysis_mm.Rmd")
```

## Download the data file for the workshop document and preview/open it

This is the the counts file generated after running [Generating counts tables](https://ucdavis-bioinformatics-training.github.io/2023-June-RNA-Seq-Analysis/data_reduction/counts).

I've also uploaded to the github repo. In the R console run the following command.
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-June-RNA-Seq-Analysis/master/datasets/rnaseq_workshop_counts.txt", "rnaseq_workshop_counts.txt")
```

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-June-RNA-Seq-Analysis/master/datasets/ensembl_mm_106.tsv", "ensembl_mm_106.tsv")
```

#### For the salmon datasets

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2023-June-RNA-Seq-Analysis/master/datasets/rnaseq_salmon_workshop_counts.txt", "rnaseq_salmon_workshop_counts.txt")
```

