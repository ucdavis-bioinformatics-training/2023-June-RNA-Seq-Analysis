library(topGO)

# Example on how to do GO enrichment when you don't have an org.Xx.eg.db database

# read in top table
infile <- "WT.C_v_WT.NC.txt"
tmp <- read.delim(infile)

# Create gene list
geneList <- tmp$P.Value
names(geneList) <- tmp$Gene.stable.ID.version

# read in list of ensembl gene ids and GO terms from ensembl (GO terms are under "external")
# Doesn't necessarily need to come from ensembl
anno <- read.delim("gene_to_GO.txt")

# Create a named list where each element corresponds to a GO term and contains all of the 
# genes corresponding to a GO term
GO2genes <- tapply(anno$Gene.stable.ID.version, anno$GO.term.accession, function(x)x)

# Create a topGOData object using our GO2gene object
# Using CC ontology just to make example run faster
GOdata <- new("topGOdata",
              ontology = "CC",
              allGenes = geneList,
              geneSelectionFun = function(x)x,
              annot = annFUN.GO2genes, GO2genes = GO2genes)

# Proceed as with previous example
resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)


