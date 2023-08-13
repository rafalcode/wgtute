library(WGCNA)
library(Cairo)

# load prev vars ...
load("FemaleLiver-01-dataInput.RData");
load("FemaleLiver-02-networkConstruction-auto.RData");

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
# sizeGrWindow(9,9)
CairoPNG("femliv_plotTOM.png", 900, 900)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

nSelect = 400

# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];

# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];

# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
CairoPNG("femliv_plotDiss.png", 900, 900)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))

# Plot the relationships among the eigengenes and the trait
# sizeGrWindow(5,7.5);
CairoPNG("femliv_eignw_rels.png", 500, 800)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
dev.off()

# Plot the dendrogram
# sizeGrWindow(6,6);
CairoPNG("femliv_eignw_dend.png", 600, 600)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
CairoPNG("femliv_eignw_hm.png", 600, 600)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()
