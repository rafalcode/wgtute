library(WGCNA)
library(Cairo)

#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv") # 3600 row (genes) 132 samples (columns)
# Read in the male liver data set
maleData = read.csv("LiverMale3600.csv") # 3600 rows - that's handy - and 143 samples.

# We work with two sets / SO, there's two sets. Fem and Male, right? Sooo obvious, but?
nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

# this data has its first 8 columns with descriptive, non expression data that we want out.
# as in:
# substanceBXH gene_symbol LocusLinkID ProteomeID cytogeneticLoc CHROMOSOME StartPosition EndPosition
multiExpr[[1]] = list(data = as.data.frame(t(femData[-c(1:8)])));
names(multiExpr[[1]]$data) = femData$substanceBXH;
rownames(multiExpr[[1]]$data) = names(femData)[-c(1:8)];
multiExpr[[2]] = list(data = as.data.frame(t(maleData[-c(1:8)])));
names(multiExpr[[2]]$data) = maleData$substanceBXH;
rownames(multiExpr[[2]]$data) = names(maleData)[-c(1:8)];
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr) # a WGCNA function, to check overal framework of sets (for consensus)

# Check that all genes and samples have sufficiently low numbers of missing values.
# standard exercise ... with real data you can expect some gruntwork here.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
if (!gsg$allOK) {
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

sampleTrees = list()
for (set in 1:nSets) {
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

CairoPNG("femmale_sampleClustering1.png", 1200, 1200)
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();

# As usual, cutting is carried out.
# he says:
# By inspection, there seems to be one outlier in the female data set, and no obvious outliers in the male set. We now remove the female outlier using a semi-automatic code that only requires a choice of a height cut. We first re-plot the two sample trees with the cut lines included

# Choose the "base" cut height for the female data set
baseHeight = 16
# Adjust the cut height for the male data set for the number of samples
cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1]);
# Re-plot the dendrograms including the cut lines
CairoPNG("femmale_sampleClustering2.png", 1200, 1200)
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets) {
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}
dev.off();

# cutting now performed. This is a basic R technique for getting rid of outliers.
for (set in 1:nSets) {
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}
collectGarbage()
# Check the size of the leftover data
exprSize = checkSets(multiExpr)

# all same as beff.
traitData = read.csv("ClinicalTraits.csv");
# dim(traitData)
# names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];

# See how big the traits are and what are the trait and sample names
# dim(allTraits)
# names(allTraits)
# allTraits$Mice

# Form a multi-set structure (ehem, two-set) that will hold the clinical traits.
# all this is, is:
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets) {
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$Mice);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();

# Define data set dimensions .. well render them into convenient variables I would say.
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

# saveRDS(multiExpr, "femmale_multiExpr.rds")
# saveRDS(Traits, "femmale_Traits.rds")
# saveRDS(nGenes, "femmale_nGenes.rds") # fairly pathetic, should go into source-able R "definitions" file.
# saveRDS(nSamples, "femmale_nSamples.rds")
# saveRDS(setLabels, "femmale_setLabels.rds")
# saveRDS(shortLabels, "femmmale_shortLabels.rds")
# saveRDS(exprSize, "femmale_exprData.rds")
# all that instead of RData file = "Consensus-dataInput.RData"
save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, file = "Consensus-dataInput.RData");
