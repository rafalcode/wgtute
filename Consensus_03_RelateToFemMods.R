library(WGCNA)

# Load the data saved in the first part
load("Consensus-dataInput.RData");
load("Consensus-NetworkConstruction-auto.RData");
load("FemaleLiver-02-networkConstruction-auto.RData")

# Rename variables to avoid conflicts
femaleLabels = moduleLabels;
femaleColors = moduleColors;
femaleTree = geneTree;
femaleMEs = orderMEs(MEs, greyName = "ME0");

load("Consensus-NetworkConstruction-auto.RData")

# Isolate the module labels in the order they appear in ordered module eigengenes
# Beware in manual version ME's already have colours, but in auto, no, they are labels ME1, mE2, etc.
femModuleLabels = substring(names(femaleMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
femModules = labels2colors(as.numeric(femModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))

# Numbers of female and consensus modules
nFemMods = length(femModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nFemMods, ncol = nConsMods);
CountTbl = matrix(0, nrow = nFemMods, ncol = nConsMods);

# Execute all pairwaise comparisons
for (fmod in 1:nFemMods)
  for (cmod in 1:nConsMods) {
    femMembers = (femaleColors == femModules[fmod]);
    consMembers = (moduleColors == consModules[cmod]);
    pTable[fmod, cmod] = -log10(fisher.test(femMembers, consMembers, alternative = "greater")$p.value);
    CountTbl[fmod, cmod] = sum(femaleColors == femModules[fmod] & moduleColors == consModules[cmod])
  }

# Truncate p values smaller than 10^{-50} to 10^{-50} 
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
pTable[pTable>50 ] = 50 ;
# Marginal counts (really module sizes)
femModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# Actual plotting
# sizeGrWindow(10,7 );
pdf(file = "ConsensusVsFemaleModules.pdf", wi = 10, he = 7);
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", femModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("Fem ", femModules, ": ", femModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of Female set-specific and Female-Male consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
dev.off();
