library(WGCNA)

# Basic settings: we work with two data sets
nSets = 2

# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Load the data saved in the first part
load("Consensus-dataInput.RData");
load("Consensus-NetworkConstruction-auto.RData");

# Create a variable weight that will hold just the body weight of mice in both sets
weight = vector(mode = "list", length = nSets);
for (set in 1:nSets) {
  weight[[set]] = list(data = as.data.frame(Traits[[set]]$data$weight_g));
  names(weight[[set]]$data) = "weight"
}

# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);

# We add the weight trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC, weight));

source("plotEigengene_rf.R")

pdf("femmale_EigengeneNetworks.pdf", width= 8, height = 10);
par(cex = 0.9)
plotEigengeneNetworks_rf(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      xLabelsAngle = 90)
                      # zlimPreservation = c(-1., 1), xLabelsAngle = 90)
dev.off();

# OOPs .. last hurdle ... getting
# Error in sum((1 - abs(dispd)) < zlimPreservation[1]) || ((1 - abs(dispd)) 
#   'length = 324' in coercion to 'logical(1)'
