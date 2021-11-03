library(WGCNA)
library(reshape2)
library(stringr)

options(stringsAsFactors = FALSE)
type = "unsigned"
corType = "bicor"
corFnc = ifelse(corType=="bicor", cor, bicor)

#setwd("/ssd/KUE_nitab/")
exprMat <- "counts.matrix"
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                       quote="", comment="", check.names=F)
dim(dataExpr)
dataExpr <- as.data.frame(dataExpr) 
dataExpr <- dataExpr[rowSums(dataExpr==0)<10, ]

dim(dataExpr)
head(dataExpr)[,1:8]
dataExpr <- as.data.frame(t(dataExpr))

gsg = goodSamplesGenes(dataExpr, verbose = 3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)

dim(dataExpr)

head(dataExpr)[,1:8]

trait <- "traits.txt"
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
head(traitData)[,1:7]

#no trait data
#sampleTree = hclust(dist(dataExpr), method = "average")
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

sampleTree2 = hclust(dist(dataExpr), method = "average")
traitColors = numbers2colors(traitData, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitData), 
                    main = "Sample dendrogram and trait heatmap")

#=====================
#=====================

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")

power = sft$powerEstimate
power

net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, corType = corType, 
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)

save(net,file = "net10.Rdata")



load("net10.Rdata")

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

table(net$colors)
table(moduleColors)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

dynamicColors <- labels2colors(net$unmergedColors)

plotDendroAndColors(net$dendrograms[[1]], cbind(dynamicColors,moduleColors),
                    c("Dynamic Tree Cut", "Module colors"),
                    dendroLabels = FALSE, hang = 0.5,
                    addGuide = TRUE, guideHang = 0.05)

MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)


sizeGrWindow(5,7.5);
par(cex = 0.9)

plotEigengeneNetworks(MEs_col, "", marDendro = c(0,3,1,4), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOM
plotTOM = dissTOM^7; 
diag(plotTOM) = NA;

#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
nSelect = 100
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#=====================
robustY = ifelse(corType=="pearson",T,F)
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
sizeGrWindow(9,9)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)


labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.3, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               #textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))

if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}
write.table(geneTraitCor,file="geneTraitCor.txt",sep="\t")
write.table(geneTraitP,file="geneTraitP.txt",sep="\t")


## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.


module = "brown"
pheno = "weight_g"
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

TOM = TOMsimilarityFromExpr(dataExpr, power = power,corType = "bicor",networkType = "unsigned",nThreads=100)


probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edgestxt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0.02,
                               nodeNames = probes, nodeAttr = moduleColors)

# Select module probes
probes = colnames(dataExpr) 
inModule = (moduleColors==module)
modProbes = probes[inModule]

modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)


cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "edges.csv",
                               nodeFile = "nodes.csv",
                               weighted = TRUE, threshold = 0.11,
                               nodeNames = modProbes, nodeAttr = moduleColors[inModule])


#Hub gene screening
#Hub genes (Membership)
datKME=signedKME(dataExpr, MEs_col)
write.table(datKME, "kME_MM_test.txt",sep="\t")

HubGenes <- chooseTopHubInEachModule(datExpr,moduleColors)
write.table (HubGenes,file = "HubGenes_of_each_module.xls",quote=F,sep='\t')



