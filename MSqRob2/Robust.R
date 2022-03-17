#	R version > 4.1.0 https://cran.r-project.org/

if(!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
BiocManager::install("msqrob2")

if(!requireNamespace("remotes", quietly = TRUE)) {BiocManager::install("remotes")}
BiocManager::install("statomics/msqrob2gui")

library(msqrob2)

# data(pe)
# pe <- aggregateFeatures(pe,i="peptide",fcol="Proteins",name="protein")
# pe <- msqrob(pe,i="protein",formula=~condition,modelColumnName="rlm")
# getCoef(rowData(pe[["protein"]])$rlm[[1]])

library(msqrob2gui)
launchMsqrob2App()

library(tidyverse)
library(limma)
library(QFeatures)
library(msqrob2)

setwd("/Users/adams/Documents/PhD/Projects/Proteomics\ Jhana/APP_NOS")

annotationFile <- "APP-NOS_experimental_annotation.xlsx"
featuresFile <- "peptides.txt"

groupBy <- "Proteins"
logTrans <- TRUE
removeRazor <- TRUE
filterColumns <- c("Reverse", "Potential.contaminant")
minObsFeat <- 2L
normMethod <- "center.median"

sumMethod <- "robust"

# form <- "~treatment*time"
# doRidge <- FALSE

# contrast <- "treatmentLN=0"
# sigLevel <- 0.05

maxPlot <- 10L

#	Data import
ecols <- grep("Intensity\\.", names(read.delim(featuresFile)))

pe <- readQFeatures(
    table = featuresFile, fnames = 1, ecol = ecols,
    name = "rawFeatures", sep = "\t")

currentAssay <- "rawFeatures"

exp_annotation <- as.data.frame(unclass(openxlsx::read.xlsx(annotationFile)))

runName <- which(
  vapply(exp_annotation, function(x) return(
    identical(
      sort(
        as.character(
          colnames(pe[[1]]))
        ),
      sort(
        as.character(x))
      )
    ), FUN.VALUE = TRUE))
if (length(runName)>0) runName <- runName[1]
rownames(exp_annotation) <- as.character(
  exp_annotation[,runName])

exp_annotation <- exp_annotation[colnames(pe[[1]]),]

for (j in colnames(exp_annotation))
  colData(pe)[[j]] <- exp_annotation[,j]

#	Preprocessing
rowData(pe[[currentAssay]])$nNonZero <- rowSums(assay(pe[[currentAssay]]) > 0)
pe <- zeroIsNA(pe, currentAssay) # convert 0 to NA
logTrans
if(logTrans) {
  pe <- logTransform(pe, base = 2, i = currentAssay, name = "logFeatures")
  currentAssay <- "logFeatures"
}

removeRazor
if (removeRazor) pe[[currentAssay]] <-
 pe[[currentAssay]][rowData(pe[[currentAssay]])[,groupBy]
 %in% smallestUniqueGroups(rowData(pe[[currentAssay]])[,groupBy]),]

filterColumns
if (length(filterColumns)>0)
  for (j in filterColumns)
  {
      rowData(pe[[currentAssay]])[is.na(rowData(pe[[currentAssay]])[,j]),j] <-""
      pe[[currentAssay]] <- pe[[currentAssay]][rowData(pe[[currentAssay]])[,j] != "+", ]
  }

minObsFeat
pe[[currentAssay]] <- pe[[currentAssay]][rowData(pe[[currentAssay]])$nNonZero >= minObsFeat, ]
nrow(pe[[currentAssay]])

normMethod
if (normMethod != "none")
{
  pe <- normalize(pe, 
                i = currentAssay, 
                name = "normFeatures", 
                method = normMethod)
  currentAssay <- "normFeatures"
}

pe[[currentAssay]] %>% assay %>% as.data.frame() %>%
	gather(sample, intensity) %>% 
	ggplot(aes(x = intensity, group = sample, color = sample)) +
	geom_density()

pe[[currentAssay]] %>% 
  assay %>%
  limma::plotMDS()

#	Summerization
sumMethod = "robust"

if (sumMethod!="none") 
{
  if (sumMethod=="robust") fun <- MsCoreUtils::robustSummary
  if (sumMethod=="sum") fun <- base::colSums
  if (sumMethod=="mean") fun <- base::colMeans
  if (sumMethod=="median") fun <- matrixStats::colMedians
  if (sumMethod=="medpolish") fun <- MsCoreUtils::medianPolish

  pe <- aggregateFeatures(pe,
  i = currentAssay,
  fcol = groupBy,
  na.rm = TRUE,
  name = "sumFeatures",
  fun = fun
  )
  currentAssay <- "sumFeatures"
}

plotMDS(assay(pe[[currentAssay]]))

write.csv(assay(pe[[currentAssay]]),"msqrob2_robust.txt", row.names = TRUE)

#	Data analysis
form = "~treatment * time"
doRidge = FALSE
as.formula(form)
pe <- msqrob(object = pe, i = currentAssay, formula = as.formula(form), ridge = doRidge)

#	Inference
library(ExploreModelMatrix)
visDesign <- VisualizeDesign(colData(pe),form)
visDesign$plotlist

contrast = "treatmentLN+ 0.5 treatmentLN:time2w + 0.5 treatmentLN:time8w=0"
L <- makeContrast(
  contrast, 
  parameterNames = colnames(visDesign[[3]])
  )

pe <- hypothesisTest(object = pe, i = currentAssay, contrast = L)
sigLevel = 0.05

rowData(pe[[currentAssay]])[,colnames(L)] %>% 
  arrange(pval) %>% 
  filter(adjPval < sigLevel) %>% 
  DT::datatable() %>% 
  DT::formatSignif(columns = 1:6,digits=3)

#	Plots
sigLevel
volcano <- ggplot(rowData(pe[[currentAssay]])[, colnames(L)],
                  aes(x = logFC, y = -log10(pval), color = adjPval < sigLevel)) +
  geom_point(cex = 2.5) +
  scale_color_manual(values = alpha(c("black", "red"), 0.5)) + theme_minimal()
volcano

sigNames <- rowData(pe[[currentAssay]])[,colnames(L)] %>%
  rownames_to_column("feature") %>%
  arrange(pval) %>% 
  filter(adjPval<sigLevel) %>%
  pull(feature)
if (length(sigNames) >1) heatmap(assay(pe[[currentAssay]])[sigNames, ]) else cat("No plots are generated because there are no significant summarized features at the", sigLevel, "FDR level")

