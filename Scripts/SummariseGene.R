rsemDir = Sys.getenv("rsemDir")
workDir = Sys.getenv("diffExpDir")
inputSRAFile = Sys.getenv("inputSRAFile")
functions = Sys.getenv("functions")
project = Sys.getenv("project")
tx2gene = Sys.getenv("tx2gene")
contrast = strsplit(Sys.getenv("contrast"), split=",")

filterGenes = Sys.getenv("filterGenes")
genefilter = Sys.getenv("genefilter")

source(functions)
setwd(workDir)
require(tximport)
require(readr)


# Prepare data based on whether you want genes filtered or not ------------
pheno<-read.table(inputSRAFile, stringsAsFactors = F)
colnames(pheno) <- c("accIDs","sampleNames","phenotype" )
tx2gene<-read.table(tx2gene, stringsAsFactors = F)


if (filterGenes == F ) {
  isoforms.results <- ".isoforms.results"
  filtering = "noF"
} else {
  #Get counts from filtered files.
  isoforms.results <- ".isoforms.results.filt"
  filtering = "geneFilter"
  geneFilter <- read.table(genefilter, header=F, stringsAsFactors = F, sep="\t")
  
  
  
  
  # Create filtered gene Count  if missing
  files <- file.path(rsemDir, pheno[,1], paste0(pheno[,1], ".genes.results" ))
  names(files) <- pheno[,1]
  
  for (i in 1:length(files) ) {
    filename.x <- paste(files[i],"filt", sep=".")
    if ( !file.exists(filename.x) ) {
      a<-read.table(files[i], sep="\t", header=T, stringsAsFactors = F)
      b<- a[!a$transcript_id %in% geneFilter$V1, ]
      write.table(b,file=filename.x, sep="\t", row.names = F)
    }
  }
  
  # Create filtered gene Count  if missing
  files <- file.path(rsemDir, pheno[,1], paste0(pheno[,1], ".isoforms.results" ))
  names(files) <- pheno[,1]
  
  for (i in 1:length(files) ) {
    filename.x <- paste(files[i],"filt", sep=".")
    if ( !file.exists(filename.x) ) {
      a<-read.table(files[i], sep="\t", header=T, stringsAsFactors = F)
      b<- a[!a$transcript_id %in% geneFilter$V1, ]
      write.table(b,file=filename.x, sep="\t", row.names = F)
    }
    
  }
}


# Gene Level Counts -------------------------------------------------------
files <- file.path(rsemDir, pheno[,1], paste0(pheno[,1], isoforms.results))
names(files) <- pheno[,1]
gene.rsem <- tximport(files
                      , type = "rsem"
                      , tx2gene = tx2gene
                      , countsFromAbundance="lengthScaledTPM"
)
genecounts<-gene.rsem$counts


# Design matrix -----------------------------------------------------
require("edgeR")
require("limma")
groups <- factor(pheno[,3],ordered=TRUE, levels=unique(pheno[,3]))
design <- model.matrix(~0 + groups); colnames(design) <- gsub("groups", "", colnames(design))
contrast.matrix <- makeContrasts(contrasts=contrast, levels=design)


#Low Count Filter -------------------------------------------------
# Calculate Estimated Cut-Off
# Length Adjusted Library Size
rawLibraryLength <- gene.rsem$length
adjustedLibraryLength <- rawLibraryLength/exp(rowMeans(log(rawLibraryLength)))

# Length Adjusted Counts
LengthAdjustedCounts <- genecounts/adjustedLibraryLength


EffectiveLibrarySize <- calcNormFactors(LengthAdjustedCounts) * colSums( LengthAdjustedCounts )
MeanLibSize <- mean(EffectiveLibrarySize) * 1e-6
MedianLibSize <- median(EffectiveLibrarySize) * 1e-6
cpmCutOff <- log2(10/MedianLibSize + 2/MeanLibSize)

#Minimum CPM cutoff
min.count = 10/MedianLibSize

#Extract Length Adjusted Raw Counts
raw<-DGEList(genecounts); 
rawCPM<-cpm(raw, log = T); colnames(rawCPM)<-pheno$sampleNames

saveRDS(raw, paste("geneCounts.dgeList.raw",filtering, project, "RDS", sep="."))

#Remove noisy low counts and normalize
keep <- filterByExpr(raw, design) #Determine which genes have stat sig large counts 


filteredCounts <- raw[keep,,keep.lib.sizes=FALSE]    				#keeps only stat sig counts
filteredCPM <- cpm(filteredCounts, log=T); colnames(filteredCPM)<-pheno$sampleNames

# TMM Normalisation
normTMMCounts <- calcNormFactors(filteredCounts, method = "TMM" ); 
normCPM<-cpm(normTMMCounts, log=T); colnames(normCPM)<-pheno$sampleNames
saveRDS(normTMMCounts, paste("geneCounts.dgeList.norm.",filtering, project, "RDS", sep="."))


# Voom Normalisation ------------------------------------------------------
v <- voom(normTMMCounts, design, plot=F, save.plot = T)
mv<-data.frame(lcpm=v$voom.xy$x,
               stdv=v$voom.xy$y)
meanLine<-v$voom.line

require(ggplot2)
mvPlot <- ggplot(mv) + 
  geom_point(aes(x=lcpm, y=stdv), size=0.01, alpha=0.5 ) +
  geom_line(aes(x=meanLine$x, y=meanLine$y), color="red") +
  theme_classic() + 
  xlab("Log2(CPM + 0.5)") +
  ylab(v$voom.xy$ylab)

mvPlotFile <- paste("mvPlot.old",filtering, project, "png", sep=".")
if (!file.exists( mvPlotFile ) || overwrite == T  )  {
  ggsave(
    mvPlotFile    
    , plot = mvPlot
    , device = "png"
    , height = 6
    , width = 8
    , scale=0.4
    , dpi=300
  )
}



# Differential Expression -------------------------------------------------
fit<-lmFit(v,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

saveRDS(fit2, paste("EbayesFit",project,"RDS", sep="."))

groupComparison <- topTable(fit2, n=Inf) 
write.table(groupComparison, file=paste("DiffExp.Global", project, "txt", sep="."), row.names = T, sep="\t")

for ( i in 1:length(contrast) ) {
  assign(contrast[[i]], topTable(fit2, n=Inf, coef=i))  
  write.table(get(contrast[[i]]), file=paste("diffExp",contrast[[i]],project, "txt",sep="."), row.names = T, sep="\t")
}


# Quality Plots
## Box Plot ----------------------------------------------------------------
raw.mat<- rawCPM; colnames(raw.mat)<-pheno$phenotype
filt.mat<-filteredCPM ; colnames(filt.mat)<-pheno$phenotype
norm.mat<-normCPM; colnames(norm.mat)<-pheno$phenotype

boxRaw <- ggBoxPlot(data = filteredCPM, title = "Unprocessed")
boxNorm <- ggBoxPlot(data = normCPM, title = "Normalised")

boxNormFile <- paste("box.norm",filtering, project, "png", sep=".")

a<- ggpubr::ggarrange(boxRaw, boxNorm
                      , ncol = 2
                      , common.legend = T
                      , legend="right") 
a<-  ggpubr::annotate_figure( a
                              , bottom=ggpubr::text_grob("Sample"
                                                         , vjust = -0.5
                                                         , hjust = 1
                                                         , size = 13)
                              , left=ggpubr::text_grob("Log2(CPM + 0.5)"
                                                       , vjust = 1
                                                       , size = 13
                                                       , rot = 90 ))

ggsave(
  paste(boxNormFile)
  , plot = a
  , device = "png"
  , height = 9
  , width = 18
  , dpi=300
  , scale=0.33
)


## Density plot -------------------------------------------------------------------------
cutoffLine <- geom_vline(xintercept = log2(min.count), color="black", linetype = "dashed" )

#Plot Density for unfiltered
unfilteredPlot <- densityPlot(rawCPM, phenoData = pheno, title="Raw", fill="phenotype") + cutoffLine
densityRaw <- densityPlot(filteredCPM, phenoData = pheno, title="Low_Count_Filtered", fill="phenotype")+ cutoffLine
densityNorm <- densityPlot(normCPM, phenoData = pheno, title="Normalised", fill="phenotype")+ cutoffLine

density2<- ggpubr::ggarrange(unfilteredPlot, densityRaw, densityNorm
                             , ncol = 3
                             , common.legend = T
                             , legend="right") 
density2<-  ggpubr::annotate_figure( density2
                                     , bottom=ggpubr::text_grob("Log2(CPM + 0.5)"
                                                                , vjust = -0.5
                                                                , hjust = 1
                                                                , size = 13)
                                     , left=ggpubr::text_grob("Density"
                                                              , vjust = 1
                                                              , size = 13
                                                              , rot = 90 ))


densityFile <- paste("density.norm.raw",filtering, project, "svg", sep=".")

ggsave(
  paste(densityFile)
  , plot = density2
  , device = "svg"
  , height = 2
  , width = 6
)


## PCA ---------------------------------------------------------------------

#PCA Plot
filtPCA <- pcaPlot(data=filt.mat, title = "Unprocessed",  shape.manual = T, shape.value = c(15,16,17,18,19,20,9), shape = colnames(filt.mat), color = colnames(filt.mat))
normPCA <- pcaPlot(data=norm.mat, title = "Normalised", shape.manual = T, shape.value = c(15,16,17,18,19,20,9),  shape = colnames(norm.mat), color = colnames(norm.mat))

pca2<- ggpubr::ggarrange(filtPCA, normPCA
                         , ncol = 2
                         , common.legend = T
                         , legend = "right"
) 
pcaFile <- paste("pca",filtering, project, "svg", sep=".")

ggsave(
  paste(pcaFile)
  , plot = pca2
  , device = "svg"
  , height = 3
  , width = 6
)
