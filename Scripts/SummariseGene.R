rsemDir = Sys.getenv("rsemDir")
workDir = Sys.getenv("diffExpDir")
inputSRAFile = Sys.getenv("inputSRAFile")
functions = Sys.getenv("functions")
project = Sys.getenv("project")
tx2gene = Sys.getenv("tx2gene")
contrast = strsplit(Sys.getenv("contrast"), split=",")
filterGenes = Sys.getenv("filterGenes")


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
  genefilter = Sys.getenv("genefilter")
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

raw<-DGEList(genecounts); 
raw <- calcNormFactors(raw); 

#Remove noisy low counts and normalize
keep <- filterByExpr(raw, design)					    #Determine which genes have stat sig large counts 
dge <- raw[keep,,keep.lib.sizes=FALSE]    				#keeps only stat sig counts
filteredExpr <- cpm(dge, log=T); colnames(filteredExpr)<-pheno$sampleNames

saveRDS(dge, paste("geneCounts.filtered", project, "RDS", sep="."))


# Voom Normalisation ------------------------------------------------------
v <- voom(dge, design, plot=T, save.plot = T)
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
va <-v$E; colnames(va) <- pheno$sampleNames
boxRaw <- ggBoxPlot(data = filteredExpr, title = "A")
boxNorm <- ggBoxPlot(data = va, title = "B")

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
                              , left=ggpubr::text_grob("Log2(gene count + 0.5)"
                                                       , vjust = 1
                                                       , size = 13
                                                       , rot = 90 ))

if (!file.exists( boxNormFile ) )  {
  ggsave(
    paste(boxNormFile)
    , plot = a
    , device = "png"
    , height = 3
    , width = 6
  )
}


## Density plot -------------------------------------------------------------------------
unfilteredExpr <- cpm(raw, log=T); colnames(unfilteredExpr)<-pheno$sampleNames

#Plot Density for unfiltered
unfilteredPlot <- densityPlot(unfilteredExpr, phenoData = pheno, title="A", fill="phenotype")
densityRaw <- densityPlot(filteredExpr, phenoData = pheno, title="B", fill="phenotype")
densityNorm <- densityPlot(va, phenoData = pheno, title="C", fill="phenotype")

density2<- ggpubr::ggarrange(unfilteredPlot, densityRaw, densityNorm
                             , ncol = 3
                             , common.legend = T
                             , legend="right") 
density2<-  ggpubr::annotate_figure( density2
                                     , bottom=ggpubr::text_grob("Log2(gene count + 0.5)"
                                                                , vjust = -0.5
                                                                , hjust = 1
                                                                , size = 13)
                                     , left=ggpubr::text_grob("Density"
                                                              , vjust = 1
                                                              , size = 13
                                                              , rot = 90 ))


densityFile <- paste("density.norm.raw",filtering, project, "svg", sep=".")
if (!file.exists( densityFile ) )  {
  ggsave(
    paste(densityFile)
    , plot = density2
    , device = "svg"
    , height = 3
    , width = 9
  )
}


## PCA ---------------------------------------------------------------------

#PCA Plot
raw.mat<-filteredExpr ; colnames(raw.mat)<-pheno$phenotype
norm.mat<-va; colnames(va)<-pheno$phenotype

rawPCA <- pcaPlot(data=raw.mat, title = "A",  shape.manual = T, shape.value = c(15,16,17,18,19,20,9), shape = colnames(raw.mat), color = colnames(raw.mat))
normPCA <- pcaPlot(data=norm.mat, title = "B", shape.manual = T, shape.value = c(15,16,17,18,19,20,9),  shape = colnames(raw.mat), color = colnames(raw.mat))

pca2<- ggpubr::ggarrange(normPCA, rawPCA
                         , ncol = 2
                         , common.legend = T
                         , legend = "right"
) 
pcaFile <- paste("pca",filtering, project, "svg", sep=".")
if (!file.exists( pcaFile ) )  {
  ggsave(
    paste(pcaFile)
    , plot = pca2
    , device = "svg"
    , height = 3
    , width = 6
  )
}
