library(MatrixEQTL)
library(minfi)
library(MEAL)
library(mixOmics)
library(org.Hs.eg.db)
library(clusterProfiler)
library(matrixStats)
library(bestNormalize)
library(car)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(plyr)
library(dplyr)
library(stringr)
library(reshape2)
library(UpSetR)

# We can load all the available datasets
load(file = "dex.stimulation.MultiAssay.rda", verbose = T)  
#subset for eQTL experiment
dex.stimulation.genotype.transcriptome.MultiAssay <- dex.stimulation.MultiAssay[, , c(3, 1)]
rm(dex.stimulation.MultiAssay)

dex.stimulation.genotype.transcriptome.colData <- as.data.frame(colData(dex.stimulation.genotype.transcriptome.MultiAssay))
transcriptome.colData <- as.data.frame(colData(experiments(dex.stimulation.genotype.transcriptome.MultiAssay)[[1]])) 

transcriptome.rowRanges <- as.data.frame(rowRanges(experiments(dex.stimulation.genotype.transcriptome.MultiAssay)[[1]])) 
transcriptome.rowRanges <- transcriptome.rowRanges[transcriptome.rowRanges$adj.P.Val <0.01,]

transcriptome.assay <- as.matrix(assays(dex.stimulation.genotype.transcriptome.MultiAssay)[[1]])
transcriptome.assay <- transcriptome.assay[match(rownames(transcriptome.rowRanges), rownames(transcriptome.assay)), match(dex.stimulation.genotype.transcriptome.colData$Stimulated_RNA_ID, colnames(transcriptome.assay))]
colnames(transcriptome.assay) <- dex.stimulation.genotype.transcriptome.colData$Sample_ID

transcriptome.slices <- SlicedData$new();
transcriptome.slices$CreateFromMatrix(transcriptome.assay) 
transcriptome.slices$fileSliceSize <- 1000; # read file in slices of 1,000 rows
transcriptome.positions <- data.frame(transcriptome=rownames(transcriptome.rowRanges), chr=transcriptome.rowRanges$seqnames, start=transcriptome.rowRanges$start, end=transcriptome.rowRanges$end)

covariates <- t(transcriptome.colData[transcriptome.colData$Dex==1, c("Status")]) 
colnames(covariates) <- transcriptome.colData$Sample_ID 
rownames(covariates)  <- c("Status")
covariates.slices = SlicedData$new();
covariates.slices$CreateFromMatrix(covariates) 
covariates.slices$fileSliceSize <- 1; # read file in slices of 1 row
covariates.slices$fileSkipRows = 1;          # one row of column labels
covariates.slices$fileSkipColumns = 1;       # one column of row labels

genotype.colData <- as.data.frame(colData(experiments(dex.stimulation.genotype.transcriptome.MultiAssay)[[2]])) 

genotype.rowRanges <- as.data.frame(rowRanges(experiments(dex.stimulation.genotype.transcriptome.MultiAssay)[[2]])) 

genotype.assay <- as.matrix(assays(dex.stimulation.genotype.transcriptome.MultiAssay)[[2]])
genotype.assay <- genotype.assay[match(rownames(genotype.rowRanges), rownames(genotype.assay)), match(dex.stimulation.genotype.transcriptome.colData$DNA_ID, colnames(genotype.assay))]

colnames(genotype.assay) <- dex.stimulation.genotype.transcriptome.colData$Sample_ID

genotype.categories <- length(unique(as.vector(genotype.assay)))

genotype.slices <- SlicedData$new() 
genotype.slices$CreateFromMatrix(genotype.assay) 
genotype.slices$fileSliceSize <- 1000  # read file in slices of 1,000 rows
genotype.positions <- data.frame(snp=rownames(genotype.rowRanges), chr=genotype.rowRanges$seqnames, pos=genotype.rowRanges$start)


options(MatrixEQTL.ANOVA.categories=genotype.categories) 
cases.vs.controls.stimulation.eQTL <- Matrix_eQTL_main(
  snps = genotype.slices,
  gene = transcriptome.slices,
  cvrt = covariates.slices,
  output_file_name = "cases.vs.controls.stimulated.trans.eQTL",
  pvOutputThreshold = 1e-5,
  useModel = modelANOVA,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = "cases.vs.controls.stimulates.cis.eQTL",
  pvOutputThreshold.cis = 1e-2,
  snpspos = genotype.positions,
  genepos = transcriptome.positions,
  cisDist = 1e6,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE)

cases.vs.controls.stimulation.cis.eqtls <- cases.vs.controls.stimulation.eQTL$cis$eqtls
names(cases.vs.controls.stimulation.cis.eqtls)[1] <- "snp_id"
names(cases.vs.controls.stimulation.cis.eqtls)[2] <- "probe_id"
names(transcriptome.rowRanges)[6] <- "probe_id"
genotype.rowRanges$snp_id <- rownames(genotype.rowRanges)
cases.vs.controls.stimulation.cis.eqtls <- merge(cases.vs.controls.stimulation.cis.eqtls, transcriptome.rowRanges, by="probe_id", all.x=T, sort=F)
names(cases.vs.controls.stimulation.cis.eqtls)[6:38] <- paste0("gene_",names(cases.vs.controls.stimulation.cis.eqtls)[6:38])
cases.vs.controls.stimulation.cis.eqtls <- merge(cases.vs.controls.stimulation.cis.eqtls, genotype.rowRanges, by="snp_id", all.x=T, sort=F)
names(cases.vs.controls.stimulation.cis.eqtls)[39:47] <- paste0("snp_",names(cases.vs.controls.stimulation.cis.eqtls)[39:47])

cases.vs.controls.stimulation.trans.eqtls <- cases.vs.controls.stimulation.eQTL$trans$eqtls
names(cases.vs.controls.stimulation.trans.eqtls)[1] <- "snp_id"
names(cases.vs.controls.stimulation.trans.eqtls)[2] <- "probe_id"

cases.vs.controls.stimulation.trans.eqtls <- merge(cases.vs.controls.stimulation.trans.eqtls, transcriptome.rowRanges, by="probe_id", all.x=T, sort=F)
names(cases.vs.controls.stimulation.trans.eqtls)[6:38] <- paste0("gene_",names(cases.vs.controls.stimulation.trans.eqtls)[6:38])
cases.vs.controls.stimulation.trans.eqtls <- merge(cases.vs.controls.stimulation.trans.eqtls, genotype.rowRanges, by="snp_id", all.x=T, sort=F)
names(cases.vs.controls.stimulation.trans.eqtls)[39:47] <- paste0("snp_",names(cases.vs.controls.stimulation.trans.eqtls)[39:47])

rm(dex.stimulation.genotype.transcriptome.MultiAssay,transcriptome.assay,genotype.assay)

load(file = "dex.stimulation.MultiAssay.rda", verbose = T)  

dex.stimulation.transcriptome.methylation.MultiAssay <- intersectColumns(dex.stimulation.MultiAssay[, , c(3, 4)])
rm(dex.stimulation.MultiAssay)

#phenotypes <- read.table("dex_stimulation_phenotypes.txt", header=T, sep="\t", check.names=F, stringsAsFactors = F)
dex.stimulation.transcriptome.methylation <- as.data.frame(colData(dex.stimulation.transcriptome.methylation.MultiAssay))
transcriptome.colData <- as.data.frame(colData(experiments(dex.stimulation.transcriptome.methylation.MultiAssay)[[1]])) 
rownames(transcriptome.colData) <- dex.stimulation.transcriptome.methylation$Meth_ID_Base

transcriptome.rowRanges <- as.data.frame(rowRanges(experiments(dex.stimulation.transcriptome.methylation.MultiAssay)[[1]])) 
names(transcriptome.rowRanges)[1] <- "chromosome"
transcriptome.assay <- as.matrix(assays(dex.stimulation.transcriptome.methylation.MultiAssay)[[1]])
 
colnames(transcriptome.assay) <- dex.stimulation.transcriptome.methylation$Meth_ID_Base
 
transcriptome.experiment <-  ExpressionSet(transcriptome.assay, featureData =  new("AnnotatedDataFrame", data=transcriptome.rowRanges), phenoData = new("AnnotatedDataFrame", data=transcriptome.colData))
methylation.experiment <-  mapToGenome(experiments(dex.stimulation.transcriptome.methylation.MultiAssay)[[2]])

transcriptome.methylation.experiment <- createMultiDataSet()
transcriptome.methylation.experiment <- add_genexp(transcriptome.methylation.experiment, transcriptome.experiment)
transcriptome.methylation.experiment <- add_methy(transcriptome.methylation.experiment, methylation.experiment)

transcriptome.rowRanges <- transcriptome.rowRanges[transcriptome.rowRanges$adj.P.Val <0.05,]
transcriptome.methylation.correlations <- data.frame(stringsAsFactors = F)
for (i in 1:nrow(transcriptome.rowRanges)) {
  targetRange <-GenomicRanges::GRanges(
    seqnames=transcriptome.rowRanges[i,c("chromosome")],
    ranges=IRanges(start=transcriptome.rowRanges[i,c("start")]-10000, end=transcriptome.rowRanges[i,c("end")]+10000),
  ) 
  names(targetRange) <- transcriptome.rowRanges[i,c("id")]
  
  ranged.transcriptome.methylation.experiment <- transcriptome.methylation.experiment[, , targetRange]
  transcriptome.methylation.correlations <- rbind(transcriptome.methylation.correlations,  correlationMethExprs(ranged.transcriptome.methylation.experiment))
  
}

rm(dex.stimulation.transcriptome.methylation.MultiAssay,transcriptome.assay, transcriptome.colData, transcriptome.experiment, methylation.experiment, transcriptome.methylation.experiment, ranged.transcriptome.methylation.experiment)



load(file = "dex.stimulation.MultiAssay.rda", verbose = T)  

normalize.assays <- function(dex.stimulation.MultiAssay, genotype.pvalue, transcriptome.pvalue, methylation.pvalue) {
  
  dex.stimulation.MultiAssay.colData <- as.data.frame(colData(dex.stimulation.MultiAssay))
  
  genotype.colData <- as.data.frame(colData(experiments(dex.stimulation.MultiAssay)[[1]])) 
  
  genotype.rowRanges <- as.data.frame(rowRanges(experiments(dex.stimulation.MultiAssay)[[1]])) 
  genotype.rowRanges <- genotype.rowRanges[genotype.rowRanges$pvalue < genotype.pvalue,]
  
  genotype.assay <- as.matrix(assays(dex.stimulation.MultiAssay)[[1]])
  
  genotype.assay <- genotype.assay[match(rownames(genotype.rowRanges), rownames(genotype.assay)), match(dex.stimulation.MultiAssay.colData$DNA_ID, colnames(genotype.assay))]
  colnames(genotype.assay) <- dex.stimulation.MultiAssay.colData$Sample_ID

  normalized.genotype.assay <- matrix(as.numeric(genotype.assay), nrow=nrow(genotype.assay), ncol=ncol(genotype.assay))
  
  colnames(normalized.genotype.assay) <- colnames(genotype.assay)
  rownames(normalized.genotype.assay) <- paste0("geno", ":", rownames(genotype.assay))
   
  unstimulated.transcriptome.colData <- as.data.frame(colData(experiments(dex.stimulation.MultiAssay)[[2]])) 
  
  unstimulated.transcriptome.rowRanges <- as.data.frame(rowRanges(experiments(dex.stimulation.MultiAssay)[[2]])) 
  unstimulated.transcriptome.rowRanges <- unstimulated.transcriptome.rowRanges[unstimulated.transcriptome.rowRanges$P.Value < transcriptome.pvalue,]
  
  unstimulated.transcriptome.assay <- as.matrix(assays(dex.stimulation.MultiAssay)[[2]])
  unstimulated.transcriptome.assay <- unstimulated.transcriptome.assay[match(rownames(unstimulated.transcriptome.rowRanges), rownames(unstimulated.transcriptome.assay)), match(dex.stimulation.MultiAssay.colData$Unstimulated_RNA_ID, colnames(unstimulated.transcriptome.assay))]
  all(unstimulated.transcriptome.colData$RNA_ID==colnames(unstimulated.transcriptome.assay))
  unstimulated.transcriptome.assay.orderNorm <- orderNorm(as.vector(unstimulated.transcriptome.assay))
  normalized.unstimulated.transcriptome.assay <- matrix(unstimulated.transcriptome.assay.orderNorm[[1]], nrow = nrow(unstimulated.transcriptome.assay), ncol = ncol(unstimulated.transcriptome.assay),  dimnames = list(rownames(unstimulated.transcriptome.assay), colnames(unstimulated.transcriptome.assay)))
  
  colnames(normalized.unstimulated.transcriptome.assay) <- dex.stimulation.MultiAssay.colData$Sample_ID
  rownames(normalized.unstimulated.transcriptome.assay) <- paste0("unstim.trans", ":", rownames(normalized.unstimulated.transcriptome.assay))
 
  print(paste0("Consistent genotype and unstimulated transcriptome Samples: " , all(colnames(normalized.unstimulated.transcriptome.assay)==colnames(normalized.genotype.assay))))
 
  stimulated.transcriptome.colData <- as.data.frame(colData(experiments(dex.stimulation.MultiAssay)[[3]])) 
  
  stimulated.transcriptome.rowRanges <- as.data.frame(rowRanges(experiments(dex.stimulation.MultiAssay)[[3]])) 
  stimulated.transcriptome.rowRanges <- stimulated.transcriptome.rowRanges[stimulated.transcriptome.rowRanges$P.Value < transcriptome.pvalue,]
  
  stimulated.transcriptome.assay <- as.matrix(assays(dex.stimulation.MultiAssay)[[3]])
  stimulated.transcriptome.assay <- stimulated.transcriptome.assay[match(rownames(stimulated.transcriptome.rowRanges), rownames(stimulated.transcriptome.assay)), match(dex.stimulation.MultiAssay.colData$Stimulated_RNA_ID, colnames(stimulated.transcriptome.assay))]
  all(stimulated.transcriptome.colData$RNA_ID == colnames(stimulated.transcriptome.assay))
  
  stimulated.transcriptome.assay.orderNorm <- orderNorm(as.vector(stimulated.transcriptome.assay))
  normalized.stimulated.transcriptome.assay <- matrix(stimulated.transcriptome.assay.orderNorm[[1]], nrow = nrow(stimulated.transcriptome.assay), ncol = ncol(stimulated.transcriptome.assay),  dimnames = list(rownames(stimulated.transcriptome.assay), colnames(stimulated.transcriptome.assay)))
  
  colnames(normalized.stimulated.transcriptome.assay) <- dex.stimulation.MultiAssay.colData$Sample_ID
  rownames(normalized.stimulated.transcriptome.assay) <- paste0("stim.trans", ":", rownames(normalized.stimulated.transcriptome.assay))
  
  print(paste0("Consistent unstimulated and stimulated transcriptome Samples: " , all(colnames(normalized.unstimulated.transcriptome.assay)==colnames(normalized.stimulated.transcriptome.assay))))

  methylation.colData <- as.data.frame(colData(experiments(dex.stimulation.MultiAssay)[[4]])) 
  
  methylation.rowData <- as.data.frame(rowData(experiments(dex.stimulation.MultiAssay)[[4]])) 
  methylation.rowData <- methylation.rowData[methylation.rowData$pval < methylation.pvalue,]
  
  methylation.assay <- as.matrix(assays(dex.stimulation.MultiAssay)[[4]])
  methylation.assay <- methylation.assay[match(rownames(methylation.rowData), rownames(methylation.assay)), match(dex.stimulation.MultiAssay.colData$Meth_ID_Base, colnames(methylation.assay))]
  all(methylation.colData$Meth_ID_Base == colnames(methylation.assay))
  
  methylation.assay.orderNorm <- orderNorm(as.vector(methylation.assay))
  normalized.methylation.assay <- matrix(methylation.assay.orderNorm[[1]], nrow = nrow(methylation.assay), ncol = ncol(methylation.assay),  dimnames = list(rownames(methylation.assay), colnames(methylation.assay)))
  
  colnames(normalized.methylation.assay) <- dex.stimulation.MultiAssay.colData$Sample_ID
  rownames(normalized.methylation.assay) <- paste0("meth", ":", rownames(normalized.methylation.assay))
  
  print(paste0("Consistent stimulated transcriptome and methylation Samples: " , all(colnames(normalized.methylation.assay)==colnames(normalized.stimulated.transcriptome.assay))))
  list(normalized.genotype.assay, normalized.unstimulated.transcriptome.assay, normalized.stimulated.transcriptome.assay, normalized.methylation.assay)
}

unsupervisedSPLS <- function(normalized.genotype.assay, normalized.methylation.assay, normalized.transcriptome.assay, ncomp) {
  
  X <- list(geno=t(normalized.genotype.assay), meth=t(normalized.methylation.assay))
  design = matrix(1, ncol = length(X), nrow = length(X), dimnames = list(names(X), names(X)))
  diag(design) =  0
  list.keepX = list(geno = rep(25, ncomp), meth = rep(25, ncomp))
  
  Y <-  t(normalized.transcriptome.assay)
  list.keepY = c(rep(100, ncomp))
  
  block.pls.res <- block.spls(X = X, Y = Y , ncomp = ncomp, keepX = list.keepX, keepY = list.keepY, design = design)
  
}
# Functional Enticements of latent structures models

# KEGG Pathways Enrichment
enrichKEGGPathways <- function(expressions.factor) {
  entrez.ids.factor <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=expressions.factor, column="ENTREZID", keytype="SYMBOL")  
  kegg.factor <- enrichKEGG(gene = entrez.ids.factor, organism='hsa', keyType='ncbi-geneid', pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 2)
  return (kegg.factor@result)
}

getKEGGPathwaysResult <- function(block.pls.res, ncomp, transcriptome.rowRanges) {
  kegg.enriched.pathways.statistics <- data.frame(factor=integer(), pathways.count=integer(), pathways.ids=character(),
                                                  genes.count=integer(), genes.symbols=character(), stringsAsFactors=F)
  
  kegg.enriched.pathways.summary <-  data.frame(factor=integer(), pathway.id=character(), pathway.description=character(), pvalue=numeric(), adjusted.pvalue=numeric(),
                                                gene.ratio=character(), genes.symbols=character(), stringsAsFactors=F)
  kegg.enriched.pathways <- list() 
  
  for(i in 1:ncomp) {
    print(paste0("Processing factor:",i))
    probe.ids <- ldply(strsplit(selectVar(block.pls.res, block =  'Y', comp =i)[[1]][[1]],":"))[,2]
    probe.gene.symols <-unique(transcriptome.rowRanges[rownames(transcriptome.rowRanges) %in% probe.ids, "symbol"])
    kegg.factor.result <- enrichKEGGPathways(probe.gene.symols)
    kegg.factor.significant.result <- kegg.factor.result[kegg.factor.result$p.adjust <0.05,]
    if(nrow(kegg.factor.significant.result) > 0) {
      for(j in 1:nrow(kegg.factor.significant.result)) {
        
        pathways.genes.ids <- unlist(strsplit(kegg.factor.significant.result[j,8], "/"))
        if(length(pathways.genes.ids) >0) {
          pathway.genes.symbols <-  AnnotationDbi::mapIds(org.Hs.eg.db, keys=pathways.genes.ids, column="SYMBOL", keytype="ENTREZID")  
        } else {
          pathway.genes.symbols <- character()
        }
        pathway.genes.symbols  <- paste(pathway.genes.symbols, collapse = "/")
        
        kegg.enriched.pathways.summary[nrow(kegg.enriched.pathways.summary) + 1,] = c(i, kegg.factor.significant.result[j,1], kegg.factor.significant.result[j,2],
                                                                                      kegg.factor.significant.result[j,5], kegg.factor.significant.result[j,6],
                                                                                      kegg.factor.significant.result[j,3], pathway.genes.symbols)
      }
    }
    kegg.pathways.ids <- paste(kegg.factor.significant.result$ID, collapse = "/")
    kegg.pathways.genes.count <- length(unique(strsplit(paste(kegg.factor.significant.result$geneID, collapse = "/") , "/") [[1]])) 
    kegg.pathways.genes.ids <-  paste(unique(strsplit(paste(kegg.factor.significant.result$geneID, collapse = "/") , "/") [[1]]), collapse = "/") 
    factor.genes.ids <- unique(unlist(strsplit(kegg.factor.significant.result$geneID,"/")))
    
    if(length(factor.genes.ids) > 0) {
      factor.genes.symbols <-  AnnotationDbi::mapIds(org.Hs.eg.db, keys=factor.genes.ids, column="SYMBOL", keytype="ENTREZID")  
    } else {
      factor.genes.symbols <- character()
    }
    kegg.factor.pathways.genes.symbols <-  paste(factor.genes.symbols, collapse = "/")
    
    kegg.enriched.pathways.statistics[i,] <- c(i, nrow(kegg.factor.significant.result), kegg.pathways.ids, kegg.pathways.genes.count, kegg.factor.pathways.genes.symbols)
    kegg.enriched.pathways[[i]] <- kegg.factor.significant.result
  }
  return (list(kegg.enriched.pathways.statistics, kegg.enriched.pathways.summary, kegg.enriched.pathways))
}



dex.stimulation.MultiAssay <-  intersectColumns(dex.stimulation.MultiAssay)
normalized.dex.stimulation.MultiAssay <- normalize.assays(dex.stimulation.MultiAssay, 0.05, 0.05, 0.05)

stimulated.spls <- unsupervisedSPLS(normalized.dex.stimulation.MultiAssay[[1]], normalized.dex.stimulation.MultiAssay[[4]], normalized.dex.stimulation.MultiAssay[[3]], 20)

stimulated.transcriptome.rowRanges <- as.data.frame(rowRanges(experiments(dex.stimulation.MultiAssay)[[3]])) 

stimulated.kegg.enriched.pathways.result <- getKEGGPathwaysResult(stimulated.spls, 20, stimulated.transcriptome.rowRanges)
stimulated.kegg.enriched.pathways.statistics <- stimulated.kegg.enriched.pathways.result[[1]]
stimulated.kegg.enriched.pathways.summary <- stimulated.kegg.enriched.pathways.result[[2]]
stimulated.kegg.enriched.pathways.list<- stimulated.kegg.enriched.pathways.result[[3]]

stimulated.kegg.enriched.pathways.statistics$factor <- as.integer(stimulated.kegg.enriched.pathways.statistics$factor)
stimulated.kegg.enriched.pathways.statistics$pathways.count <- as.integer(stimulated.kegg.enriched.pathways.statistics$pathways.count)
stimulated.kegg.enriched.pathways.statistics$genes.count <- as.integer(stimulated.kegg.enriched.pathways.statistics$genes.count)

stimulated.kegg.enriched.pathways.statistics.figure <- ggbarplot(stimulated.kegg.enriched.pathways.statistics, x = "factor", y = "pathways.count", fill =  "genes.count", label = TRUE, label.pos = "out", 
                                                                           legend.title = "Expressed Genes", xlim = c(1, 20), xticks.by=1, xlab ="Factor", ylab = "KEGG Pathways")

print(stimulated.kegg.enriched.pathways.statistics.figure)

omics.factor.functional.genes <- str_split_fixed(stimulated.kegg.enriched.pathways.statistics[3, c("genes.symbols")], "/", Inf) 

# extract correlations per component
extractFactorCorrelation <- function (omics, block.pls.res, comp) {
  # here you may normalize and centralize the omics matrix  
  # extract loading of a specific component
  omics.var.names <- character()
  factor.vars <- selectVar(block.pls.res,  comp = comp) 
  
  for (i in 1:(length(factor.vars)-1)) {
    omics.var.names <- c(omics.var.names, factor.vars[[i]]$name)
    
  }
  
  omics.factor.vars <- omics[, colnames(omics) %in% omics.var.names]
  dim(omics.factor.vars)
  omics.factor.vars.pairs <- data.frame(t(combn( colnames(omics.factor.vars), 2)), stringsAsFactors = T ) 
  # apply the correlations functions for the selected assays
  omics.factor.vars.correlations <- do.call(rbind, mapply(function(var1, var2, data) {
    result = cor.test(data[,var1], data[,var2], method="pearson")
    data.frame(var1, var2, result[c("estimate","p.value","statistic")], stringsAsFactors=FALSE)
  }, omics.factor.vars.pairs[,1], omics.factor.vars.pairs[,2], MoreArgs=list(data=omics.factor.vars), SIMPLIFY=FALSE))
  omics.factor.vars.correlations[,c("assay1", "var1")] <- str_split_fixed(omics.factor.vars.correlations$var1, ":", 2) 
  omics.factor.vars.correlations[,c("assay2", "var2")] <- str_split_fixed(omics.factor.vars.correlations$var2, ":", 2) 
  omics.factor.vars.correlations <- omics.factor.vars.correlations[order(omics.factor.vars.correlations$p.value), ]
  
}
plotFactorCorrelations <- function(factor.correlations, comp) {
  
  factor.correlations.matrix <- dcast(factor.correlations, var1~var2, value.var='estimate')
  factor.correlations.matrix[is.na(factor.correlations.matrix)] <- 0
  rownames(factor.correlations.matrix) <- factor.correlations.matrix[,1]
  factor.correlations.matrix <- t(as.matrix(factor.correlations.matrix[-1]))
  dim(factor.correlations.matrix)
  
  annotation.row <- data.frame(factor.correlations %>% group_by(var2)  %>%  summarise(omics=unique(assay2)))
  rownames(annotation.row) <- annotation.row$var2
  
  main.title <- paste("Factor", comp, "Omics Correlations")
  pheatmap::pheatmap(factor.correlations.matrix,
                     cluster_cols = TRUE,
                     cluster_rows = TRUE,
                     annotation_row = annotation.row,
                     annotation_legend = FALSE,
                     show_colnames=FALSE,
                     main = main.title,
                     clustering_distance_rows = "manhattan",
                     fontsize_row = 1)
}

omics <- cbind(t(normalized.dex.stimulation.MultiAssay[[1]]), t(normalized.dex.stimulation.MultiAssay[[4]]), t(normalized.dex.stimulation.MultiAssay[[3]]))
 
#Extract and plot factor correlations for factor 3.
omics.factor.vars.correlations <- extractFactorCorrelation(omics, stimulated.spls, 3)
omics.factor.vars.correlations$var1.genomic.region <- NA
omics.factor.vars.correlations$var2.genomic.region <- NA
plotFactorCorrelations(omics.factor.vars.correlations, 3)


