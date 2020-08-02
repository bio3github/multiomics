library(MultiAssayExperiment)
library(snpStats)
library(SNPRelate)
library(limma)
library(minfi)
library(ggpubr)

# Loading clinical data
phenotypes <- read.table("dex_stimulation_phenotypes.txt", header=T, sep="\t", check.names=F, stringsAsFactors = F)
rownames(phenotypes) <- phenotypes$Sample_ID

## Genotype Preprocessing
# Loading genotyping data
dex.imputed.genotypes <- read.plink("dex_imputed_bggt_all_qc.bed", "dex_imputed_bggt_all_qc.bim", "dex_imputed_bggt_all_qc.fam", na.strings = ("-9"))

# Loading pedigree information 
genotype.colData <- dex.imputed.genotypes$fam
available.dna.samples <- phenotypes[!is.na(phenotypes$DNA_ID), c("DNA_ID")]
genotype.colData <- genotype.colData[match(available.dna.samples, genotype.colData$member),]

# Obtain the genotype matrix
genotype <- dex.imputed.genotypes$genotype
genotype.summary <- col.summary(genotype)
# filter out SNPs with MAF less than 0.05 and call rate less than 0.95
genotype <- genotype[, with(genotype.summary, (!is.na(MAF) & MAF > 0.05) & Call.rate >= 0.95)]
# prune using ld threshold 0.9
snpgdsBED2GDS("dex_imputed_bggt_all_qc.bed",  "dex_imputed_bggt_all_qc.fam","dex_imputed_bggt_all_qc.bim", "dex_imputed_bggt_all_qc.gds")
genotype.gds <- snpgdsOpen( "dex_imputed_bggt_all_qc.gds", readonly = FALSE)
genotype.pruned <- snpgdsLDpruning(genotype.gds, ld.threshold = 0.9, snp.id = colnames(genotype)) # Only analyze the filtered SNPs
genotype.pruned.ids <- unlist(genotype.pruned, use.names=FALSE)
snpgdsClose(genotype.gds)
 
genotype.assay <-  t(genotype@.Data)[genotype.pruned.ids, available.dna.samples]

#associate SNPs with the status of patients to sort them by p-value
available.dna.samples.phenotype <- phenotypes[phenotypes$DNA_ID %in%available.dna.samples, ]
genotype.association <- data.frame(SNP=character(), P.Value=numeric(), stringsAsFactors=F)

for(i in 1:nrow(genotype.assay)) {
  model <- glm(Status ~ genotype + Sex + Age, family = binomial(), data=cbind(genotype=as.numeric(genotype.assay[i,]), available.dna.samples.phenotype))
  genotype.association[nrow(genotype.association)+1,] <- c(rownames(genotype.assay)[i], as.numeric(summary(model)[["coefficients"]][2,4]))
}
genotype.association$P.Value <- as.numeric(genotype.association$P.Value)
genotype.association$adj.P.Val <- p.adjust(genotype.association$P.Value, method = "BH", n = nrow(genotype.association))
genotype.association <- genotype.association[order(genotype.association$P.Value),]
genotype.assay <- genotype.assay[match(genotype.association$SNP, rownames(genotype.assay)),]
 
#Obtain the SNP information from geno list
genotype.rowData <- dex.imputed.genotypes$map[match(genotype.association$SNP, rownames(dex.imputed.genotypes$map)),]
colnames(genotype.rowData) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
genotype.rowData <- merge(genotype.rowData, genotype.association, by="SNP", sort=F)
rownames(genotype.rowData) <- genotype.rowData$SNP

genotype.rowRanges <-  GenomicRanges::GRanges(
  seqnames=genotype.rowData$chr,
  ranges=IRanges(start=genotype.rowData$position, end=genotype.rowData$position),
  A=genotype.rowData$A1,
  B=genotype.rowData$A2,
  pvalue=genotype.rowData$P.Value,
  adj.pvalue=genotype.rowData$adj.P.Val
) 
names(genotype.rowRanges) <- genotype.rowData$SNP

genotype.experiment <- SummarizedExperiment(genotype.assay,  rowRanges=genotype.rowRanges, colData=genotype.colData)
genotype.map <- data.frame(primary=rownames(phenotypes[!is.na(phenotypes$DNA_ID),]), colname=colnames(genotype.assay), stringsAsFactors = FALSE) 

## Transcriptom Preprocessing
load(file="ESet_all_residuals.rda") 

#Obtain the clinical data
transcriptome.colData <- pData(ESetNorm)

transcriptome.assay <- exprs(ESetNorm)

transcriptome.annotations <- read.csv(file = "GPL6947-13512.txt", sep = "\t", dec = ".", comment.char = "#",  stringsAsFactors = F) 
transcriptome.annotations <- transcriptome.annotations[match(rownames(transcriptome.assay), rownames(transcriptome.annotations)),]
transcriptome.rowData <- transcriptome.annotations[!is.na(transcriptome.annotations$start), ]
transcriptome.assay <- transcriptome.assay[match(rownames(transcriptome.rowData), rownames(transcriptome.assay)),]

transcriptome.rowRanges <- makeGRangesFromDataFrame(transcriptome.rowData, keep.extra.columns=TRUE)

dex.stimulation.transcriptome.experiment <- ExpressionSet(transcriptome.assay, featureData =  new("AnnotatedDataFrame", data=transcriptome.rowData), phenoData = new("AnnotatedDataFrame", data=transcriptome.colData))

stimulation <- factor(transcriptome.colData$Dex, levels =c(0,1))

design <- model.matrix(~0+stimulation)
colnames(design) <- c("Unstimulated", "Stimulated")

fit <- lmFit(dex.stimulation.transcriptome.experiment, design)
fit <- eBayes(fit)

contrast.matrix <- makeContrasts(StimulatedvsUnstimulatedT=Stimulated-Unstimulated, levels=design)
contrast.fit <- contrasts.fit(fit, contrast.matrix)
contrast.fit <- eBayes(contrast.fit)

dex.stimulation.dea.results <- topTable(contrast.fit, number = Inf, adjust.method="BH")

transcriptome.rowRanges <- makeGRangesFromDataFrame(dex.stimulation.dea.results, keep.extra.columns=TRUE)

dex.unstimulated.colData <- transcriptome.colData[transcriptome.colData$Dex==0,]
dex.unstimulated.colData <- dex.unstimulated.colData[match(phenotypes$Sample_ID, dex.unstimulated.colData$Sample_ID),]

dex.stimulated.colData <- transcriptome.colData[transcriptome.colData$Dex==1,]
dex.stimulated.colData <- dex.stimulated.colData[match(phenotypes$Sample_ID, dex.stimulated.colData$Sample_ID),]

dex.unstimulated.transcriptome.assay <- transcriptome.assay[match(dex.stimulation.dea.results$id, rownames(transcriptome.assay)), dex.unstimulated.colData$RNA_ID]
dex.stimulated.transcriptome.assay <- transcriptome.assay[match(dex.stimulation.dea.results$id, rownames(transcriptome.assay)), dex.stimulated.colData$RNA_ID]

unstimulated.transcriptome.experiment <- SummarizedExperiment(dex.unstimulated.transcriptome.assay,  rowRanges=transcriptome.rowRanges, colData=dex.unstimulated.colData)
unstimulated.transcriptome.map <- data.frame(primary=dex.unstimulated.colData$Sample_ID, colname=dex.unstimulated.colData$RNA_ID, stringsAsFactors = FALSE) 

stimulated.transcriptome.experiment <- SummarizedExperiment(dex.stimulated.transcriptome.assay,  rowRanges=transcriptome.rowRanges, colData=dex.stimulated.colData)
stimulated.transcriptome.map <- data.frame(primary=dex.stimulated.colData$Sample_ID, colname=dex.stimulated.colData$RNA_ID, stringsAsFactors = FALSE) 

## Methylation Preprocessing
load(file="Meth_baseline_qc_230.rda") 

methylation.colData <- pData(MethSet_base_final)
available.methylation.ids <- phenotypes[!is.na(phenotypes$Meth_ID_Base), c("Meth_ID_Base")]
methylation.colData <- methylation.colData[match(available.methylation.ids, methylation.colData$Meth_ID_Base),]

#Obtain the beta values matrix
methylation.assay <- exprs(MethSet_base_final)
methylation.assay <- methylation.assay[, available.methylation.ids]

methylation.dmp <- dmpFinder(methylation.assay, pheno = methylation.colData$Status, type = "categorical")

methylation.assay <- methylation.assay[match(rownames(methylation.dmp), rownames(methylation.assay)), ]

methylation.annotation <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
names(methylation.annotation) <- c("array", "annotation")
methylation.experiment <- RatioSet(Beta = methylation.assay, annotation=methylation.annotation, colData=methylation.colData, rowData=methylation.dmp)
methylation.map <- data.frame(primary=rownames(phenotypes[!is.na(phenotypes$Meth_ID_Base),]), colname=colnames(methylation.assay), stringsAsFactors = FALSE) 

dex.stimulation.experiments <- list("genotype"=genotype.experiment, "unstimulated.transcriptome"=unstimulated.transcriptome.experiment, "stimulated.transcriptome"=stimulated.transcriptome.experiment, "methylation"=methylation.experiment)

dex.stimulation.colData <- phenotypes

dex.stimulation.listmap <- list(genotype.map, unstimulated.transcriptome.map, stimulated.transcriptome.map, methylation.map)
names(dex.stimulation.listmap) <- c("genotype", "unstimulated.transcriptome", "stimulated.transcriptome", "methylation")
dex.stimulation.sampleMap <- listToMap(dex.stimulation.listmap)

dex.stimulation.prep.MultiAssay <- prepMultiAssay(dex.stimulation.experiments, dex.stimulation.colData, dex.stimulation.sampleMap)

dex.stimulation.MultiAssay <- MultiAssayExperiment(dex.stimulation.prep.MultiAssay$experiments, dex.stimulation.prep.MultiAssay$colData,
                                                   dex.stimulation.prep.MultiAssay$sampleMap, dex.stimulation.prep.MultiAssay$metadata)

summary(complete.cases(dex.stimulation.MultiAssay))
upsetSamples(dex.stimulation.MultiAssay)

save(dex.stimulation.MultiAssay, file = "dex.stimulation.MultiAssay.rda")

