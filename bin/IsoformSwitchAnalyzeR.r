                                            #########################
                                            # IsoformSwitchAnalyzeR #
                                            #########################


##########################################################################################################

# Setup R Environment

### Load library
library(IsoformSwitchAnalyzeR, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)

### Check that you're using yhe ritght version
packageVersion("IsoformSwitchAnalyzeR")
packageVersion("dplyr")

# Set working directory
setwd("/groups/dog/stage/enora/isoformswitch/")

##########################################################################################################

args = commandArgs(trailingOnly=TRUE)

design          <- strsplit(grep('--tag*', args, value = TRUE), split = '=')[[1]][[2]]
counts          <- strsplit(grep('--ncore*', args, value = TRUE), split = '=')[[1]][[2]]
exon_annot      <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]
transcriptome   <- strsplit(grep('--fasta*', args, value = TRUE), split = '=')[[1]][[2]]

# Experimental Data
design <- as.data.frame.matrix(design)
cancerTypes <- unique(design$cancer)
design <- unite(design, condition, cancer, condition, sep = "_")
design <- subset(design, select = c(sample_id, condition))
names(design)[names(design) == "sample_id"] <- "sampleID"

# Header : samples names to match between the design file and the counts file
samples <- (c(design[,1]))
samples=sort(samples)
samples[1:2] = rev(samples[1:2])

# Isoforms counts
counts <- as.data.frame.matrix(counts)
counts <- subset(counts, select = -c(GENEID))
colnames(counts) <- c("isoform_id",samples)


# Load Data


# Experimental Data
design <- read.table("/groups/dog/nanopore/lncrna_resist_cgo/secondary/design_lncrnaResist_allSamples.csv", header=TRUE, sep=',')
counts <- read.table("/groups/dog/stage/enora/isoformswitch/Data/counts_transcript_filter.full.txt", header=TRUE)
exon_annot <- "/groups/dog/stage/enora/isoformswitch/Data/extended_annotations_filter.full.gtf"
transcriptome <- "/groups/dog/stage/enora/isoformswitch/Data/transcriptome_lncrna-resist.fa"


# Comparisons table
comparisons <- data.frame(
    condition_1 = c("melanoma_sensitive", "glioblastoma_sensitive", "prostate_cancer_sensitive", "lung_cancer_sensitive"),
    condition_2 = c("melanoma_resistant", "glioblastoma_resistant", "prostate_cancer_resistant", "lung_cancer_resistant")
    )

# Transcriptome (FASTA)

message('Importing data')
#if (file.exists("SwitchList_lncrna_resist.R")==TRUE){
#    SwitchList <- readRDS(file = "SwitchList_lncrna_resist.Rds")
#}
#else {
    # Import data
SwitchList <- importRdata(
    ### Core arguments
    isoformCountMatrix = counts,
    isoformRepExpression = NULL,
    designMatrix = design,
    isoformExonAnnoation= exon_annot,
    isoformNtFasta = transcriptome,
    comparisonsToMake = comparisons,
    showProgress = FALSE,
    quiet = TRUE
)
saveRDS(SwitchList, file = "SwitchList_lncrna_resist.Rds")

SwitchList

isoforms = SwitchList$isoformFeatures[!duplicated(SwitchList$isoformFeatures$isoform_id), ]
biotypes_total <- data.frame(table(isoforms$iso_biotype))
colnames(biotypes_total) <- c('biotype', 'frequency')
biotypes_total$biotype <- factor(biotypes_total$biotype, levels = biotypes_total$biotype[order(biotypes_total$frequency)])
plotbiotypes <- ggplot(
    data=biotypes_total,
    aes(x=biotype, y=frequency, fill=biotype)) +
    geom_bar(stat="identity") +
    coord_flip() +
    geom_text(
        aes(label = frequency),
        position = position_stack(vjust = 0.5),
        size=3) +
    scale_fill_hue(direction=-1) +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle=90),
        legend.position="none")



for (cancer in cancerTypes) {

    message(paste0(cancer, 'subseting'))
    subsetSwitchGB <- subsetSwitchAnalyzeRlist(
        switchAnalyzeRlist=SwitchList,
        SwitchList$isoformFeatures$condition_1 %in% c(
            paste0(cancer, "_sensitive"),
            paste0(cancer, "_resistant")
        ) &
        SwitchList$isoformFeatures$condition_2 %in% c(
            paste0(cancer, "_sensitive"),
            paste0(cancer, "_resistant")
        )
    )

    message(paste0(cancer, 'filtering'))
    subsetSwitchGB <- preFilter(
        switchAnalyzeRlist= subsetSwitchGB,
        geneExpressionCutoff = 1,
        isoformExpressionCutoff = 0,
        quiet=TRUE
    )
    saveRDS(subsetSwitchGB, file = paste0(cancer,"_lncrna_resist.Rds"))

    subsetSwitchGB <- isoformSwitchTestDEXSeq(
        ### Core arguments
        switchAnalyzeRlist= subsetSwitchGB,
        alpha = 0.05,
        dIFcutoff = 0.1,
        quiet=TRUE
    )
    saveRDS(subsetSwitchGB, file = paste0(cancer,"_DEXSeq_lncrna_resist.Rds"))
}
