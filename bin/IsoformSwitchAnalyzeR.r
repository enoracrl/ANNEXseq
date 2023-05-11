                                            #########################
                                            # IsoformSwitchAnalyzeR #
                                            #########################


##########################################################################################################


# Setup R Environment

### Load library
suppressMessages(library(dplyr))
suppressMessages(library(IsoformSwitchAnalyzeR))
suppressMessages(library(tidyr))
suppressMessages(library(data.table))
suppressMessages(library(ggvenn))
suppressMessages(library(cowplot))
suppressMessages(library(UpSetR))
suppressMessages(library(GenomicFeatures))

##########################################################################################################

args = commandArgs(trailingOnly=TRUE)

design          <- strsplit(grep('--design*', args, value = TRUE), split = '=')[[1]][[2]]           # in input
counts          <- strsplit(grep('--count_matrix*', args, value = TRUE), split = '=')[[1]][[2]]     
exon_annot      <- strsplit(grep('--annotation*', args, value = TRUE), split = '=')[[1]][[2]]
transcriptome   <- strsplit(grep('--transcriptome*', args, value = TRUE), split = '=')[[1]][[2]]
comparisons     <- strsplit(grep('--comparisons*', args, value = TRUE), split = '=')[[1]][[2]]      # in input

# Experimental Data
design      <- as.data.frame.matrix(read.table(design, header=TRUE, sep=','))
types       <- unique(design$condition)
design      <- unite(design, condition, condition, condition, sep = "_")
design      <- subset(design, select = c(sample_id, condition))
names(design)[names(design) == "sample_id"] <- "sampleID"

# Header : samples names to match between the design file and the counts file
samples         <- (c(design[,1]))
samples         <- sort(samples)
samples[1:2]    <- rev(samples[1:2])

# Isoforms counts
counts              <- as.data.frame.matrix(read.table(counts, header=TRUE))
counts              <- subset(counts, select = -c(GENEID))
colnames(counts)    <- c("isoform_id",samples)

# Comparisons table
comparisons <- data.frame(read.table(comparisons, header=TRUE))


txdb <- makeTxDbFromGFF(exon_annot, format="gtf")
tx <- data.frame(transcripts(txdb))
names(tx)[names(tx) == 'tx_name'] <- 'isoform_id'
tx <- tx[,-6]



# ##################### #
#  STEP 1 - SOME STATS  #
# ##################### #

samples_mean=c()
for (i in seq(2,25)){
    samples_mean=c(samples_mean, mean(counts_matrice[,i]))
}
df = data.frame(samples,samples_mean)

# Mean number of isoform per sample
ggplot(df,aes(x=samples,y=samples_mean)) +
    geom_col(aes(fill=samples_mean)) +
    coord_flip() +
    geom_text(data=df, aes(x=samples,y=samples_mean/2,label=round(samples_mean,2)),vjust=0,angle=0,color="white") +
    theme(axis.text.x = element_text(angle=90),
    legend.position = "none")

samples_sum=c()
for (i in seq(2,25)){
    samples_sum=c(samples_sum, sum(counts_matrice[, i]))
}
df3 = data.frame(samples, samples_sum)

# Sum of isoform number per sample
ggplot(df3 , aes(x=samples , y=samples_sum)) +
    geom_col(aes(fill=samples_sum)) +
    coord_flip() +
    geom_text(data=df, aes(x=samples, y=samples_sum/2, label=round(samples_sum, 2)), vjust=0, angle=0, color="white") +
    theme(axis.text.x = element_text(angle=90),
    legend.position = "none")

iso_sup=c()
for (i in seq(2,25)){
    iso_sup=c(iso_sup, nrow(counts_matrice[counts_matrice[,i] == 0,])/nrow(counts_matrice)*100)
}
df2 = data.frame(samples, iso_sup)

# Number of isoform that are not present per sample
ggplot(df2, aes(x=samples, y=iso_sup)) +
    geom_col(aes(fill=iso_sup)) +
    coord_flip() +
    geom_text(data=df2, aes(x=samples, y=iso_sup/2, label=round(iso_sup, 2)), vjust=0,  angle=0, color="white") +
    theme(axis.text.x = element_text(angle=90),
    legend.position = "none")




# ##################### #
# STEP 2 - LOADING DATA #
# ##################### #

message('Importing data')
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
saveRDS(SwitchList, file = "SwitchList.Rds")

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


counts      = c()
countsDE    = c()
downRegTx   = list()
upRegTx     = list()
# Loop 

for (condition in groups) {

    message(paste0(condition, 'subseting'))
    subsetSwitch <- subsetSwitchAnalyzeRlist(
        switchAnalyzeRlist=SwitchList,
        SwitchList$isoformFeatures$condition_1 %in% c(
            paste0(condition, "_sensitive"),
            paste0(condition, "_resistant")
        ) &
        SwitchList$isoformFeatures$condition_2 %in% c(
            paste0(condition, "_sensitive"),
            paste0(condition, "_resistant")
        )
    )

    message(paste0(condition, 'filtering'))
    subsetSwitch <- preFilter(
        switchAnalyzeRlist= subsetSwitch,
        geneExpressionCutoff = 1,
        isoformExpressionCutoff = 0,
        quiet=TRUE
    )
    saveRDS(subsetSwitch, file = paste0(condition,"_lncrna_resist.Rds"))

    
    subsetSwitchDE <- isoformSwitchTestDEXSeq(
        ### Core arguments
        switchAnalyzeRlist= subsetSwitch,
        alpha = 0.05,
        dIFcutoff = 0.1,
        quiet=TRUE
    )

    subsetSwitchDE <- subsetSwitchDE$isoformFeatures[abs(
        subsetSwitchDE$isoformFeatures$dIF) > 0.1  & 
        subsetSwitchDE$isoformFeatures$isoform_switch_q_value < 0.05, ]
    subsetSwitchDE <- subsetSwitchDE %>% mutate(
        isoform_annot = case_when(
            startsWith(isoform_id, "BambuTx") ~ "NOVEL",
            startsWith(isoform_id, "ENST") ~ "KNOWN"
            )) %>%
        mutate(gene_annot = case_when(
            startsWith(gene_id, "BambuGene") ~ "NOVEL",
            startsWith(gene_id, "ENSG") ~ "KNOWN"
        ))
    subsetSwitchDE <- merge(subsetSwitchDE, tx, by = "isoform_id")
    saveRDS(subsetSwitch, file = paste0(condition,"_DEXSeq_lncrna_resist.Rds"))

    counts = c(
        counts, 
        type = c(
            length(unique(subsetSwitch$isoformFeatures$isoform_id[grep("^[ENST].*", subsetSwitch$isoformFeatures$isoform_id)])),
            length(unique(subsetSwitch$isoformFeatures$isoform_id[grep("^[BambuTx].*", subsetSwitch$isoformFeatures$isoform_id)])),
            length(unique(subsetSwitch$isoformFeatures$gene_id[grep("^[ENSG].*", subsetSwitch$isoformFeatures$gene_id)])),
            length(unique(subsetSwitch$isoformFeatures$gene_id[grep("^[BambuGene].*", subsetSwitch$isoformFeatures$gene_id)]))
            )
        )
    countsDE = c(
        countsDE, 
        type = c(
            length(unique(subsetSwitchDE$isoform_id[grep("^[ENST].*", subsetSwitchDE$isoform_id)])),
            length(unique(subsetSwitchDE$isoform_id[grep("^[BambuTx].*", subsetSwitchDE$isoform_id)])),
            length(unique(subsetSwitchDE$gene_id[grep("^[ENSG].*", subsetSwitchDE$gene_id)])),
            length(unique(subsetSwitchDE$gene_id[grep("^[BambuGene].*", subsetSwitchDE$gene_id)]))
        ))

    downRegTx   <- setNames(lapply(groups, function(condition) { subsetSwitchDE$isoform_id[subsetGB_DE$dIF < 0.1]}), groups)
    upRegTx     <- setNames(lapply(groups, function(condition) { subsetSwitchDE$isoform_id[subsetGB_DE$dIF > 0.1]}), groups)


}

# Plots
palette = c('#85AE5B', '#5C9EC1', '#F4B266', '#e7c6ff', '#A7333F','#B6A6CA', '#adc178','#ade8f4','#ffc2d1','#EB9784','#C44536','#f79d65','#F6D379','#F2B79F','#006633')
pdf(paste("plots", ".pdf", sep = ""), width = 20, height = 12)



# Histogram of the number of novel and known genes/isoforms per type
groupsDF <- data.frame(
    group   = rep(groups, each = length(groups)), 
    type    = rep(c("isoforms","genes"), each =2, times=length(groups)),
    prop    = rep(c("known", "novel"), times=2*length(groups)),
    count   = counts
)

groupsDF_DE <- data.frame(
    group   = rep(groups, each = length(groups)),
    type    = rep(c("isoforms", "genes"), each=2, times=length(groups)),
    prop    = rep(c("known", "novel"), times=2*length(groups)),
    count   = countsDE
)

p1 <- ggplot(groupsDF, aes(fill=type, y=count, x=type, group=type)) +
    geom_bar(stat="identity", position="stack", color="black", aes(alpha=prop)) +
    facet_wrap( ~ group, nrow=1) +
    labs(x="Sample type", y="Counts", subtitle = "After filtering") +
    scale_alpha_manual(values=c(0.9, 0.3)) +
    geom_text(aes(label = count, y = ifelse(count > 1500, count, count + 3500)), position = position_stack(0.5), angle = 0, size= 4, color = "black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.subtitle = element_text(hjust = 0.5, size = 15), text = element_text(size = 13)) +
    scale_fill_manual(values=c("#99E2FF", "#93C8A9"))

p2 <- ggplot(groupsDF_DE, aes(fill=type, y=count, x=type, group=type)) +
    geom_bar(stat="identity", position="stack", color="black", aes(alpha=prop)) +
    facet_wrap( ~ group, nrow=1) +
    labs(x="Sample type", y="Counts", subtitle = "After DEXSeq analysis (alpha < 0.05, abs(dIF) > 0.1)") +
    scale_alpha_manual(values=c(0.9, 0.3)) +
    geom_text( aes(label = count, y = ifelse(count > 30, count, count + 50)), position = position_stack(0.5), angle = 0, size= 4, color = "black") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.subtitle = element_text(hjust = 0.5, size = 15), text = element_text(size = 13)) +
    scale_fill_manual(values=c("#99E2FF", "#93C8A9"))

plots <- list(p1, p2)
prow <- plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), align = 'vh', hjust = -1, nrow = 1)
legend <- get_legend(p1 + theme(legend.box.margin = margin(0, 0, 0, 12)))
prow2 <- plot_grid(prow, legend, rel_widths = c(3, .4))
title <- ggdraw() + draw_label("Number of known and novel isoforms and genes per group", fontface='bold', size=20)
pdf(paste("counts_group", ".pdf", sep = ""), width = 12, height = 7)
plot_grid(title, prow2, nrow=2, rel_heights=c(0.1, 0.7,0.4,0.1))
dev.off()

# Down-regulated genes
pdf(paste("upsetplot_downDTU", ".pdf", sep = ""), width = 12, height = 7)
UpSetR::upset(
    fromList(downRegTx), 
    keep.order=TRUE, 
    sets.bar.color=palette[1:length(groups)],
    point.size = 4, line.size = 1, text.scale=2
    )
dev.off()

# Up-regulated genes
pdf(paste("upsetplot_upDTU", ".pdf", sep = ""), width = 12, height = 7)
UpSetR::upset(
    fromList(upRegTx), 
    keep.order=TRUE, 
    sets.bar.color=palette[1:length(groups)],
    point.size = 4, line.size = 1, text.scale = 2
    )
dev.off()