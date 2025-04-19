# Load required libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(stringr)

# Initialize objects
subsets <- list()
Muestra_filtro <- list()
rango <- list()
paratest.control <- list()
paratest.ppt <- list()
resultados <- data.frame()
shapiro.c <- list()
shapiro.p <- list()
testT <- list()
maximo <- list()
minimo <- list()
df <- list()
testWilcoxon <- list()

# Rename first column of Muestras to match
colnames(Muestras)[1] <- "Muestra"

# Loop over each gene (assumed to be 12 total)
for (i in 1:12) {
  
  # Get the sample with the fewest CpGs sequenced for this gene
  Muestra_filtro[[i]] <- DFrames[[i]] %>% 
    group_by(Muestra) %>%
    tally() %>% 
    slice_min(n) %>%
    select(Muestra)
  
  # Save min and max number of CpGs per sample for this gene
  minimo[[i]] <- DFrames[[i]] %>%
    group_by(Muestra) %>%
    tally() %>%
    slice_min(n) %>%
    pull(n)
  
  maximo[[i]] <- DFrames[[i]] %>%
    group_by(Muestra) %>%
    tally() %>%
    slice_max(n) %>%
    pull(n)
  
  # Select CpG IDs from the sample with the fewest sequenced CpGs
  subsets[[i]] <- DFrames[[i]] %>%
    filter(Muestra == Muestra_filtro[[i]]) %>%
    distinct(id)
  
  # Extract methylation percentages for preterm samples
  paratest.ppt[[i]] <- DFrames[[i]] %>%
    filter(id %in% subsets[[i]]$id) %>%
    left_join(Muestras, by = "Muestra") %>%
    filter(Condicion == "Prematuro") %>%
    pull(Perc)
  
  # Extract methylation percentages for control samples
  paratest.control[[i]] <- DFrames[[i]] %>%
    filter(id %in% subsets[[i]]$id) %>%
    left_join(Muestras, by = "Muestra") %>%
    filter(Condicion == "control") %>%
    pull(Perc)
  
  # Shapiro-Wilk normality test for both groups
  shapiro.c[[i]] <- shapiro.test(paratest.control[[i]])$p.value
  shapiro.p[[i]] <- shapiro.test(paratest.ppt[[i]])$p.value
  
  # T-test (note: not used if data is non-normal)
  testT[[i]] <- t.test(paratest.control[[i]], paratest.ppt[[i]])$p.value
  
  # Wilcoxon test (non-parametric)
  testWilcoxon[[i]] <- wilcox.test(x = paratest.control[[i]], y = paratest.ppt[[i]], 
                                   alternative = "two.sided", mu = 0,
                                   paired = FALSE, correct = FALSE, conf.int = 0.95)$p.value
  
  # Collect all results in a single row for this gene
  resultados <- rbind.data.frame(
    resultados,
    c(
      names(DFrames)[i],
      minimo[[i]],
      maximo[[i]],
      mean(DFrames[[i]]$Perc, na.rm = TRUE),
      median(DFrames[[i]]$Perc, na.rm = TRUE),
      shapiro.c[[i]],
      shapiro.p[[i]],
      testT[[i]],
      testWilcoxon[[i]]
    ),
    stringsAsFactors = FALSE
  )
  
  # Prepare data for boxplot
  df[[i]] <- melt(
    DFrames[[i]] %>%
      filter(id %in% subsets[[i]]$id) %>%
      left_join(Muestras, by = "Muestra") %>%
      select(Perc, Condicion),
    id = "Condicion"
  )
}

# Assign column names to the result table
colnames(resultados) <- c(
  "GeneName", "MinCpG", "MaxCpG", "Mean", "Median", 
  "Shapiro_Control", "Shapiro_Premature", "TTest", "WilcoxonTest"
)

# Adjust p-values using FDR correction
resultados$FDR <- p.adjust(resultados$WilcoxonTest, method = "fdr")

# Label lists with gene names
names(paratest.control) <- names(DFrames)
names(paratest.ppt) <- names(DFrames)
names(df) <- names(DFrames)

# Log2FC table from publication Pereyra et al 2019
log2FC_table <- data.frame(
  Gene = c("ACCS", "ANKRD24", "BIRC3", "CXCL2", "EGR3", "GK", "MAMDC2", "MIR155HG", "NRN1", "PIK3AP1", "STEAP1", "WNT1"),
  LogFC = c(-2.79, -2.21, 3.36, 4.96, 3.67, 5.04, -2.23, 2.66, -2.19, 4.37, 4.17, -3.29)
)

# Filter data for significantly differentially methylated genes
datosGruped <- melt(df) %>%
  filter(L1 %in% resultados %>% filter(WilcoxonTest < 0.05) %>% pull(GeneName))

# Merge with Log2FC data for annotation
merged_dataGrouped <- merge(datosGruped, log2FC_table, by.x = "L1", by.y = "Gene", all.x = TRUE)
merged_dataGrouped$facet_title <- paste0(merged_dataGrouped$L1, "\n", "Log2FC: ", round(merged_dataGrouped$LogFC, 2))

# Boxplot for methylation percentages
grupedCpG <- ggplot(merged_dataGrouped, aes(x = variable, y = value, fill = Condicion)) + 
  geom_boxplot() + 
  facet_grid(~ L1, scales = "free") +
  ylab("Methylation percentage") +
  xlab("Grouped CpGs") +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold"), 
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.key.size = unit(0.8, "lines"),
    text = element_text(size = 14), 
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.spacing = unit(0.2, "lines"),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "grey95", size = 1),
    strip.text = element_text(face = "italic"),
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    plot.background = element_blank(), 
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 10)
  ) +
  scale_y_continuous(limits = c(0, 101), expand = c(0, 0)) +
  scale_fill_discrete(labels = c('Term control', 'Severe preterm')) +
  guides(fill = guide_legend(title = "Condition"))

# RNA-seq boxplot (dds_norm_long_filtered is a data.frame obtained from Deseq analysis with Pereyra et al 2019 results)
p2 <- ggplot(dds_norm_long_filtered, aes(x = Condicion, y = count, fill = Condicion)) + 
  geom_boxplot() +
  facet_wrap(~ facet_title, scales = "free_y", nrow = 1) +
  theme_bw() +
  ylab("Transcriptome \n counts") +
  xlab("Condition") +
  scale_fill_manual(values = c("control" = "#F8766D", "Prematuro" = "#66BBBB")) +
  scale_x_discrete(labels = c("control" = "Term \n controls", "Prematuro" = "Severe \n preterms")) +
  guides(fill = guide_legend(title = "Condition")) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.8, vjust = 0.5),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.key.size = unit(0.8, "lines"),
    legend.position = "none",
    text = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.spacing = unit(0.2, "lines"),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(color = "black", fill = "grey95", size = 1),
    strip.text = element_text(face = "italic"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Combine both plots into a single figure
combinadotodo <- ggarrange(
  grupedCpG, p2,
  ncol = 1,
  nrow = 2,
  heights = c(0.6, 0.4),
  labels = c("A", "B")
)

# Save output image
ggsave("Combined_GroupedCpG_RNAseq.png", plot = combinadotodo, width = 9, height = 7, units = "in", dpi = 300)