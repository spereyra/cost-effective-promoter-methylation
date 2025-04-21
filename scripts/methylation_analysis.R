# Main Script: `methylation_analysis.R`
# Title: "Cost-Effective Promoter Methylation Analysis via Target Long-Read Bisulfite Sequencing"
# Authors: Silvana Pereyra1*, Angela Sardina1, Rita Neumann2, Celia May2, Rossana Sapiro3, Bernardo Bertoni1, Mónica Cappetta
# Institutions:1 Unidad Académica de Genética, Facultad de Medicina, Universidad de la República, Montevideo, Uruguay. 2 Department of Genetics & Genome Biology, College of Life Sciences, University of Leicester, UK. 3 Unidad Académica de Histología y Embriología, Facultad de Medicina, Universidad de la República. Montevideo, Uruguay.
# Date: April 2025
# Purpose: This script performs the analysis of DNA methylation data with a targeted sequencing approach using Oxford Nanopore. It includes filtering of CpG sites, statistical tests, and visualization of results in boxplots and RNA-seq expression comparisons.

# Load required libraries
library(readr)
library(dplyr)
library(GenomicRanges)

# Load sample design metadata (replace with actual paths)
design_files <- list(
  Run1 = "Documents/ProjectFolder/Run1.csv",
  Run2 = "Documents/ProjectFolder/Run2.csv"
)

design_data <- lapply(design_files, read.csv)

# Define paths and load all count files
data_paths <- list(
  Run1 = "Documents/ProjectFolder/ConteosRun1",
  Run2 = "Documents/ProjectFolder/ConteosRun2"
)

# Define filename patterns for each run
patterns <- list(
  Run1 = "trimmed-barcode", # Adjust per your filename, if neccesary
  Run2 = "trimmed-barcode" # Adjust per your filename, if neccesary
)

# Helper function to read and format raw count files
read_counts <- function(path, pattern) {
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  df_list <- lapply(files, function(f) {
    read.delim(f, col.names = c("Chr", "start", "end", "Perc", "CountMet", "CountUnMet"))
  })
  df <- bind_rows(df_list, .id = "barcode")
  df$CountMet <- as.numeric(df$CountMet)
  df$CountUnMet <- as.numeric(df$CountUnMet)
  return(df)
}

# Read and process all datasets
df_list <- Map(read_counts, data_paths, patterns)
names(df_list) <- names(data_paths)

# Add sample names to each dataframe using corresponding design file
for (run_name in names(df_list)) {
  df <- df_list[[run_name]]
  design <- design_data[[run_name]]
  df$Muestra <- NA
  for (i in 1:12) {
    df$Muestra[df$barcode == i] <- rep(
      as.character(design[design$Barcode == i, "Ind_PPT"]),
      sum(df$barcode == i)
    )
  }
  df$Source <- run_name
  df_list[[run_name]] <- df
}

# Combine all runs into one dataframe
df_all <- bind_rows(df_list)

# Filter complete cases, collapse duplicates, and calculate methylation percentage
df_clean <- df_all %>%
  filter(complete.cases(.)) %>%
  group_by(Muestra, Chr, start) %>%
  summarise(across(c(CountMet, CountUnMet), sum), .groups = "drop") %>%
  mutate(
    end = start,
    Perc = 100 * CountMet / (CountMet + CountUnMet)
  ) %>%
  relocate(end, .after = start) %>%
  relocate(Perc, .before = CountMet)

# Keep only rows with minimum read coverage (CountMet + CountUnMet >= 30)
df_filtered <- df_clean[rowSums(df_clean[, c("CountMet", "CountUnMet")]) > 29, ]

# Convert to GRanges for genomic analysis
gr <- makeGRangesFromDataFrame(
  df_filtered,
  seqnames.field = "Chr",
  start.field = "start",
  end.field = "end",
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = TRUE
)
genome(gr) <- "hg38"

# Load BED file with regions of interest
genes_bed <- read.csv(
  "Documents/ProjectFolder/GenesBED_PPT.bed.csv",
  header = TRUE,
  col.names = c("chrom", "chromStart", "chromEnd", "name", "strand", "geneName")
)

# Ready to intersect GRanges with genes of interest or proceed with analysis
# Convert data frame to GRanges
require(rtracklayer)
GenesSecs.gr <- makeGRangesFromDataFrame(GenesSecs,
                                         seqnames.field = c("X.chrom"),
                                         start.field = c("chromStart"),
                                         end.field = c("chromEnd"),
                                         keep.extra.columns = TRUE,
                                         starts.in.df.are.0based = TRUE)
genome(GenesSecs.gr) <- "hg38"

# Subset overlapping regions between CpG data and candidate genes
df.gr_PPT <- subsetByOverlaps(df.gr, GenesSecs.gr)

# Generate list of candidate genes
genes_PPT <- unique(GenesSecs$geneName)

# Selected indexes
genes_idx <- c(1:7, 9:11, 13:14)
genes_sel <- genes_PPT[genes_idx]

# Create a list of GRanges objects, one per gene
df_LIST <- setNames(lapply(genes_sel, function(gen) {
  subsetByOverlaps(df.gr, GenesSecs.gr[GenesSecs.gr$geneName == gen, ])
}), genes_sel)

# Count the number of CpGs per gene
NoCpG <- setNames(data.frame(lapply(df_LIST, function(x) length(unique(x)))), 
                  genes_sel)

# Compute average sequencing coverage per gene
resultados <- do.call(rbind, lapply(seq_along(df_LIST), function(i) {
  cobertura_gen <- as.data.frame(df_LIST[[i]]) %>%
    mutate(Coverage = rowSums(.[8:9])) %>%
    summarize(Mean = mean(Coverage, na.rm = TRUE), 
              Max = max(Perc, na.rm = TRUE), 
              Median = median(Perc, na.rm = TRUE), 
              SD = sd(Perc, na.rm = TRUE))
  cobertura_gen$Gen <- names(df_LIST)[i]
  return(cobertura_gen)
}))

# Add gene names and export coverage table
write.csv(resultados %>% select(Mean, SD, Gene), 
          file = "CoberturaXgen.csv", 
          row.names = FALSE)

# Wilcoxon test: example for one gene (PIK3AP1)
# Note: This example shows the procedure for the first gene only.
# The same approach should be applied iteratively to each gene of interest.

library(Repitools)
DFrames <- lapply(df_LIST, function(x) annoGR2DF(x))
DFrames.wide <- list()
Wilcoxon <- list()

# Prepare data for gene 1 (PIK3AP1)
DFrames[[1]]$id <- paste(DFrames[[1]]$chr, ":", DFrames[[1]]$start, sep = "")
DFrames[[1]] <- DFrames[[1]][, -c(1:4, 7:8)]  # Keep only 'id', 'Muestra', 'Perc'
DFrames.wide[[1]] <- reshape(DFrames[[1]], idvar = "id", timevar = "Muestra", direction = "wide")

# Perform Wilcoxon test per CpG (row-wise)
Wilcoxon[[1]] <- apply(DFrames.wide[[1]][complete.cases(DFrames.wide[[1]]), ], 1, function(x) {
  control <- as.numeric(x[c(2:3)])   # Adjust indexes based on control samples
  patient <- as.numeric(x[c(4:5)])   # Adjust indexes based on patient samples
  wilcox.test(patient, control, na.action = "na.exclude")$p.value
})

# Find CpGs with significant differences
which(Wilcoxon[[1]] < 0.051)

# NOTE: Repeat the above steps for each gene in the list to perform a full Wilcoxon analysis.


# HEATMAP for MIR155HG ------------------------------------

library(plyr)
library(viridis)
library(heatmaply)
library(RColorBrewer)

# Extract data for MIR155HG
df_HeatMap_MIR155 <- bind_rows(DFrames[["MIR155HG"]])
df_HeatMap_MIR155$id <- paste0(df_HeatMap_MIR155$chr, ":", df_HeatMap_MIR155$start)

# Remove unnecessary columns
df_HeatMap_MIR155 <- df_HeatMap_MIR155[, !(colnames(df_HeatMap_MIR155) %in% c("chr", "start", "end", "strand", "no", "total"))]

# Reshape the data to wide format (CpGs as rows, samples as columns)
DataHeatMIR <- reshape(df_HeatMap_MIR155, idvar = "id", timevar = "Muestra", direction = "wide")
colnames(DataHeatMIR) <- sub("Perc.", "", colnames(DataHeatMIR))
rownames(DataHeatMIR) <- DataHeatMIR$id
DataHeatMIR <- DataHeatMIR[, -1]
DataHeatMIR <- data.matrix(DataHeatMIR)
x <- DataHeatMIR[complete.cases(DataHeatMIR), ]

# Labels for columns and rows
labColMIR <- Muestras$Codigos[match(colnames(x), Muestras$ID)]
labrowMIR <- graficarMIR155HG.long$relPos[match(rownames(x), graficarMIR155HG.long$id)]

# Colors for sample conditions
MuestrasMIR <- Muestras[match(colnames(DataHeatMIR), Muestras$ID), ]
cond_colors <- ifelse(MuestrasMIR$Condicion == "control", "#F8766D", "#66BBBB")

# Color palette
my_palette <- colorRampPalette(brewer.pal(9, "YlGnBu"))(101)

# Plot
png(filename = "Figure2.HeatmapMIR.png", width = 800, height = 1000)
heatmap.2(x,
          scale = "none",
          col = my_palette,
          labCol = FALSE,
          labRow = paste0("+", labrowMIR),
          cexRow = 0.8,
          breaks = seq(0, 100, 1),
          dendrogram = 'column',
          Rowv = FALSE,
          trace = "none",
          density.info = "none",
          key.xlab = "Methylation percentage",
          ColSideColors = cond_colors,
          key.title = "")
dev.off()
