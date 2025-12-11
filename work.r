# Set working environment
rm(list = ls())
gc()
setwd('/home/workdir/project_for_test/HT20220905701/')
set.seed(666)

# Load required R packages
library(tidyverse)
library(data.table)
library(ggplot2)
library(PCAtools)
library(colorspace)
library(pheatmap)

# Read input data
data_path <- '01.raw_data/'
df_pca_res <- fread(paste0(data_path, 'tpm_group_pca.txt'))
df_tpm_all <- fread(paste0(data_path, 'gene_tpm_all_samples_anno.xls'))
df_intestine <- fread(paste0(data_path, 'P_I_vs_N_I_all_anno.xls'))
df_kidney <- fread(paste0(data_path, 'P_K_vs_N_K_all_anno.xls'))
df_spleen <- fread(paste0(data_path, 'P_S_vs_N_S_all_anno.xls'))
df_liver <- fread(paste0(data_path, 'P_L_vs_N_L_all_anno.xls'))
df_sample_cor <- read.table(paste0(data_path, 'tpm_corr.txt'), row.names = 1)

######################### Step 1: PCA Analysis #####################
# Data preprocessing for PCA
df_pca <- column_to_rownames(df_tpm_all |> dplyr::select(-2), var = 'gene_id')

# Construct metadata for PCA
df_pca_metadata <- tibble(
  sample_id = colnames(df_pca),
  Treatment = colnames(df_pca) |> str_sub(start = 1, end = 1),
  organ = colnames(df_pca) |> str_sub(start = 3, end = 3)
)

df_pca_metadata <- column_to_rownames(df_pca_metadata, var = 'sample_id')

# Assign human-readable organ names
df_pca_metadata <- df_pca_metadata |>
  mutate(organ_new = if_else(
    organ == "I", "Intestine",
    if_else(organ == "K", "Kidney",
            if_else(organ == "L", "Liver",
                    if_else(organ == "S", "Spleen", "Unknown")))))

# Assign human-readable treatment labels
df_pca_metadata <- df_pca_metadata |> 
  mutate(Treatment_new = if_else(Treatment == 'N', 'Normal', 'Parasite'))

df_pca_metadata <- df_pca_metadata[, c(3,4)]
names(df_pca_metadata) <- c('Organ', 'Treatment')

# Perform PCA
pca_res <- pca(df_pca, metadata = df_pca_metadata, removeVar = 0.1)

# PCA biplot colored by organ type
biplot(pca_res, 
       colby = 'Organ', 
       legendPosition = 'right',
       encircle = TRUE,
       max.overlaps = 40,
       encircleLineCol = 'black', 
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       borderWidth = 0.5,
       borderColour = 'black')

# Extract PCA rotated coordinates
pca_df <- data.frame(
  PC1 = pca_res$rotated[,1],
  PC2 = pca_res$rotated[,2],
  Organ = pca_res$metadata$Organ,
  Sample = rownames(pca_res$rotated)
)

library(ggforce)

# Generate PCA plot using ggplot2
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Organ, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_mark_ellipse(aes(fill = Organ, label = NULL), alpha = 0.2, color = NA) +
  geom_text_repel(show.legend = FALSE, max.overlaps = 40, size = 3) +
  theme_bw(base_size = 14) +
  labs(
    x = paste0("PC1, ", round(pca_res$variance[1], 2), "% variation"),
    y = paste0("PC2, ", round(pca_res$variance[2], 2), "% variation"),
    title = ""
  ) +
  theme(legend.position = "right") +
  coord_cartesian(
    xlim = c(min(pca_df$PC1) - 5000, max(pca_df$PC1) + 23000),
    ylim = c(min(pca_df$PC2) - 5000, max(pca_df$PC2) + 10000),
    clip = "off"
)

# Alternative PCA plot with legend inside the plot area
ggplot(pca_df, aes(x = PC1, y = PC2, color = Organ, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_mark_ellipse(aes(fill = Organ, label = NULL), 
                    alpha = 0.2, color = NA, show.legend = FALSE) +
  geom_text_repel(show.legend = FALSE, max.overlaps = 40, size = 3) +
  labs(
    x = paste0("PC1, ", round(pca_res$variance[1], 2), "% variation"),
    y = paste0("PC2, ", round(pca_res$variance[2], 2), "% variation"),
    title = NULL
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(
    xlim = c(min(pca_df$PC1) - 5000, max(pca_df$PC1) + 23000),
    ylim = c(min(pca_df$PC2) - 5000, max(pca_df$PC2) + 10000),
    clip = "off"
  ) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "gray70"),
    legend.key = element_rect(fill = NA, color = NA),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold")
  )

# Display PCA plot
print(p)

# PCA categorized by treatment condition
biplot(pca_res, 
       colby = 'Treatment',
       legendPosition = 'right',
       encircle = FALSE,
       max.overlaps = 40,
       encircleLineCol = 'black', 
       gridlines.major = FALSE,
       gridlines.minor = FALSE,
       borderWidth = 0.5,
       borderColour = 'black')
######################### Step 2: Comparison of Average TPM Expression #####################

df_tpm_all <- df_tpm_all[,-2]   # Remove the 2nd column (annotation column)

# Convert to long format
df_long <- df_tpm_all %>%
  pivot_longer(cols = 2:25, names_to = "SampleID", values_to = "TPM")

# Extract GroupID information (e.g., N_I, P_L, etc.)
df_long <- df_long %>%
  mutate(GroupID = str_extract(SampleID, "^[A-Z]_[A-Z]"))

# Log10-transform TPM values
df_long <- df_long %>%
  mutate(logTPM = log10(TPM + 0.01))  # Add 0.01 to avoid log(0)

# Reorder group names for visualization
df_long$GroupID_new <- df_long$GroupID
df_long$GroupID_new <- df_long$GroupID_new |>
  str_replace_all("N_I", "I_N") |>
  str_replace_all("N_K", "K_N") |>
  str_replace_all("N_L", "L_N") |>
  str_replace_all("N_S", "S_N") |>
  str_replace_all("P_I", "I_P") |>
  str_replace_all("P_K", "K_P") |>
  str_replace_all("P_L", "L_P") |>
  str_replace_all("P_S", "S_P")

df_long$GroupID_new <- factor(
  df_long$GroupID_new, 
  levels = c("I_N","I_P","K_N","K_P","S_N","S_P","L_N","L_P")
)

# Draw boxplot of expression level distribution
ggplot(df_long, aes(x = GroupID_new, y = logTPM, fill = GroupID_new)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  labs(x = "GroupID", y = "log10(TPM)", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  guides(fill = guide_legend(title = "Group"))

######################### Step 3: Sample Correlation #####################

# Display available HCL palettes
hcl_palettes(plot = TRUE)

# Generate sequential HCL palette
sequential_hcl(10, palette = 'Heat2')

# Draw sample correlation heatmap
pheatmap(
  df_sample_cor,
  color = sequential_hcl(10, palette = 'Heat2') |> rev()
)
######################### Step 4: Bar Plot of Differential Gene Counts #####################
# Differential genes are defined as those with abs(logFC) > 1 and FDR < 0.05
list_degs <- list(df_intestine, df_kidney, df_spleen, df_liver)
names(list_degs) <- c('Intestine', 'Kidney', 'Spleen', 'Liver')

list_degs_1 <- lapply(list_degs, function(x){
  x <- x |> dplyr::select(gene_id, logFC, pvalue, FDR, level)
  x$level <- str_replace_all(x$level, 
                             pattern = 'Decreased',
                             replacement = 'Down')
  x$level <- str_replace_all(x$level, 
                             pattern = 'Increased',
                             replacement = 'Up')
  x$level <- str_replace_all(x$level, 
                             pattern = 'nonsignificant',
                             replacement = 'Nonsig')
  return(x)
})

# Read pre-calculated DEG counts data for plotting
df_plot_degs <- fread('degs_number_plot.txt')
names(df_plot_degs)
df_plot_degs$Group <- factor(df_plot_degs$Group, levels = c('Intestine','Liver','Spleen','Kidney'))
df_plot_degs$Level <- factor(df_plot_degs$Level, levels = c('Down', 'Up','Total' ))

# Plot bar chart of differential gene numbers
ggplot(df_plot_degs, aes(x = Group, y = Count, fill = Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.6) +
  geom_text(aes(label = Count), 
            position = position_dodge2(width = 0.6, preserve = "single"),
            hjust = -0.2,
            size = 3,
            color = "black") +
  ylim(0,3100) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Up" = "#F38181",     # light coral red
    "Down" = "#95E1D3",   # soft green
    "Total" = "#A8D8B9"   # light cyan-blue
  )) +
  labs(title = "",
       x = "Organ", 
       y = "Gene Count") +
  theme_minimal() +
  theme(
    legend.position = c(0.95, 0.05),              # Legend in the bottom right inside the plot
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "gray70"),
    legend.key = element_rect(fill = NA, color = NA),
    legend.title = element_blank(),
    axis.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  )

######################### Step 4 (continued): Volcano / MA Plot #####################
library(extrafont)

plot_ma <- function(df){
  df <- list_degs[[1]] |> dplyr::select(-58)
  
  df <- df |> 
    dplyr::distinct(gene_id,.keep_all = TRUE) |> 
    dplyr::select('gene_id', 'gene_name', 'logFC', 'logCPM', 'FDR', 'level')
  
  df$level <- df$level |>
    str_replace_all('Decreased', 'Down') |>
    str_replace_all('Increased', 'Up') |>
    str_replace_all('nonsignificant', 'Nonsig') 
  
  # Extract top 5 up-regulated and top 5 down-regulated genes (FDR < 0.05)
  dat_rep <- df %>%
    filter(FDR < 0.05) %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 5) %>%
    bind_rows(
      df %>%
        filter(FDR < 0.05) %>%
        arrange(logFC) %>%
        slice_head(n = 5)
    )
  
  # Color mapping for MA plot
  color_map <- c("Up" = "#7A4D7B",   # purple tone for up-regulated
                 "Down" = "#008080", # teal for down-regulated
                 "Nonsig" = "#808080") # grey for non-significant
  
  # Draw MA plot
  ggplot(df, aes(x = logCPM, y = logFC, color = level)) +
    geom_point(aes(size = abs(logFC)), alpha = 0.4, na.rm = TRUE) +
    scale_size_continuous(range = c(1.5, 6), guide = "none") +
    scale_alpha_continuous(range = c(0.3, 0.8), guide = "none") +
    scale_color_manual(
      values = color_map,
      breaks = c("Up", "Down", "Nonsig")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    labs(
      x = "logCPM",
      y = "log2 Fold Change",
      title = ""
    ) +
    theme_minimal(base_size = 13) +
    geom_label_repel(
      data = dat_rep,
      aes(label = gene_id),
      max.overlaps = 20,
      size = 4,
      box.padding = unit(0.5, "lines"),
      min.segment.length = 0,
      point.padding = unit(0.8, "lines"),
      segment.color = "black",
      show.legend = FALSE
    ) +
    theme(
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.title = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}

# Draw volcano plot
ggplot(data = df,
       aes(x= logFC, y = -log10(FDR), color = level)) +
  scale_color_manual(values = c('#008080','#7A4D7B', '#808080'))+ # custom color scheme
  geom_point(size = 2.4, alpha = 0.4, na.rm = TRUE) +
  theme_bw(base_size = 12) +
  geom_vline(xintercept = c(-1, 1),
             linetype = 4,
             color = "gray40",
             size = 0.6) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 4,
             color = "gray40",
             size = 0.6) +
  ylim(0,80) +
  theme_clean() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold", color = "black", family = "Times New Roman", size = 13),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(face = "bold", color = "black", size = 15),
    axis.text.y = element_text(face = "bold", color = "black", size = 15),
    axis.title.x = element_text(face = "bold", color = "black", size = 15),
    axis.title.y = element_text(face = "bold", color = "black", size = 15)
  ) +
  geom_label_repel(
    data = dat_rep,
    aes(label = gene_id),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE
  ) +
  labs(
    x = "log2FoldChange",
    y = "-log10(FDR)",
    title = ''
  )

# Extract differentially expressed genes for MA plot
names(list_degs_1)
list_degs_1[[1]] |> plot_ma()
list_degs_1[[2]] |> plot_ma()
list_degs_1[[3]] |> plot_ma()
list_degs_1[[4]] |> plot_ma()

######################### Step 5: GO Enrichment Analysis and Plotting #####################
library(DOSE)
library(clusterProfiler)
library(org.Dr.eg.db)
library(GOSemSim)

# Load custom organism database
org.my.db <- loadDb("Mono_albus.OrgDb")
columns(org.my.db)

test.2 <- bitr(geneID = test_deg$gene_id,
               fromType = 'GENENAME',
               toType = 'ENTREZID',
               OrgDb = org.my.db)

keys(org.my.db, keytype = "ENTREZID")

test_deg <- df |> filter(level != 'Nonsig')

go_res <- enrichGO(gene=test_deg$ENTREZID,  # gene list (converted IDs)
                   keyType="ENTREZID",      # gene ID type
                   OrgDb=org.SMonopterusAlbus.eg.db, # organism-specific OrgDb object
                   ont="ALL",              # CC, MF, BP or ALL
                   pvalueCutoff=0.05,      # P-value cutoff; may need adjustment if too few significant terms
                   pAdjustMethod="fdr",    # multiple testing correction
                   minGSSize=10,           # minimum gene set size
                   maxGSSize=500,          # maximum gene set size
                   qvalueCutoff=0.05,      # q-value cutoff; may need adjustment if too few terms
                   readable=TRUE)          # convert gene IDs to gene symbols

# Extract GO result table
res <- go_res@result
dim(res)

# Assign GO category (BP/CC/MF) based on Aspect column
df_tpm <- res |> 
  left_join(df_genome_go |> 
              dplyr::select(GO_ID, Aspect), by = c('ID'='GO_ID')) |>
  distinct(ID, .keep_all = TRUE) 
  dplyr::select(ID, Aspect)
dim(df_tpm)

df_tpm <- df_tpm |> 
  mutate(ONTOLOGY = if_else(
    Aspect == 'C', 'CC', 
    if_else(Aspect == 'F', 'MF', 'BP')
  ))

# Read GO enrichment result files for each organ (filtered by adjusted P < 0.05)
df_intestine_go <- fread('01.raw_data/Intestine_go_enrich_all.xls') |> filter(p.adjust < 0.05)
df_kidney_go    <- fread('01.raw_data/kidney_go_enrich_all.xls')   |> filter(p.adjust < 0.05)
df_spleen_go    <- fread('01.raw_data/Spleen_go_enrich_all.xls')   |> filter(p.adjust < 0.05)
df_liver_go     <- fread('01.raw_data/Liver_go_enrich_all.xls')    |> filter(p.adjust < 0.05)

list_degs_go <- list(df_intestine_go, df_kidney_go, df_spleen_go, df_liver_go)
names(list_degs_go) <- names(list_degs_1)

# Check table dimensions for each GO result
lapply(list_degs_go, dim)

# Standardize ONTOLOGY labels
list_degs_go_1 <- lapply(list_degs_go, function(x) {
  x$ONTOLOGY <- str_replace_all(
    x$ONTOLOGY,
    c(
      'biological_process'   = 'BP',
      'cellular_component'   = 'CC',
      'molecular_function'   = 'MF'
    )
  )
  return(x)
})

# GO visualization wrapper function for all organ groups
plot_GO_all <- function(go_list, go_res, out_dir = "GO_Plots", my_colors) {
  
  dir.create(out_dir, showWarnings = FALSE)
  
  library(clusterProfiler)
  library(simplifyEnrichment)
  library(ggplot2)
  library(stringr)
  library(dplyr)
  
  for (i in seq_along(go_list)) {
    df <- go_list[[i]]
    sample_id <- names(go_list)[i]  # Get corresponding group name
    
    message("Start processing sample: ", sample_id)
    
    # ① Build GO semantic similarity tree and draw treeplot
    go_res@result <- df
    edx2 <- pairwise_termsim(go_res)
    
    p1 <- treeplot(edx2,
                   hclust_method = 'average',
                   nCluster = 5,
                   group_color = my_colors)
    
    ggsave(filename = file.path(out_dir, paste0("GO_", sample_id, "_treeplot.pdf")),
           plot = p1, width = 8, height = 6)
    
    # ② simplifyGO plots for BP/CC/MF separately
    mat_bp <- GO_similarity(df$ID, ont = "BP")
    mat_cc <- GO_similarity(df$ID, ont = "CC")
    mat_mf <- GO_similarity(df$ID, ont = "MF")
    
    pdf(file.path(out_dir, paste0("GO_", sample_id, "_BP.pdf")), width = 8, height = 6)
    simplifyGO(mat_bp)
    dev.off()
    
    pdf(file.path(out_dir, paste0("GO_", sample_id, "_CC.pdf")), width = 8, height = 6)
    simplifyGO(mat_cc)
    dev.off()
    
    pdf(file.path(out_dir, paste0("GO_", sample_id, "_MF.pdf")), width = 8, height = 6)
    simplifyGO(mat_mf)
    dev.off()
    
    # ③ Simplify enriched GO terms and construct network visualization
    sim_res <- clusterProfiler::simplify(go_res, cutoff = 0.7, by = "p.adjust", measure = "Wang")
    
    gosim_res <- data.frame(sim_res@result)
    gosim_res$Description <- str_to_sentence(gosim_res$Description)
    gosim_res <- gosim_res %>%
      mutate(negLogFDR = -log10(p.adjust))
    
    p3 <- enrichmentNetwork(gosim_res,
                            drawEllipses = TRUE,
                            colorBy = 'negLogFDR',
                            colorByType = 'negLogFDR',
                            nodeSize = 'Count',
                            label_repel = TRUE) +
      scale_color_gradientn(
        colors = my_colors,
        name = "-log10(P)"
      ) +
      guides(color = guide_colorbar(title = "-log10(P)", barwidth = 1, barheight = 10),
             size = guide_legend(title = "Gene Count")) +
      theme(text = element_text(size = 10))
    
    ggsave(filename = file.path(out_dir, paste0("GO_", sample_id, "_network.pdf")),
           plot = p3, width = 8, height = 6)
    
    message("✅ All GO plots saved for: ", sample_id)
  }
}

plot_GO_all(go_list = list_degs_go_1,
            go_res = go_res,
            my_colors = diverge_hcl(5, 'Cyan-Magenta'))

######################### Step 6: KEGG Enrichment and Visualization ######################### 
df_intestine_kegg <- fread('01.raw_data/intestine_kegg_enrich_all.xls') |>
  filter(pvalue < 0.05)
df_spleen_kegg <- fread('01.raw_data/Spleen_kegg_enrich_all.xls') |>
  filter(pvalue < 0.05)
df_kidney_kegg <- fread('01.raw_data/kidney_kegg_enrich_all.xls') |>
  filter(pvalue < 0.05)
df_liver_kegg <- fread('01.raw_data/Liver_kegg_enrich_all.xls') |>
  filter(pvalue < 0.05)

list_kegg <- list(df_intestine_kegg, df_spleen_kegg, df_kidney_kegg, df_liver_kegg)
names(list_kegg) <- c('intestine', 'spleen', 'kidney', 'liver')
lapply(list_kegg, dim)

library(ggridges)

# Function to plot KEGG ridge plots per pathway
plot_kegg_ridge <- function(kegg_df, deg_df, output_file = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggridges)
  library(stringr)
  
  # Remove duplicated gene IDs from DEG table
  deg_df <- deg_df %>% distinct(gene_id, .keep_all = TRUE)
  
  # Expand geneID column, map to logFC, and prepare data for ridge plot
  df_ridges <- kegg_df %>%
    filter(p.adjust < 0.05) %>%
    mutate(log10P = -log10(p.adjust)) %>%
    separate_rows(geneID, sep = "/") %>%
    left_join(deg_df %>% select(gene_id, logFC), by = c("geneID" = "gene_id")) %>%
    select(Description, geneID, log10P, logFC)
  
  # For each pathway (Description), use average log10P as the fill color value
  df_ridges <- df_ridges %>%
    group_by(Description) %>%
    mutate(log10P_fill = mean(log10P, na.rm = TRUE)) %>%
    ungroup()
  
  # Custom color palette
  custom_colors <- colorRampPalette(c("#3daeb7", "#eeeeee", "#8075ad"))(100)
  
  p <- ggplot(df_ridges, aes(
    x = logFC,
    y = str_wrap(Description, 40),
    fill = log10P_fill
  )) +
    geom_density_ridges(alpha = 0.8, scale = 1.5, rel_min_height = 0.01,
                        color = "gray30", linewidth = 0.2) +
    labs(x = "log2 Fold Change", y = "") +
    scale_fill_gradientn(colors = custom_colors) +
    scale_y_discrete(position = "right") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    )
  
  # Save plot if file path provided
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = 10, height = 6)
  }
  
  return(p)
}

# Draw KEGG ridge plots for each tissue
plot_kegg_ridge(kegg_df = df_intestine_kegg,
                deg_df = list_degs_1[[1]],
                output_file = "KEGG_ridge_intestine.pdf")
plot_kegg_ridge(kegg_df = df_kidney_kegg |> head(20),
                deg_df = list_degs_1[[2]],
                output_file = "KEGG_ridge_kidney.pdf")
plot_kegg_ridge(kegg_df = df_spleen_kegg |> head(20),
                deg_df = list_degs_1[[3]],
                output_file = "KEGG_ridge_spleen.pdf")
plot_kegg_ridge(kegg_df = df_liver_kegg |> head(20),
                deg_df = list_degs_1[[4]],
                output_file = "KEGG_ridge_liver.pdf")

######################### Step 7: Venn and UpSet Plots #####################
names(list_degs_1)

# Keep only significant DEGs (exclude "Nonsig") and remove gene ID duplicates
list_degs_2 <- lapply(list_degs_1, function(x){
  df <- x |> filter(level != 'Nonsig')
  df <- df |> distinct(gene_id, .keep = TRUE)
  return(df)
})
lapply(list_degs_2, dim)

# Extract gene lists for each tissue
gene_list <-  list(
  Intestine = list_degs_2$Intestine$gene_id,
  Spleen    = list_degs_2$Spleen$gene_id,
  Liver     = list_degs_2$Liver$gene_id,
  Kidney    = list_degs_2$Kidney$gene_id
)

# Draw Venn diagram
ggvenn(gene_list, columns = names(gene_list), fill_color = c("#FFC566",  "#C1D9A1",  "#BDD7EE", "#FFCDBA")) 

# Prepare data for UpSet plot
gene_long <- enframe(gene_list, name = "Group", value = "Gene") %>%
  unnest(Gene) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Group, values_from = value, values_fill = 0)

library(ComplexUpset)
library(tidyverse)

# Convert to long then back to wide for UpSet input
gene_long_df <- gene_long %>%
  pivot_longer(-Gene, names_to = "Tissue", values_to = "Presence") %>%
  filter(Presence == 1) %>%
  select(-Presence)

gene_upset_df <- gene_long_df %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Tissue, values_from = value, values_fill = 0)

# Color settings for each tissue
set_colors <- c("Intestine" = "#E0F2E6", 
                "Spleen"    = "#8BC1AF", 
                "Liver"     = "#398E88", 
                "Kidney"    = "#005D67") 

library(ComplexUpset)
library(colorspace)
custom_palette <- colorRampPalette(c('#ffb3b5', "#ffd9da"))

# Draw UpSet plot (using custom set sizes color and intersection fill)
upset(
  gene_upset_df,
  intersect = names(group_colors),
  name = 'Gene combinations',
  width_ratio = 0.2,
  sort_sets = 'descending',
  sort_intersections = 'descending',
  stripes = 'white',
  set_sizes = upset_set_size(
    geom = geom_bar(fill = group_colors, width = 0.8)
  ),
  matrix = intersection_matrix(
    geom = geom_point(size = 4)
  ),
  queries = list(
    upset_query(set = 'Intestine', fill = group_colors["Intestine"], color = group_colors["Intestine"]),
    upset_query(set = 'Spleen', fill = group_colors["Spleen"], color = group_colors["Spleen"]),
    upset_query(set = 'Liver', fill = group_colors["Liver"], color = group_colors["Liver"]),
    upset_query(set = 'Kidney', fill = group_colors["Kidney"], color = group_colors["Kidney"])
  ),
  base_annotations = list(
    'Intersection size' = intersection_size(
      bar_number_threshold = 1,
      text_colors = c(on_background = 'black', on_bar = 'white'),
      fill = custom_palette(15) |> rev()
    )
  )
)

# Compute intersection shared by all four tissues
names(gene_list)
result <- Reduce(intersect, list(gene_list$Intestine,
                                 gene_list$Spleen,
                                 gene_list$Liver,
                                 gene_list$Kidney))

# Convert intersecting gene IDs and annotate
result <- tibble(gene_id = result |> as.vector())
result <- result |> 
  left_join(df_genome_gff, by = c('gene_id' = 'gene_ID'))
dim(result)

# GO enrichment analysis for intersecting DEGs
go_degs_intersect <- enrichGO(gene=result$ENTREZID, 
                              keyType="ENTREZID", 
                              OrgDb=org.SMonopterusAlbus.eg.db, 
                              ont="ALL",
                              pvalueCutoff=0.05, 
                              pAdjustMethod="fdr", 
                              minGSSize=10,
                              maxGSSize=500, 
                              qvalueCutoff=0.05,
                              readable=TRUE)

library(aPEAR)
library(stringr)
gosim_res <- data.frame(go_degs_intersect@result)
gosim_res$Description <- str_to_sentence(gosim_res$Description)
gosim_res <- gosim_res %>%
  mutate(negLogFDR = -log10(p.adjust))

aPEAR::enrichmentNetwork(gosim_res,
                         drawEllipses = TRUE,
                         colorBy = 'negLogFDR',
                         colorByType = 'negLogFDR',
                         nodeSize = 'Count',
                         label_repel = TRUE,
                         requireTheme = TRUE)

######################### Step 8: multiWGCNA Analysis #####################
 加载R包
library(tidyverse)
library(data.table)
library(ggplot2)
library(PCAtools)
library(colorspace)
library(pheatmap)
library(multiWGCNA)

# Read expression and PCA data
data_path <- '01.raw_data/'
df_pca_res <- fread(paste0(data_path, 'tpm_group_pca.txt'))
df_tpm_all <- fread(paste0(data_path, 'gene_tpm_all_samples_anno.xls'))

# Data preprocessing for WGCNA
df <- df_tpm_all |> dplyr::select(-2)
df <- df |> column_to_rownames("gene_id") |> as.data.frame()

metadata <- data.frame(
  Sample = colnames(df),
  Status = if_else(grepl('^P_', colnames(df)), 'Treatment', 'Controls'),
  Tissue = if_else(grepl('_I_', colnames(df)), 'Intestine',
                   if_else(grepl('_K_', colnames(df)), 'Kidney',
                           if_else(grepl('_L_', colnames(df)), 'Liver', 'Spleen')))
)

# Filter out genes with zero variance and keep top 75% most variable genes
df_filter <- df[apply(df, 1, var) != 0, ]
gene_variances <- apply(df_filter, 1, var)
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:round(length(gene_variances) * 0.75)]
df_filter <- df_filter[top_genes, ]

# Log-transform expression matrix
experdf <- log1p(df_filter)
sampleTable <- metadata
conditions1 <- unique(metadata[,2])
conditions2 <- unique(metadata[,3])

# Perform network construction, module eigengene calculation, and module–trait correlation
networks = constructNetworks(experdf, 
                             sampleTable, 
                             conditions1, 
                             conditions2,
                             networkType = "unsigned", 
                             power = 10,
                             minModuleSize = 40, 
                             maxBlockSize = 25000,
                             reassignThreshold = 0, 
                             minKMEtoStay = 0.7,
                             mergeCutHeight = 0.10, 
                             numericLabels = TRUE,
                             pamRespectsDendro = FALSE, 
                             verbose=3,
                             saveTOMs = FALSE)

# Compare modules by overlap across conditions
results = list()
results$overlaps = iterate(networks, overlapComparisons, plot=TRUE)

# Check reciprocal best matches between treatment and control networks
head(results$overlaps$combined_vs_Treatment)

# Run differential module expression analysis (DME) on combined networks
results$diffModExp = runDME(networks[["combined"]], 
                            sampleTable, 
                            p.adjust="fdr", 
                            refCondition="Tissue", 
                            testCondition="Status",
                            plot=TRUE, 
                            out="combined_DME.pdf")

# Calculate preservation statistics for disease-related networks
results$preservation_disease = iterate(networks[conditions1], 
                                       preservationComparisons, 
                                       write=FALSE, 
                                       plot=TRUE, 
                                       nPermutations=10)

results$diffModExp |> View()
networks$Treatment@datExpr$dynamicLabels |> table()

# Calculate preservation statistics for tissue-specific networks
results$preservation_tissue = iterate(networks[conditions2], 
                                      preservationComparisons, 
                                      write=FALSE, 
                                      plot=TRUE, 
                                      nPermutations=10)

# Extract module genes for module "Treatment_015"
module_df <- networks$Treatment@datExpr
dM15 <- module_df |> filter(dynamicLabels == 'Treatment_015') |> pull(X)

# Draw module network
library(igraph)
drawMultiWGCNAnetwork2()

# KEGG/GO enrichment for module M15
library(clusterProfiler)
library(org.SMonopterusAlbus.eg.db)

gene_id <- bitr(geneID = str_remove_all(dM15, '^gene-'),
                fromType = 'SYMBOL',
                toType = 'ENTREZID',
                OrgDb = org.SMonopterusAlbus.eg.db)

go_res <- enrichGO(gene=gene_id$ENTREZID, 
                   keyType="ENTREZID", 
                   OrgDb=org.SMonopterusAlbus.eg.db,
                   ont="ALL",
                   pvalueCutoff=0.05, 
                   pAdjustMethod="fdr", 
                   minGSSize=10,
                   maxGSSize=500, 
                   qvalueCutoff=0.05,
                   readable = TRUE)

kegg <- enrichKEGG(
  gene = gene_id$ENTREZID,
  organism = "malb",  # KEGG organism code for Monopterus albus
  keyType = "kegg",   # or "ncbi-geneid", depending on the gene ID format
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

write_csv(x = go_res@result,   file = '02.result/10.mlituwgcna/dM15_go_res.csv')
write_csv(x = kegg@result,     file = '02.result/10.mlituwgcna/dM15_kegg_res.csv')

# Font settings for plotting
library(showtext)
font_add("Arial", regular = "/home/workdir/data/fonts/ARIAL.TTF",
         bold    = "/home/workdir/data/fonts/ARIALBD.TTF",
         italic  = "/home/workdir/data/fonts/ARIALI.TTF",
         bolditalic = "/home/workdir/data/fonts/ARIALBI.TTF")
showtext_auto()

library(RColorBrewer)
base_pt <- 4.5
lab_pt <- 2.5
axis_pt <- 8
pt_range <- c(1,3)
title_pt <- 4.5
left_pt   <- -max(kegg@result$Count) / 60  # x coordinate of circle markers
left_lab  <- 0.1                            # x coordinate for left-side labels

# Color mapping for KEGG pathway categories
pathway_colors <- c(
  "Cellular Processes"                = "#8dd3c7",
  "Genetic Information Processing"    = "#ffffb3",
  "Human Diseases"                    = "#bebada",
  "Metabolism"                        = "#fb8072",
  "Organismal Systems"                = "#80b1d3",
  "Environmental Information Processing" = "#fdb462",
  "NA"                                = "#9b59b6"  # purple for missing category
)

{
  df_kegg <- kegg@result
  
  df_kegg <- df_kegg %>% 
    arrange(category, pvalue) %>% 
    mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  df_kegg_filter <- df_kegg |> 
    filter(pvalue < 0.05 & Count > 4)
  
  ggplot(df_kegg_filter) +
    # Bars for -log10(pvalue)
    geom_col(aes(x = -log10(pvalue), y = Description, fill = category)) +
    
    # Left-side pathway labels
    geom_text(aes(x = left_lab, y = Description, label = Description),
              hjust = 0, size = lab_pt) +
    
    # Circles indicating gene counts
    geom_point(aes(x = left_pt, y = Description, size = Count, fill = category),
               shape = 21, stroke = .1) +
    
    # Color and size scales
    scale_fill_manual(values = pathway_colors) +
    scale_size(range = pt_range, breaks = c(10, 20, 30, 40)) +
    
    # X-axis settings
    scale_x_continuous(expand = expansion(mult = c(0, .25))) +
    
    theme_minimal() + 
    theme(
      legend.position   = "none",
      text = element_text(family = 'Arial'),
      axis.text.x       = element_text(size = axis_pt),
      axis.title.x      = element_text(size = axis_pt),
      axis.text.y       = element_blank(),
      axis.title        = element_text(size = title_pt),
      panel.grid        = element_blank(),
      plot.margin       = margin(2, 28, 2, 2)   # top, right, bottom, left (pt)
    ) +
    labs(x = "-log10(pvalue)", y = NULL) +
    coord_cartesian(clip = "off") -> p1
  
  ggsave('02.result/10.mlituwgcna/kegg_dM15.pdf',
         p1,
         device = 'pdf',
         width = 60,
         height = 70,
         units = 'mm')
}  

# Save module gene IDs
saveRDS(dM15, '02.result/10.mlituwgcna/dM15.rds')

# Color palette for network visualization
colors <- c(
  '#e41a1c',   # red
  '#377eb8',   # blue
  '#4daf4a',   # green
  '#984ea3',   # purple
  '#ff7f00',   # orange
  '#ffff33',   # yellow
  '#00ced1'    # cyan (newly added)
)

# Draw module network and save as PDF
pdf("02.result/10.mlituwgcna/network.pdf",
    width  = 234 / 72,
    height = 212 / 72)

drawMultiWGCNAnetwork2(networks,
                       results$overlaps,
                       "Treatment_015",
                       design = sampleTable,
                       overlapCutoff = 10, 
                       padjCutoff = 0.0001, 
                       removeOutliers = TRUE, 
                       layout = NULL, 
                       hjust = 0.4, 
                       vjust = 0.3, 
                       width = 0.5,
                       colors = colors,
                       alpha = 0.05)

dev.off()

# GSVA analysis for module M15 to assess tissue- and disease-specific patterns
experdf <- as.matrix(experdf)
gene_list <- list(dm15 = dM15)

library(GSVA)
gsvaPar <- gsvaParam(exprData = experdf, geneSets = gene_list)
gsva.es <- gsva(gsvaPar, verbose=TRUE)

gsva_long <- gsva.es %>%
  as.data.frame() %>%
  rownames_to_column(var = "M15") %>%
  pivot_longer(-M15, names_to = "Sample", values_to = "GSVA_score")

gsva_long$Tissue <- str_split(gsva_long$Sample, pattern = '_', n = 4, simplify = TRUE)[,2]
gsva_long$disease <- str_split(gsva_long$Sample, pattern = '_', n = 4, simplify = TRUE)[,1]
gsva_long$disease <- factor(gsva_long$disease, levels = c('N','P'))

ggplot(gsva_long, aes(x = disease, y = GSVA_score, fill = disease)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, linewidth = 0.2) +
  theme_bw() +
  scale_fill_manual(values = c(N = "#a9d9d8", P = "#d95f02"),
                    name = "disease") +
  labs(title = "M15", x = "", y = "GSVA enrichment score") +
  theme(
    text = element_text(family = 'Arial'),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    legend.position = 'none'
  ) +
  facet_wrap(~ Tissue, nrow = 1) -> p1

ggsave('02.result/10.mlituwgcna/gsva_M15_tissue.pdf',
       p1,
       device = 'pdf',
       width = 180,
       height = 100,
       units = 'mm')

