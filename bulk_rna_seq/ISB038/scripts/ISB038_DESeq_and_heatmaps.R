library(DESeq2)
library(data.table)
library(RColorBrewer)
library(stringr)
library(pheatmap)
library(apeglm)
library(ggsci)
library(dplyr)
library(xfun)
library(ComplexHeatmap)
library(circlize)

### SETUP ###

# Project name here:
project_name <- "ISB038"

# Set the output directory (to hold the files created by this script)
output_dir <- file.path('')

# input file from partek for all samples:
combined_counts <- ".../ISB038/data/20250514_heatmap_data/Partek_ISB038-redo_Differential_analysis_filter_Ensembl_List_isb038_2 (1)/counts.txt"
### Heatmap code ###

## Setup z_score function for use in building heatmaps
cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}

# set this list of genes after examining the excel files (either directly:)
goi_list <- c( "POU5F1",
               "NANOG",
               "OTX2",
               "GBX2",
               "FGF8",
               "EN2",
               "EN1",
               "WNT1",
               "NEUROD1",
               "ATOH1",
               "BARHL1",
               "CBLN1",
               "ZIC1",
               "PAX6",
               "SOX9",
               "PAX2",
               "NFIA",
               "MAP2",
               "PTF1A",
               "LHX1",
               "SKOR2",
               "S100B",
               "GFAP",
               "GAD1",
               "GRIN1",
               "PDE1C",
               "RELN",
               "CALB2",
               "GABRA6"
)

##################################################
########  Reading in the files from partek #######
##################################################

########### Combined cell-lines:
# Read in file
input_data_2 <- read.csv( combined_counts,
                          sep="\t")
# give rownames from first column
rownames(input_data_2) <- input_data_2$Feature
# look at dims to get last column
dim(input_data_2)
# remove first column (rather, keep 2nd column through last column)
input_data_2 <- input_data_2[,2:46]

# now make the heatmap:
heatmapValues2 <- input_data_2[goi_list,]
heatmapData2 <- t(apply(heatmapValues2,1,cal_z_score))

first_heatmapData2 <- t(heatmapData2)

# Reorder rows so it's grouped by timepoint first, then GM23913, GM25256, and
# NCRM within each timepoint:
first_heatmapData2 <- first_heatmapData2[c(1,2,3,16,17,18,31,32,33,
                                           4,5,6,19,20,21,34,35,36,
                                           7,8,9,22,23,24,37,38,39,
                                           10,11,12,25,26,27,40,41,42,
                                           13,14,15,28,29,30,43,44,45),]

col_fun = colorRamp2(c(-1.45,0,5.22), c("blue","white","red"))
# Draw annotation blocks on timepoint first leaving space for cell line blocks:
Heatmap(first_heatmapData2,
        name = "z_score",
        col = col_fun,
        left_annotation = rowAnnotation(
          timepoint = anno_block( gp = gpar(fill = c( "blue", "red", "orange", "green3", "purple")),
                                  labels = c( "day 0", "day 14", "day 21", "day 35", "day 60"),
                                  labels_gp = gpar( col = c("white", "white", "black", "white", "white"),
                                                    fontsize = 10 )
          ),
          celline = anno_empty( border = FALSE )
        ),
        row_split = c(rep("day0",9),rep("day14",9), rep("day21",9), rep("day35", 9), rep("day60", 9)),
        row_order = rownames(first_heatmapData2),
        row_title = NULL,
        show_row_names = FALSE,
        column_order = colnames(first_heatmapData2),
        column_title = "All cell lines by day"
)

# Function to annotate a single timepoint block with cell line blocks:
cell_line_anno <- function( group ) {

  # Select this group
  seekViewport( group )

  # Get coords of large box
  loc1 = deviceLoc( x = unit( 1, "npc" ), y = unit( 1, "npc" ) )
  loc2 = deviceLoc( x = unit( 0, "npc" ), y = unit( 0, "npc" ) )

  # find small block height
  small_height = ( loc2$y - loc1$y ) / 3
  small_width = loc2$x - loc1$x

  seekViewport( "global")

  # first box: (GM23913)
  grid.rect( loc1$x,
             loc1$y,
             width = small_width,
             height = small_height,
             just = c("left", "bottom"),
             gp = gpar( fill = "blue")
             )
  grid.text( "GM23913",
             x = (loc1$x + loc2$x) * 0.5,
             y = ( loc1$y + loc2$y) * 0.5 - small_height,
             rot = 60,
             gp = gpar( col = "white", fontsize = 8)
             )

  # second box: (GM25256)
  grid.rect( loc1$x,
             loc1$y + small_height,
             width = small_width,
             height = small_height,
             just = c("left", "bottom"),
             gp = gpar( fill = "red")
  )
  grid.text( "GM25256",
             x = (loc1$x + loc2$x) * 0.5,
             y = ( loc1$y + loc2$y) * 0.5,
             rot = 60,
             gp = gpar( col = "white", fontsize = 8)
  )

  # third box: (NCRM)
  grid.rect( loc1$x,
             loc1$y + (2 * small_height),
             width = small_width,
             height = small_height,
             just = c("left", "bottom"),
             gp = gpar( fill = "orange")
  )
  grid.text( "NCRM",
             x = (loc1$x + loc2$x) * 0.5,
             y = ( loc1$y + loc2$y) * 0.5 + small_height,
             rot = 60,
             gp = gpar( col = "black", fontsize = 8)
  )
}

# Draw the annotations
cell_line_anno( "annotation_celline_1" )
cell_line_anno( "annotation_celline_2" )
cell_line_anno( "annotation_celline_3" )
cell_line_anno( "annotation_celline_4" )
cell_line_anno( "annotation_celline_5" )

# EXPORT NOW, then:
dev.off()


###################################################################
#
#       pca plot
#
#
library( plotly )
library( ggfortify )
library( tidyr )

pca_data <- read.csv( combined_counts, sep="\t")

# multistep transformation:
genes <- pca_data$Feature
df.data <- as.data.frame(t(pca_data[,-1]))
colnames(df.data) <- genes
df.data$sample <- gsub("\\.\\d+","",row.names(df.data))
df.data$cell_line <- gsub("_day\\d+\\.\\d+","",row.names(df.data))
df.data$timepoint <- gsub(".*(day\\d+).*","\\1",row.names(df.data))

pca_counts <- df.data[,1:45]
pca_res <- prcomp( pca_counts, scale. = FALSE )
p <- autoplot( pca_res, data = df.data, colour = "timepoint", shape = "cell_line", size = 6)
ggplotly(p)
head(df.data)

## 3d
prin_comp <- prcomp(pca_counts, rank. = 3 )
components <- prin_comp[["x"]]
components <- data.frame(components)
components$PC2 <- components$PC2
components$PC3 <- components$PC3
components = cbind(components, df.data$timepoint)
components = cbind(components, df.data$cell_line)
tevr <- summary(prin_comp)[["importance"]]['Proportion of Variance',]
tevr <- 100 * sum(tevr)

tit <- paste0("Total Explained Variance = ", tevr )
fig <- plot_ly( components,
                x = ~PC1,
                y = ~PC2,
                z = ~PC3,
                color = ~df.data$timepoint,
                colors = c("blue","red","orange","green","purple"),
                symbol = ~df.data$cell_line,
                symbols = c("cross","square","triangle-down")
) %>%
  add_markers( size = 12)

fig <- fig %>%
  layout(
    title = tit,
    scene = list(bgcolor = "grey75")
  )

fig

#################### With my normalization from DESeq2
# multistep transformation:

normalized_goi <- normalizedCounts[goi_list,]
df.data <- as.data.frame(t(normalizedCounts))
df.data$sample <- gsub("\\.\\d+","",row.names(df.data))
df.data$cell_line <- gsub("_day.*","",row.names(df.data))
df.data$timepoint <- gsub(".*(day\\d+).*","\\1",row.names(df.data))
dim(df.data)

pca_counts <- df.data[,1:37877]
pca_res <- prcomp( pca_counts, scale. = FALSE )
p <- autoplot( pca_res, data = df.data, colour = "timepoint", shape = "cell_line", size = 6)
ggplotly(p)

# with ellipses:
PC1 <- pca_res$x[,"PC1"]
PC2 <- pca_res$x[,"PC2"]
ggplot( df.data,
        aes( PC1,
             PC2,
             color = timepoint,
             shape = cell_line,
             group = timepoint
             )) +
  geom_point( size = 6) +
  theme_bw() +
  stat_ellipse() +
  xlab( "PC1 (42.9%)" ) +
  ylab( "PC2 (24.64%)")
##  for 3d, run the 3d code above at this point.