####################################
# Format files for uplod to CReSCENT
# L.Richards
# 2021-07-05
####################################

#----- PACKAGES & PARSE -----#

suppressMessages(library(Seurat))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))

option_list <- list(make_option("--seurat",
                                type = "character",
                                default = NULL,
                                help = "file name",
                                metavar= "character"
                               ),
                    make_option("--inputPath",
                                type = "character",
                                default = NULL,
                                help = "path to seurat object",
                                metavar= "character"
                               )
                    )

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

file <- opt$seurat
inputPath <- opt$inputPath


#----- LOAD DATA -----#

print("Loading data...")
print(file)

# get list of seurat objects & define output file prefix
outPrefix <- gsub("_seurat.rds", "", file)
# load data
dat <- readRDS(paste0(inputPath, "/", file))


#----- (1) MARKER GENES -----#

print("Marker genes...")

# cluster data if not already done
if ("leiden" %in% colnames(dat@meta.data)){

    dat@meta.data$leiden <- paste0("C", dat@meta.data$leiden)
    Idents(dat) <- "leiden"


} else if ("seurat_clusters" %in% colnames(dat@meta.data)){

    dat@meta.data$seurat_clusters <- paste0("C", dat@meta.data$seurat_clusters)
    Idents(dat) <- "seurat_clusters"

} else if (sum(c("seurat_clusters", "leiden") %in% colnames(dat@meta.data)) == 0){

    dat <- FindNeighbors(dat, dims = 1:20)
    dat <- FindClusters(dat, resolution = 0.5)
    dat@meta.data$leiden <- paste0("C", dat@meta.data$leiden)
    Idents(dat) <- "leiden"

}

# find gene markers for clusters & format
markers <- FindAllMarkers(dat)
markers <- markers %>%
            group_by(cluster) %>%
            top_n(n = 15, wt = avg_log2FC) # extract top 15 per cluster
markers <- data.frame(markers)
markers <- markers[ ,c("gene", "cluster", "p_val", "avg_log2FC")]

# write tsv file
write.table(markers,
            file = paste0(outPrefix, "_deMarkers.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE
            )


#----- (2) METADATA -----#

print("Metadata...")

# extract metadata and format Broad Portal style

dat@meta.data$orig.ident <- sapply(strsplit(file,"_"), `[`, 1)
dat@meta.data$IntegrationMethod <- sapply(strsplit(file,"_"), `[`, 2)
meta <- data.frame(dat@meta.data)

NAME <- rownames(meta)
meta <- cbind(NAME, meta)
row2 <- sapply(meta, is.numeric)
row2 <- gsub("TRUE", "numeric", row2)
row2 <- gsub("FALSE", "group", row2)
row2[1] <- "TYPE"
meta <- rbind(row2, meta)

# write out file file
write.table(meta,
            file = paste0(outPrefix, "_meta.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE
            )



#----- (3) COORDINATES FILE -----#

print("Coordinate File...")

# extract UMAP, unless its Conos, then it will be largevis
if (unique(dat@meta.data$IntegrationMethod) == "Conos"){

    coords <- dat@reductions$largeVis@cell.embeddings
    coords <- data.frame(coords)
    Barcode <- rownames(coords)
    coords <- cbind(Barcode, coords)

} else {

    coords <- dat@reductions$umap@cell.embeddings
    coords <- data.frame(coords)
    Barcode <- rownames(coords)
    coords <- cbind(Barcode, coords)

}

# write out file
write.table(coords,
            file = paste0(outPrefix, "_coordinates.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE
          )



#----- (4) NORMALIZED EXPRESSION MATRIX -----#

print("Normalized Expression Matrix...")

# format
exp <- dat@assays$RNA@data
exp <- as.matrix(exp)
# exp <- data.frame(exp)
# GENE <- rownames(exp)
# exp <- cbind(GENE, exp)

# write out file file
write.table(exp,
            file = paste0(outPrefix, "_normalizedExpression.tsv"),
            quote = FALSE,
            sep = "\t",
            row.names = TRUE
            )


############################
print("End")
