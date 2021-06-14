######################################
# Run STACAS on Yost-BCC
# L.Richards
# June 14, 2021
######################################

library(Seurat)
library(STACAS)
library(optparse)


########################
#### PARSE OPTIONS #####
########################
option_list <- list(make_option("--kweight",
                                type = "double",
                                default = NULL,
                                help = "k.weight parameter for IntegrateData(); default 100",
                                metavar= "double"
                               ),
                     make_option("--distpct",
                                type = "double",
                                default = NULL,
                                help = "dist.pct parameter for AnchorFiltering; default 0.8",
                                metavar= "double"
                              ),
                      make_option("--vargenes",
                                  type = "integer",
                                  default = NULL,
                                  help = "number of variable genes to integrate with; default 1000",
                                  metavar= "integer"
                                  )
                      )
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

k.weight <- opt$kweight
dist.pct <- opt$distpct
var.genes.n <- opt$vargenes

#########################################
#### HARCODED USER DEFINED VARIABLES ####
#########################################
setwd("/cluster/projects/pughlab/projects/cancer_scrna_integration/integration/")

inputSeurat <- "/cluster/projects/pughlab/projects/cancer_scrna_integration/data/Yost-BCC/Yost-BCC_seurat.rds"
filePrefix <- "Yost-BCC"
sampleCol <- "SampleID"
patientCol <- "patient"
cellCol <- "CellType"

# define STACAS running paramters
var.genes.integrated.n <- 2000
ndim <- 20

##########################
#### PRINT PARAMETERS ####
##########################
print("Yost-BCC")
print("kweight=")
print(k.weight)
print("dist.pct=")
print(dist.pct)


####################
#### RUN STACAS ####
####################
# load data
dat <- readRDS(inputSeurat)

# split by sample
ref.list <- SplitObject(dat, split.by = sampleCol)

# split up reference list by sample and normalize
for (i in 1:length(ref.list)) {

    ref.list[[i]] <- NormalizeData(ref.list[[i]], verbose = FALSE)

    ref.list[[i]] <- FindVariableFeatures(ref.list[[i]],
                                          selection.method = "vst",
                                          nfeatures = var.genes.n*2,
                                          verbose = FALSE
                                         )

    mito.genes <- grep(pattern = "^MT-", rownames(ref.list[[i]]), value = TRUE)
    ribo.genes <- grep(pattern = "^RP[LS]", rownames(ref.list[[i]]), value = TRUE)

    #ref.list[[i]]@assays$RNA@var.features <- setdiff(ref.list[[i]]@assays$RNA@var.features, cellCycle.symbol)
    ref.list[[i]]@assays$RNA@var.features <- setdiff(ref.list[[i]]@assays$RNA@var.features, mito.genes)
    ref.list[[i]]@assays$RNA@var.features <- setdiff(ref.list[[i]]@assays$RNA@var.features, ribo.genes)
    ref.list[[i]]@assays$RNA@var.features <- head( ref.list[[i]]@assays$RNA@var.features, var.genes.n)

}

# Run STACAS
start <- Sys.time()
ref.anchors <- FindAnchors.STACAS(ref.list,
                                  dims=1:ndim,
                                  anchor.features=var.genes.integrated.n
                                 )

ref.anchors.filtered <- FilterAnchors.STACAS(ref.anchors,
                                             dist.pct = dist.pct
                                            )

all.genes <- row.names(ref.list[[1]])

for (i in 2:length(ref.list)) {

    all.genes <- intersect(all.genes, row.names(ref.list[[i]]))

}

mySampleTree <- SampleTree.STACAS(ref.anchors.filtered)
print(mySampleTree)


ref.integrated <- IntegrateData(anchorset = ref.anchors.filtered,
                                dims = 1:ndim,
                                features.to.integrate = all.genes,
                                sample.tree = mySampleTree,
                                preserve.order = T,
                                k.weight = k.weight
                               )


# process and cluster
ref.integrated <- ScaleData(ref.integrated, verbose = TRUE)
ref.integrated <- RunPCA(ref.integrated,
                         features = ref.integrated@assays$integrated@var.features,
                         ndims.print = 1:5,
                         nfeatures.print = 5
                        )
ref.integrated <- RunUMAP(ref.integrated,
                          reduction = "pca",
                          dims = 1:ndim,
                          seed.use = 123,
                          n.neighbors = 30,
                          min.dist = 0.3
                         )

end <- Sys.time()
end - start #

# plot results of integration
plot.name <- paste0(filePrefix, "_kweight", k.weight, "_distpct", dist.pct, "_STACAS_UMAP.pdf")
pdf(plot.name, width = 18, height = 5)
DimPlot(ref.integrated,
        reduction = "umap",
        group.by = c(sampleCol, patientCol, cellCol),
        ncol = 3
       )
dev.off()

# save results
seurat.name <- paste0(filePrefix, "_kweight", k.weight, "_distpct", dist.pct, "_STACAS_seurat.rds")
saveRDS(ref.integrated, file = seurat.name)
