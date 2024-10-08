# load(file=file.path('objects', '05.rda'))

# Colon RNA
seu <- obj
seu <- SetIdent(seu, value = "Tissue")
seu <- subset(seu, Tissue  ==  "Colon")

## Recluster RNA
seu <- RunPCA(seu, npcs = 30, verbose = TRUE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
seu <- FindClusters(seu, resolution = 0.5)

## Recluster ADT
DefaultAssay(seu) <- 'ADT'
VariableFeatures(seu) <- rownames(seu[["ADT"]])
seu <- NormalizeData(seu, normalization.method = 'CLR', margin = 2)
seu <- ScaleData(seu)
seu <- RunPCA(seu, reduction.name = 'apca', verbose = TRUE)
seu <- RunUMAP(seu, reduction = "apca", reduction.name= "aumap",  dims = 1:2)

seu <- FindMultiModalNeighbors(
  seu, reduction.list = list("umap", "aumap"), 
  dims.list = list(1:2, 1:2), modality.weight.name = "RNA.weight",
)

# Colon RNA
DefaultAssay(seu) <- 'SCT'
DefaultDimReduc(seu) <- "umap"
scConf1 = createConfig(seu, maxLescreevels = 60)
makeShinyFiles(seu, scConf1,
             gene.mapping = TRUE,
             shiny.prefix = "sc1",
             shiny.dir = "colon_pbmc_uc_scRNAseq_GSE189040/",
             gex.assay = "SCT")

# Colon ACT
DefaultAssay(seu) <- 'ADT'
DefaultDimReduc(seu) <- "aumap"
scConf2 = createConfig(seu, maxLevels = 60)
makeShinyFiles(seu, scConf2,
             gene.mapping = TRUE,
             shiny.prefix = "sc2",
             shiny.dir = "colon_pbmc_uc_scRNAseq_GSE189040/"
             gex.assay = "SCT")

# PBMC RNA
seu <- obj
seu <- SetIdent(seu, value = "Tissue")
seu <- subset(seu, Tissue  ==  "PBMC")

## Recluster RNA
seu <- FindVariableFeatures(seu)
seu <- RunPCA(seu, npcs = 30, verbose = TRUE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:20)
seu <- FindClusters(seu, resolution = 0.5)

## Recluster ADT
DefaultAssay(seu) <- 'ADT'
VariableFeatures(seu) <- rownames(seu[["ADT"]])
seu <- NormalizeData(seu, normalization.method = 'CLR', margin = 2)
seu <- ScaleData(seu)
seu <- RunPCA(seu, reduction.name = 'apca', verbose = TRUE)
seu <- RunUMAP(seu, reduction = "apca", reduction.name= "aumap",  dims = 1:2)

seu <- FindMultiModalNeighbors(
  seu, reduction.list = list("umap", "aumap"), 
  dims.list = list(1:2, 1:2), modality.weight.name = "RNA.weight"
)

## PBMC SCT
DefaultAssay(seu) <- 'SCT'
DefaultDimReduc(seu) <- "umap"
scConf3 = createConfig(seu)
makeShinyFiles(seu, scConf3,
             gene.mapping = TRUE,
             shiny.prefix = "sc3",
             shiny.dir = "colon_pbmc_uc_scRNAseq_GSE189040/",
             gex.assay = "SCT")

# PBMC ACT
DefaultAssay(seu) <- 'ADT'
DefaultDimReduc(seu) <- "aumap"
scConf4 = createConfig(seu)
makeShinyFiles(seu, scConf4,
             gene.mapping = TRUE,
             shiny.prefix = "sc4",
             shiny.dir = "colon_pbmc_uc_scRNAseq_GSE189040/",
             gex.assay = "SCT")