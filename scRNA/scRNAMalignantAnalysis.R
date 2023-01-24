
# load library
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)
library(ggsignif)
library(circlize)
library(future)
library(infercnv)
library(corrplot)
library(clusterProfiler)
library(genomicInstability)
library(org.Hs.eg.db)
library(GSVA)
library(GSEABase)
library(msigdbr)
library(fgsea)
library(ggrepel)

global_seurat_data <- "/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds"
sample_file <- "/P04_MalignantCell/P999_Malignant_Analysis/data/sample_list.txt"
malignant_signature_file <- "/signature/Malignant/malignant_signature.txt"
outdir <- "/P04_MalignantCell/P999_Malignant_Analysis"


cat("Loading Global Seurat...\n")
global_sce <- readRDS(global_seurat_data)
cat("Complete.\n")

plan("multicore", workers = 16)



tumor_sample_colors <- c(
    ZP01T = "#E5D2DD",
    ZP02T = "#F1BB72",
    ZP03T = "#D6E7A3",
    ZP04T = "#476D87",
    ZP05T = "#E59CC4",
    ZP06T = "#23452F",
    ZP07T = "#8C549C",
    ZP08T = "#9FA3A8",
    ZP09T = "#5F3D69",
    ZP10T = "#58A4C3",
    ZP11T = "#F7F398",
    ZR01T = "#E63863",
    ZR02T = "#C1E6F3",
    ZR03T = "#91D0BE",
    ZR04T = "#712820",
    ZR05T = "#CCE0F5",
    ZR06T = "#625D9E",
    ZR07T = "#3A6963"
)


sample_colors <- c(
    ZP01T = "#E5D2DD",
    ZP01N = "#53A85F",
    ZP02T = "#F1BB72",
    ZP02N = "#F3B1A0",
    ZP03T = "#D6E7A3",
    ZP03N = "#57C3F3",
    ZP04T = "#476D87",
    ZP04N = "#E95C59",
    ZP05T = "#E59CC4",
    ZP05N = "#AB3282",
    ZP06T = "#23452F",
    ZP06N = "#BD956A",
    ZP07T = "#8C549C",
    ZP07N = "#585658",
    ZP08T = "#9FA3A8",
    ZP08N = "#E0D4CA",
    ZP09T = "#5F3D69",
    ZP09N = "#C5DEBA",
    ZP10T = "#58A4C3",
    ZP10N = "#E4C755",
    ZP11T = "#F7F398",
    ZP11N = "#AA9A59",
    ZR01T = "#E63863",
    ZR01N = "#E39A35",
    ZR02T = "#C1E6F3",
    ZR02N = "#6778AE",
    ZR03T = "#91D0BE",
    ZR03N = "#B53E2B",
    ZR04T = "#712820",
    ZR04N = "#DCC1DD",
    ZR05T = "#CCE0F5",
    ZR05N = "#CCC9E6",
    ZR06T = "#625D9E",
    ZR06N = "#68A180",
    ZR07T = "#3A6963",
    ZR07N = "#968175"
)

group_cols <- c(
    PHN = "#ef9020",
    PHT = "#00af3e",
    LRN = "#0081b4",
    LRT = "#CD1E24"
)


group_cols <- c(
    PHT = "#00af3e",
    LRT = "#CD1E24"
)

# Figure 1: inferCNV by sample
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure1")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("EPCAM+")]
    
    malignant_sce <- subset(malignant_sce, group_id %in% c("PHT", "LRT"))

    counts <- as.matrix(malignant_sce@assays$RNA@counts)

    # gene name translate
    gene <- row.names(counts)
    gene.df <- bitr(gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                    toType = c("ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                    OrgDb = org.Hs.eg.db)

    entrezid_array <- gene.df$ENTREZID
    names(entrezid_array) <- gene.df$SYMBOL
    counts <- counts[gene.df$SYMBOL,]
    row.names(counts) <- as.character(entrezid_array[row.names(counts)])

    # plot gis density plot by gis score
    cnv <- inferCNV(counts, species = c("human"))
    cnv_gis <- genomicInstabilityScore(cnv)
    malignant_sce <- AddMetaData(malignant_sce,
        cnv_gis$gis,
        col.name = "cnv_gis"
    )

    sample_array <- as.character(names(table(malignant_sce$orig.ident)[table(malignant_sce$orig.ident) >= 100]))


    for (sample_id in sample_array) {
        cat("Processing", sample_id, "\n")
        sub_malignant_sce <- malignant_sce[, malignant_sce$orig.ident %in% c(sample_id)]
        sub_malignant_sce <- subset(sub_malignant_sce,cnv_gis > -0.5)
        cat("Sample:", sample_id, nrow(sub_malignant_sce@meta.data),"\n")

        if (nrow(sub_malignant_sce@meta.data) > 100) {
            sample_dir <- file.path(fig_outdir,sample_id)
            if (!file.exists(sample_dir)) {
                dir.create(sample_dir)
            }
            setwd(sample_dir)

            merge_sce <- sub_malignant_sce

            sce.subcluster.Fibroblast <- global_sce[, global_sce$clusters %in% c("HSC")]
            # 提取500个Fibroblast细胞
            if (length(colnames(sce.subcluster.Fibroblast)) > 200) {
                sce.subcluster.Fibroblast <- sce.subcluster.Fibroblast[, sample(colnames(sce.subcluster.Fibroblast), size = 200, replace = F)]
            }

            sce.subcluster.BCell <- global_sce[, global_sce$clusters %in% c("Plasma/B")]
            # 提取500个BCell细胞
            if (length(colnames(sce.subcluster.BCell)) > 200) {
                sce.subcluster.BCell <- sce.subcluster.BCell[, sample(colnames(sce.subcluster.BCell), size = 200, replace = F)]
            }

            infercnv_info_array <- merge_sce@meta.data$orig.ident
            merge_sce <- AddMetaData(merge_sce,
                infercnv_info_array,
                col.name = "infercnv_info"
            )

            infercnv_info_array <- sce.subcluster.Fibroblast@meta.data$clusters
            sce.subcluster.Fibroblast <- AddMetaData(sce.subcluster.Fibroblast,
                infercnv_info_array,
                col.name = "infercnv_info"
            )

            infercnv_info_array <- sce.subcluster.BCell@meta.data$clusters
            sce.subcluster.BCell <- AddMetaData(sce.subcluster.BCell,
                infercnv_info_array,
                col.name = "infercnv_info"
            )

            # 合并Epithelial细胞和Fibroblast细胞
            sce_merge <- merge(merge_sce,
                y = sce.subcluster.Fibroblast,
                add.cell.ids = NULL, project = "merge"
            )

            sce_merge <- merge(sce_merge,
                y = sce.subcluster.BCell,
                add.cell.ids = NULL, project = "merge"
            )

            counts <- as.matrix(sce_merge@assays$RNA@counts)
            write.table(counts, file = paste0(sample_id, ".merge_express.matrix.txt"), sep = "\t", quote = FALSE)
            write.table(sce_merge[["infercnv_info"]],
                file = paste0(sample_id,".merge_cell.cluster.txt"),
                sep = "\t",
                quote = FALSE,
                col.names = FALSE
            )
        }
    }
}

# By Sample cluster_by_groups = T
if (F) {
    library(Seurat)
    library(infercnv)

    fig_outdir <- paste0(outdir, "/", "Figure1")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    for (sample_id in list.files(fig_outdir,"Z")) {
        setwd(file.path(fig_outdir,"cluster_by_groups_is_TRUE", sample_id))
        cat("Parsing", sample_id,"\n")
        ######### infercnv分析##########
        infercnv_obj <- CreateInfercnvObject(
            raw_counts_matrix = paste0(sample_id,".merge_express.matrix.txt"),
            annotations_file = paste0(sample_id,".merge_cell.cluster.txt"),
            delim = "\t",
            gene_order_file = "/seurat_jerry/hcc_project/P04_MalignantCell/P999_Malignant_Analysis/data/gencode_v19_gene_pos.txt",
            ref_group_names = c("HSC")
        )

        cat("Running infercnv...\n")
        infercnv_obj <- infercnv::run(
            infercnv_obj,
            cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
            out_dir = file.path(fig_outdir, sample_id), # dir is auto-created for storing outputs
            # cluster_by_groups=FALSE,   # cluster
            cluster_by_groups = T,
            denoise = TRUE,
            HMM = F,
            num_threads = 28,
            analysis_mode = "subclusters"
        )
    }
}


# Figure 1-2 inferCNV heatmap
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure1-2")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("EPCAM+")]

    group_id_array <- paste0(
        substr(malignant_sce@meta.data$sample_id, 1, 2),
        substr(malignant_sce@meta.data$sample_id, 5, 5)
    )
    malignant_sce <- AddMetaData(malignant_sce,
        group_id_array,
        col.name = "group_id"
    )

    malignant_sce <- subset(malignant_sce, group_id %in% c("ZPT", "ZRT"))

    counts <- as.matrix(malignant_sce@assays$RNA@counts)

    # gene name translate
    gene <- row.names(counts)
    gene.df <- bitr(gene, fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                    toType = c("ENTREZID"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                    OrgDb = org.Hs.eg.db)

    entrezid_array <- gene.df$ENTREZID
    names(entrezid_array) <- gene.df$SYMBOL
    counts <- counts[gene.df$SYMBOL,]
    row.names(counts) <- as.character(entrezid_array[row.names(counts)])

    # plot gis density plot by gis score
    cnv <- inferCNV(counts, species = c("human"))
    cnv_gis <- genomicInstabilityScore(cnv)
    malignant_sce <- AddMetaData(malignant_sce,
        cnv_gis$gis,
        col.name = "cnv_gis"
    )

    sample_array <- as.character(names(table(malignant_sce$orig.ident)[table(malignant_sce$orig.ident) >= 100]))


    for (sample_id in sample_array) {
        cat("Processing", sample_id, "\n")
        sub_malignant_sce <- malignant_sce[, malignant_sce$orig.ident %in% c(sample_id)]
        sub_malignant_sce <- subset(sub_malignant_sce,cnv_gis > -0.5)
        cat("Sample:", sample_id, nrow(sub_malignant_sce@meta.data),"\n")

        if (nrow(sub_malignant_sce@meta.data) > 100) {
            sample_dir <- file.path(fig_outdir,sample_id)
            if (!file.exists(sample_dir)) {
                dir.create(sample_dir)
            }
            setwd(sample_dir)

            merge_sce <- sub_malignant_sce

            sce.subcluster.Fibroblast <- global_sce[, global_sce$clusters %in% c("HSC")]
            # 提取500个Fibroblast细胞
            if (length(colnames(sce.subcluster.Fibroblast)) > 200) {
                sce.subcluster.Fibroblast <- sce.subcluster.Fibroblast[, sample(colnames(sce.subcluster.Fibroblast), size = 200, replace = F)]
            }

            infercnv_info_array <- merge_sce@meta.data$orig.ident
            merge_sce <- AddMetaData(merge_sce,
                infercnv_info_array,
                col.name = "infercnv_info"
            )

            infercnv_info_array <- sce.subcluster.Fibroblast@meta.data$clusters
            sce.subcluster.Fibroblast <- AddMetaData(sce.subcluster.Fibroblast,
                infercnv_info_array,
                col.name = "infercnv_info"
            )

            # 合并Epithelial细胞和Fibroblast细胞
            sce_merge <- merge(merge_sce,
                y = sce.subcluster.Fibroblast,
                add.cell.ids = NULL, project = "merge"
            )


            counts <- as.matrix(sce_merge@assays$RNA@counts)
            write.table(counts, file = paste0(sample_id, ".merge_express.matrix.txt"), sep = "\t", quote = FALSE)
            write.table(sce_merge[["infercnv_info"]],
                file = paste0(sample_id,".merge_cell.cluster.txt"),
                sep = "\t",
                quote = FALSE,
                col.names = FALSE
            )
        }
    }
}


# By Sample cluster_by_groups = F
if (T) {
    library(Seurat)
    library(infercnv)

    fig_outdir <- paste0(outdir, "/", "Figure1-2")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    for (sample_id in list.files(fig_outdir,"Z")) {
        setwd(file.path(fig_outdir, sample_id))
        cat("Parsing", sample_id,"\n")
        ######### infercnv分析##########
        infercnv_obj <- CreateInfercnvObject(
            raw_counts_matrix = paste0(sample_id,".merge_express.matrix.txt"),
            annotations_file = paste0(sample_id,".merge_cell.cluster.txt"),
            delim = "\t",
            gene_order_file = "/seurat_jerry/hcc_project/P04_MalignantCell/P999_Malignant_Analysis/data/gencode_v19_gene_pos.txt",
            ref_group_names = c("HSC")
        )

        cat("Running infercnv...\n")
        infercnv_obj <- infercnv::run(
            infercnv_obj,
            cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
            out_dir = file.path(fig_outdir, sample_id), # dir is auto-created for storing outputs
            # cluster_by_groups=FALSE,   # cluster
            cluster_references = FALSE,
            cluster_by_groups = F,
            denoise = TRUE,
            HMM = F,
            num_threads = 28,
            analysis_mode = "subclusters"
        )
    }
}



# Figure 2 inferCNV by group
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure2")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    if (F) {
    malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("Epithelial", "Malignant")]
    table(malignant_sce@meta.data$orig.ident)
    group_id_array <- paste0(substr(malignant_sce@meta.data$orig.ident, 1, 2), substr(malignant_sce@meta.data$orig.ident, 5, 5))
    malignant_sce <- AddMetaData(malignant_sce,
        group_id_array,
        col.name = "group_id"
    )
    malignant_sce <- subset(malignant_sce, group_id %in% c("ZPT","ZRT","ZPN","ZRN"))
    group_id_array <- malignant_sce$group_id
    sce_array <- c()
    cell_ids_array <- c()

    for (group_id in unique(group_id_array)) {
        cat("Processing", group_id, "\n")
        sub_malignant_sce <- malignant_sce[, malignant_sce$group_id %in% c(group_id)]
        if (length(colnames(sub_malignant_sce)) > 300) {
            sub_malignant_sce <- sub_malignant_sce[, sample(colnames(sub_malignant_sce),
                size = 300, replace = F
            )]
        }
        sce_array <- c(sce_array, sub_malignant_sce)
        cell_ids_array <- c(cell_ids_array, cluster)
    }

    merge_sce <- sce_array[[1]]
    for (sub_sce in sce_array){
        merge_sce <- merge(merge_sce, sub_sce, add.cell.ids = NULL,project = "merge")
    }

    sce.subcluster.Fibroblast <- global_sce[, global_sce$clusters %in% c("HSC")]
    # 提取1000个Fibroblast细胞
    if (length(colnames(sce.subcluster.Fibroblast)) > 300) {
        sce.subcluster.Fibroblast <- sce.subcluster.Fibroblast[, sample(colnames(sce.subcluster.Fibroblast), size = 300, replace = F)]
    }

    sce.subcluster.BCell <- global_sce[, global_sce$clusters %in% c("Plasma/B")]
    # 提取1000个BCell细胞
    if (length(colnames(sce.subcluster.BCell)) > 300) {
        sce.subcluster.BCell <- sce.subcluster.BCell[, sample(colnames(sce.subcluster.BCell), size = 300, replace = F)]
    }

    infercnv_info_array <- merge_sce@meta.data$group_id
    merge_sce <- AddMetaData(merge_sce,
        infercnv_info_array,
        col.name = "infercnv_info"
    )

    infercnv_info_array <- sce.subcluster.Fibroblast@meta.data$clusters
    sce.subcluster.Fibroblast <- AddMetaData(sce.subcluster.Fibroblast,
        infercnv_info_array,
        col.name = "infercnv_info"
    )

    infercnv_info_array <- sce.subcluster.BCell@meta.data$clusters
    sce.subcluster.BCell <- AddMetaData(sce.subcluster.BCell,
        infercnv_info_array,
        col.name = "infercnv_info"
    )


    # 合并Epithelial细胞和Fibroblast细胞
    sce_merge <- merge(merge_sce,
        y = sce.subcluster.Fibroblast,
        add.cell.ids = NULL, project = "merge"
    )

    sce_merge <- merge(sce_merge, y = sce.subcluster.BCell, add.cell.ids = NULL,project = "merge")


    counts <- as.matrix(sce_merge@assays$RNA@counts)
    write.table(counts, file = "merge_express.matrix.txt", sep = "\t", quote = FALSE)
    write.table(sce_merge[["infercnv_info"]],
        file = "merge_cell.cluster.txt",
        sep = "\t",
        quote = FALSE,
        col.names = FALSE
    )
    }
    ######### infercnv分析##########
    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = "merge_express.matrix.txt",
        annotations_file = "merge_cell.cluster.txt",
        delim = "\t",
        gene_order_file = "/seurat_jerry/hcc_project/P04_MalignantCell/P999_Malignant_Analysis/data/gencode_v19_gene_pos.txt",
        ref_group_names = c("HSC")
    )
    saveRDS(infercnv_obj, file = "infercnv.rds", compress = F)

    infercnv_obj <- infercnv::run(
        infercnv_obj,
        cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
        out_dir = fig_outdir, # dir is auto-created for storing outputs
        # cluster_by_groups=FALSE,   # cluster
        cluster_by_groups = T,
        denoise = TRUE,
        HMM = TRUE,
        num_threads = 32,
        analysis_mode = "subclusters"
    )

    infercnv::plot_cnv(infercnv_obj, #上两步得到的infercnv对象
                   plot_chr_scale = T, #画染色体全长，默认只画出（分析用到的）基因
                   output_filename = "final_plot",
                   output_format = "pdf", #保存为pdf文件
                   custom_color_pal =  color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2))) #改颜色

}



# FigureS4H-Figure02 Malignant Signature violinPlot between Groups
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure02")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    # get epcam+ seurat object
    malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("EPCAM+")]
    # samples_list <- read.table("/P04_MalignantCell/P999_Malignant_Analysis/data/sample_list.txt", header = F, sep = "\t")
    # group_array <- samples_list[, "V2"]
    # names(group_array) <- samples_list[, "V1"]
    # malignant_sce$group <- group_array[malignant_sce$orig.ident]
    malignant_sce <- subset(malignant_sce, group_id %in% c("PHT","LRT"))

    # sce_array <- c()
    # for (sample_id in unique(malignant_sce$orig.ident)) {
    #     sub_sce <- subset(malignant_sce, orig.ident %in% c(sample_id))
    #     print(length(Cells(sub_sce)))
    #     if (length(Cells(sub_sce)) < 200) {
    #         next
    #     } else {
    #         cell_id_array <- sample(Cells(sub_sce), 200, replace = FALSE)
    #     }
    #     sub_sce <- subset(sub_sce, cells = cell_id_array)
    #     sce_array <- c(sce_array, sub_sce)
    # }

    sub_malignant_sce <- sce_array[[1]]
    for (i in 2:length(sce_array)) {
        sub_malignant_sce <- merge(sub_malignant_sce, y = sce_array[[i]])
    }
    malignant_sce <- sub_malignant_sce
    malignant_signature_frame <- read.table(file = "/P04_MalignantCell/P999_Malignant_Analysis/data/malignant_signature.txt", sep = "\t",header = F,stringsAsFactors = F)

    malignant_signature_frame <- malignant_signature_frame %>% filter(V1 %in% c("Proliferation","Immune_surveillance"))

    split_signature_list = split(malignant_signature_frame, malignant_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        malignant_sce <- AddModuleScore(
            object = malignant_sce,
            features = list(gene_array),
            name = signature
        )
        colnames(malignant_sce@meta.data) <- gsub(colnames(malignant_sce@meta.data),
           pattern = paste0(signature, 1), replacement = signature
        )
    }
    malignant_sce$group_id <- factor(malignant_sce$group_id,levels = names(group_cols))
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
        compaired <- list(
            c("PHT", "LRT")
        )

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'group_id', log=TRUE,... ) +
            geom_boxplot(width = 0.2) +
            xlab("") + 
            ylab("") + 
            ggtitle(feature) +
            theme(legend.position = "none",
            #axis.text.x = element_text(size = rel(1), angle = 90),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(),
            axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),              
            plot.margin = plot.margin ) + 
            geom_signif(
                comparisons = compaired,
                step_increase = 0.3, 
                map_signif_level = T,
                #test.args = "greater", 
                test = wilcox.test
            )
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.25, 0, -0.25, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
                plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
                theme(axis.text.x=element_text(), axis.ticks.x = element_line())
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 2)
        return(p)
    }

    p <- StackedVlnPlot(malignant_sce,
        names(split_signature_list),
        pt.size = 0, 
        cols = group_cols
    )


    ggsave(filename = "Figure3.png", plot = p, width = 6, height = 5, units = c("in"))
    ggsave(filename = "Figure3.pdf", plot = p, width = 6, height = 5, units = c("in"))

}



#VIP Figure2B Malignant Markers Split ViolinPlot Between Groups
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure4")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    markers <- c("CD47","HLA-A","HLA-B","HLA-C","B2M")

    #orig.ident
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
        compaired <- list(
            c("PHT", "LRT")
        )

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'group_id', log=TRUE,... ) +
            xlab("") + 
            ylab("") + 
            ggtitle(feature) +
            theme(legend.position = "none",
            #axis.text.x = element_text(size = rel(1), angle = 90),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(),
            axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),              
            plot.margin = plot.margin ) + 
            geom_signif(
                comparisons = compaired,
                step_increase = 0.3, 
                map_signif_level = T,
                #test.args = "greater", 
                test = wilcox.test
            )
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
                plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
                theme(axis.text.x=element_text(), axis.ticks.x = element_line())
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 5)
        return(p)
    }

    p <- StackedVlnPlot(malignant_sce,
        markers,
        pt.size = 0, 
        cols = group_cols
    )


    ggsave(filename = "Figure4.png", plot = p, width = 10, height = 6, units = c("in"))
    ggsave(filename = "Figure4.pdf", plot = p, width = 10, height = 6, units = c("in"))

    # updated split violin plot version
    malignant_sce = AddMetaData(malignant_sce, malignant_sce[["RNA"]]@data["HLA-A", ],col.name = "HLAA")
    malignant_sce = AddMetaData(malignant_sce, malignant_sce[["RNA"]]@data["HLA-B", ],col.name = "HLAB")
    malignant_sce = AddMetaData(malignant_sce, malignant_sce[["RNA"]]@data["HLA-C", ],col.name = "HLAC")
    malignant_sce = AddMetaData(malignant_sce, malignant_sce[["RNA"]]@data["B2M", ],col.name = "B2M")
    malignant_sce = AddMetaData(malignant_sce, malignant_sce[["RNA"]]@data["HLA-DRA", ],col.name = "HLADRA")
    malignant_sce = AddMetaData(malignant_sce, malignant_sce[["RNA"]]@data["HLA-DRB1", ],col.name = "HLADRB1")

    malignant_frame <- malignant_sce@meta.data[,c("group_id","HLAA","HLAB","HLAC","B2M","HLADRA","HLADRB1")]
    malignant_frame <- malignant_frame %>% tidyr::gather("geneName","expression","HLAA","HLAB","HLAC","B2M","HLADRA","HLADRB1")
    malignant_frame$group_id <- factor(malignant_frame$group_id,levels = c("PHT","LRT"))

    compaired <- list(c("PHT", "LRT"))
    library(ggpubr)
    library(scales)
    library(introdataviz)

    gp <- ggplot(malignant_frame, aes(x = geneName, y = expression, fill = group_id)) +
        geom_split_violin(trim = T, colour = T) +
        stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
               position = position_dodge(.175)) +
        stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
        scale_fill_manual(values = group_cols)+
        theme_bw() +
        ylab("RNA exression") +
        xlab("") +
        theme(axis.text.x = element_text(size = 16),
            axis.text.y = element_blank()) +
        scale_x_discrete(breaks = c("HLAA","HLAB","HLAC","B2M","HLADRA","HLADRB1"), 
                 labels = c("HLA-A","HLA-B","HLA-C","B2M","HLA-DRA","HLA-DRB1")) +
        stat_compare_means(aes(group=group_id),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif")

    ggsave(file = "ImmuneExpressionVlnPlot.png",plot = gp, w = 18, h = 6)
    ggsave(file = "ImmuneExpressionVlnPlot.pdf",plot = gp, w = 18, h = 6)


}



# Figure 5 Malignant Markers ViolinPlot Among Samples
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure5")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    markers <- c("MYC")

    #orig.ident
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'orig.ident', log=TRUE,... ) +
            xlab("") + 
            ylab("") + 
            ggtitle(feature) +
            theme(legend.position = "none",
            #axis.text.x = element_text(size = rel(1), angle = 90),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(),
            axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),              
            plot.margin = plot.margin ) + 
            geom_signif(
                comparisons = compaired,
                step_increase = 0.3, 
                map_signif_level = T,
                test.args = "greater", 
                test = wilcox.test
            )
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
                plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
                theme(axis.text.x=element_text(), axis.ticks.x = element_line())
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
        return(p)
    }

    p <- StackedVlnPlot(malignant_sce,
        markers,
        pt.size = 0, 
        cols = sample_colors
    )


    ggsave(filename = "Figure5.png", plot = p, width = 10, height = 4, units = c("in"))
    ggsave(filename = "Figure5.pdf", plot = p, width = 10, height = 4, units = c("in"))

}



# Figure 6 Malignant and Epithelial Classification
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure6")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)

    # extract EPCAM plus cells from global seurat object
    malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("EPCAM+")]

    # get "ZPT", "ZRT" array for group info column from the meta data's orig.ident
    group_id_array <- paste0(
        substr(malignant_sce@meta.data$sample_id, 1, 2),
        substr(malignant_sce@meta.data$sample_id, 5, 5)
    )

    malignant_sce <- AddMetaData(malignant_sce,
        group_id_array,
        col.name = "group_id"
    )

    # according to the genomicInstability requirement, get the count matrix from seurat object
    counts <- as.matrix(malignant_sce@assays$RNA@counts)

    # gene name translate
    gene <- row.names(counts)

    # fromType 是指你的数据ID类型是属于哪一类的
    # toType 是指你要转换成哪种ID类型，可以写多种，也可以只写一种
    gene.df <- bitr(gene, fromType = "SYMBOL", 
                    toType = c("ENTREZID"), 
                    OrgDb = org.Hs.eg.db)

    # counts matrix行名称转成entrizid id
    entrezid_array <- gene.df$ENTREZID
    names(entrezid_array) <- gene.df$SYMBOL
    counts <- counts[gene.df$SYMBOL,]
    row.names(counts) <- as.character(entrezid_array[row.names(counts)])

    # plot gis density plot by gis score
    # 注意这里inferCNV是genomicInstability包中的function
    cnv <- inferCNV(counts, species = c("human"))
    cnv_gis <- genomicInstabilityScore(cnv)

    par(mai = c(.8, .8, .2, .2))
    png("gis_density.png", width = 1024, height = 1024)
    plot(density(cnv_gis$gis), lwd=2, xlab="GIS", main="")
    dev.off()

    # add likelihood curve to the density plot
    par(mai=c(.8, .8, .2, .8))
    png("distribution_inference.png",width = 1024,height = 1024)
    giDensityPlot(cnv_gis, ylim = c(0, .8),col = "white")
    dev.off()

    # Adding the likelihood data and second-axis
    cnv_gis <- giLikelihood(cnv_gis, recompute = F, normal = 1, tumor = 2)
    pdf("likelihood.pdf",width = 8, height = 8)
    giDensityPlot(cnv_gis, ylim = c(0, .8),col = "white")
    pos <- order(cnv_gis$gis)
    lines(cnv_gis$gis[pos],
        cnv_gis$gi_likelihood[pos] * .8,
        lwd = 2
    )
    axis(4, seq(0, .8, length=6), seq(0, 1, .2))
    axis(4, .4, "Relative likelihood", tick=FALSE, line=1.5)
    pos5 <- which.min((.5-cnv_gis$gi_likelihood)^2)
    lines(c(
        rep(cnv_gis$gis[pos5], 2),
        max(cnv_gis$gis * 1.05)
    ), c(0, rep(cnv_gis$gi_likelihood[pos5] * .8, 2)), lty = 3)

    dev.off()


    cnv_score_array <- cnv_gis$gis
    cnv_score_array <- as.numeric(cnv_score_array[row.names(malignant_sce@meta.data)])
    # cnv_status_array <- ifelse(cnv_score_array > -0.8905025, "T", "N")
    cnv_status_array <- c()
    for (score in cnv_score_array) {
        if (score >= -1.5) {
            status <- "T"
        } else if (score < -1.5) {
            status <- "N"
        } 
        cnv_status_array <- c(cnv_status_array, status)
    }

    malignant_sce <- AddMetaData(malignant_sce,
        cnv_status_array,
        col.name = "cnv_status"
    )

    malignant_sce <- subset(malignant_sce, group_id %in% c("ZPT","ZRT"))
    sce_array <- c()

    for (cnv_status in unique(cnv_status_array)) {
        cat("Processing", cnv_status, "\n")
        sub_malignant_sce <- malignant_sce[, malignant_sce$cnv_status %in% c(cnv_status)]
        if (length(colnames(sub_malignant_sce)) > 500) {
            sub_malignant_sce <- sub_malignant_sce[, sample(colnames(sub_malignant_sce),
                size = 500, replace = F
            )]
        }
        sce_array <- c(sce_array, sub_malignant_sce)
    }

    merge_sce <- sce_array[[1]]
    for (sub_sce in sce_array){
        merge_sce <- merge(merge_sce, sub_sce, add.cell.ids = NULL,project = "merge")
    }

    sce.subcluster.Fibroblast <- global_sce[, global_sce$clusters %in% c("HSC")]
    # 提取1000个Fibroblast细胞
    if (length(colnames(sce.subcluster.Fibroblast)) > 500) {
        sce.subcluster.Fibroblast <- sce.subcluster.Fibroblast[, sample(colnames(sce.subcluster.Fibroblast), size = 500, replace = F)]
    }

    sce.subcluster.BCell <- global_sce[, global_sce$clusters %in% c("Plasma/B")]
    # 提取1000个BCell细胞
    if (length(colnames(sce.subcluster.BCell)) > 500) {
        sce.subcluster.BCell <- sce.subcluster.BCell[, sample(colnames(sce.subcluster.BCell), size = 500, replace = F)]
    }

    infercnv_info_array <- merge_sce@meta.data$cnv_status
    merge_sce <- AddMetaData(merge_sce,
        infercnv_info_array,
        col.name = "infercnv_info"
    )

    infercnv_info_array <- sce.subcluster.Fibroblast@meta.data$clusters
    sce.subcluster.Fibroblast <- AddMetaData(sce.subcluster.Fibroblast,
        infercnv_info_array,
        col.name = "infercnv_info"
    )

    infercnv_info_array <- sce.subcluster.BCell@meta.data$clusters
    sce.subcluster.BCell <- AddMetaData(sce.subcluster.BCell,
        infercnv_info_array,
        col.name = "infercnv_info"
    )


    # 合并Epithelial细胞和Fibroblast细胞
    sce_merge <- merge(merge_sce,
        y = sce.subcluster.Fibroblast,
        add.cell.ids = NULL, project = "merge"
    )

    sce_merge <- merge(sce_merge, y = sce.subcluster.BCell, add.cell.ids = NULL,project = "merge")


    counts <- as.matrix(sce_merge@assays$RNA@counts)
    write.table(counts, file = "merge_express.matrix.txt", sep = "\t", quote = FALSE)
    write.table(sce_merge[["infercnv_info"]],
        file = "merge_cell.cluster.txt",
        sep = "\t",
        quote = FALSE,
        col.names = FALSE
    )

    ######### infercnv分析##########
    infercnv_obj <- CreateInfercnvObject(
        raw_counts_matrix = "merge_express.matrix.txt",
        annotations_file = "merge_cell.cluster.txt",
        delim = "\t",
        gene_order_file = "/seurat_jerry/hcc_project/P04_MalignantCell/P999_Malignant_Analysis/data/gencode_v19_gene_pos.txt",
        ref_group_names = c("HSC")
    )

    infercnv_obj <- infercnv::run(
        infercnv_obj,
        cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
        out_dir = fig_outdir, # dir is auto-created for storing outputs
        # cluster_by_groups=FALSE,   # cluster
        cluster_by_groups = T,
        denoise = TRUE,
        HMM = TRUE,
        num_threads = 32,
        analysis_mode = "subclusters"
    )




















    if (F) {
    # plot UMAP by CNV score

    cnv_score_array <- cnv_gis$gis
    cnv_score_array <- as.numeric(cnv_score_array[row.names(malignant_sce@meta.data)])
    # cnv_status_array <- ifelse(cnv_score_array > -0.8905025, "T", "N")
    cnv_status_array <- c()
    for (score in cnv_score_array) {
        if (score >= -1.2) {
            status <- "Malignant"
        } else if (score < -1.2) {
            status <- "Epithelial"
        }
        cnv_status_array <- c(cnv_status_array, status)
    }

    malignant_sce <- AddMetaData(malignant_sce,
        cnv_status_array,
        col.name = "cnv_status"
    )
    }

}



# Figure 7 Malignant GSEA By Group (ZPT vs ZRT)
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure7")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)

    malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("EPCAM+")]

    group_id_array <- paste0(substr(malignant_sce@meta.data$orig.ident, 1, 2), substr(malignant_sce@meta.data$orig.ident, 5, 5))
    malignant_sce <- AddMetaData(malignant_sce,
        group_id_array,
        col.name = "group_id"
    )
    malignant_sce <- subset(malignant_sce, group_id %in% c("PHT","LRT"))

    #hallmark_geneset <- getGmt("/seurat_jerry/hcc_project/P04_MalignantCell/P999_Malignant_Analysis/data/h.all.v6.2.symbols.gmt")

    deg <- FindMarkers(malignant_sce,
        ident.1 = "ZRT", 
        group.by = "group_id", 
        assay = "RNA", 
        slot = "counts", 
        logfc.threshold = 0, 
        min.pct = 0
    )

    saveRDS(deg,file = "ZRTvsZPT_Malignant.rds")
    # GSEA BarPlot

    deg <- readRDS("/P04_MalignantCell/P999_Malignant_Analysis/Figure7/ZRTvsZPT_Malignant.rds")
    ##通路基因集
    msgdC2 = msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
    fgsea_sets = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_description)


    deg$genes = rownames(deg)
    cluster0.genes<- deg %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
    ranks<- deframe(cluster0.genes)

    fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

    gp <- ggplot(fgseaRes %>% 
        as_tibble() %>% 
        arrange(desc(NES)) %>% 
        filter(pval < 0.05) %>% 
        head(n= 5), aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = NES)) +
            coord_flip() +
            labs(x = "", y = "Normalized Enrichment Score", 
            title = "KEGG gene set enrichment in LRT") +
            theme_bw()

    ggsave(filename = "kegg_barplot.png", width = 8, height = 6, units = c("in"))
    ggsave(filename = "kegg_barplot.pdf", width = 8, height = 6, units = c("in"))

    
    mdb_h <- msigdbr(species = "Homo sapiens", category = "H")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    deg$genes = rownames(deg)
    cluster0.genes<- deg %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
    ranks<- deframe(cluster0.genes)

    fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

    gp <- ggplot(fgseaRes %>% 
        as_tibble() %>% 
        arrange(desc(NES)) %>% 
        filter(pval < 0.05) %>% 
        head(n= 40), aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = NES)) +
            coord_flip() +
            labs(x = "Hallmark", y = "Normalized Enrichment Score", 
            title = "HALLMARK gene sets NES from GSEA") +
            theme_bw()

    ggsave(filename = "hallmark_barplot.png", width = 8, height = 6, units = c("in"))
    ggsave(filename = "hallmark_barplot.pdf", width = 8, height = 6, units = c("in"))



    # Pseudo Volcalno Plot
    log2FC_threshold = 1
    p_value_threshold = 0.01
    deg[which(deg$p_val_adj < 0.05 & deg$avg_log2FC <= -log2FC_threshold),'sig'] <- 'ZRT'
    deg[which(deg$p_val_adj < 0.05 & deg$avg_log2FC >= log2FC_threshold),'sig'] <- 'ZPT'
    deg[which(deg$p_val_adj >= 0.05 | abs(deg$avg_log2FC) < log2FC_threshold),'sig'] <- 'None'

    # 横轴 log2FC，纵轴 -log10(adj.P.Val)，颜色表示差异
    deg$Gene <- row.names(deg)
    deg$percent_difference <- deg$pct.1 - deg$pct.2
    p <- ggplot(deg, aes(x = percent_difference, y = avg_log2FC, color = sig)) +
        geom_point(alpha = 0.8, size = 0.6) +
        scale_colour_manual(values = c("#00af3e", "#CD1E24", "gray"), limits = c("ZPT", "ZRT", "None")) +
        theme(panel.grid = element_blank(), 
            panel.background = element_rect(color = "black", fill = "transparent"), 
            plot.title = element_text(hjust = 0.5)) +
        theme(legend.key = element_rect(fill = "transparent"), 
            legend.background = element_rect(fill = "transparent"), 
            legend.position = c(0.9, 0.93)) +
        xlim(-2, 2) +
        ylim(-4, 4) +
        labs(
            x = "\nPercentage Difference",
            y = "Log2Fold Change\n", 
            color = "", 
            title = "ZPT vs ZRT\n"
        )

    up <- subset(deg, sig == 'ZPT')
    up <- up[order(up$p_val_adj), ][1:60, ]
    down <- subset(deg, sig == 'ZRT')
    down <- down[order(down$p_val_adj), ][1:60, ]

    p1 <- p + theme(legend.position = 'right') +
            geom_label_repel(data = rbind(up, down), 
                    aes(x = percent_difference, y = avg_log2FC, label = Gene),
                    size = 3,
                    box.padding = unit(0.5, 'lines'), 
                    segment.color = 'black',
                    alpha = 0.7,
                    max.overlaps = Inf, 
                    show.legend = FALSE)

    ggsave('ZPT_ZRT_VolcanoPlot.png', p1, width = 6 * 2, height = 4 * 2)
    ggsave('ZPT_ZRT_VolcanoPlot.pdf', p1, width = 6 * 2, height = 4 * 2)


}




# Figure 8 Malignant and Epithelial Volcano Plot
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure8")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    epi_malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("EPCAM+")]

    group_id_array <- substr(epi_malignant_sce@meta.data$orig.ident, 5, 5)
    epi_malignant_sce <- AddMetaData(epi_malignant_sce,
        group_id_array,
        col.name = "TorN"
    ) 
    
    group_id_array <- paste0(substr(epi_malignant_sce@meta.data$orig.ident, 1, 2), substr(epi_malignant_sce@meta.data$orig.ident, 5, 5))
    epi_malignant_sce <- AddMetaData(epi_malignant_sce,
        group_id_array,
        col.name = "group_id"
    )

    markers <- FindMarkers(epi_malignant_sce,
        ident.1 = "T", 
        group.by = "TorN", 
        assay = "RNA", 
        slot = "counts", 
        logfc.threshold = 0, 
        min.pct = 0
    )

    saveRDS(markers,file = "malignant_vs_epithelial_markers.rds")

    epithelial_markers <- markers %>%
        filter(p_val_adj < 0.01, avg_log2FC < -1) %>%
        row.names() %>%
        head(n = 50)
    malignant_markers <- markers %>%
        filter(p_val_adj < 0.01, avg_log2FC > 1) %>%
        row.names() %>%
        head(n = 50)

    epi_malignant_sce <- AddModuleScore(
        object = epi_malignant_sce,
        features = list(epithelial_markers), 
        name = "epithelial_signature"
    )
    epi_malignant_sce <- AddModuleScore(
        object = epi_malignant_sce,
        features = list(malignant_markers), 
        name = "malignant_signature"
    )

    sig_frame <- epi_malignant_sce@meta.data[, c("epithelial_signature1","malignant_signature1","TorN","group_id","orig.ident")]

    gp <- sig_frame %>%
        ggplot(aes(x = epithelial_signature1, 
        y = malignant_signature1, 
        color = group_id)) +
        geom_point(size = 0.3) +
        scale_color_manual(values = group_cols) +
        theme_bw()

    ggsave(filename = "MalignantEpithelialScatterPlot.png", plot = gp, width = 8, height = 8, units = c("in"))
    ggsave(filename = "MalignantEpithelialScatterPlot.pdf", plot = gp, width = 8, height = 8, units = c("in"))

}








# Figure 9 Malignant and Epithelial differentially genes volcano Plot and GSEA analysis
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure9")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    markers <- readRDS("/seurat_jerry/hcc_project/P04_MalignantCell/P999_Malignant_Analysis/Figure8/malignant_vs_epithelial_markers.rds")

    # GSEA analysis

    geneList= markers$avg_log2FC 
    names(geneList)= toupper(rownames(markers))
    geneList=sort(geneList,decreasing = T)
    head(geneList)
    library(ggplot2)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    #选择gmt文件（MigDB中的全部基因集）
    gmtfile ='/seurat_jerry/hcc_project/P04_MalignantCell/P999_Malignant_Analysis/data/h.all.v6.2.symbols.gmt'
    # 31120 个基因集
    #GSEA分析
    library(GSEABase) # BiocManager::install('GSEABase')
    geneset <- read.gmt( gmtfile )  
    length(unique(geneset$term))
    egmt <- GSEA(geneList,
        TERM2GENE = geneset,
        minGSSize = 1,
        pvalueCutoff = 0.99,
        verbose = FALSE
    )

    head(egmt)
    egmt@result 
    gsea_results_df <- egmt@result 
    rownames(gsea_results_df)
    write.csv(gsea_results_df,file = 'gsea_results_df.csv')
    library(enrichplot)

    pdf("HALLMARK_OXIDATIVE_PHOSPHORYLATION.pdf",width = 12,height = 6)
    gseaplot2(egmt, geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", pvalue_table = F)
    dev.off()

    pdf("HALLMARK_ADIPOGENESIS.pdf",width = 12,height = 6) 
    gseaplot2(egmt,geneSetID = 'HALLMARK_ADIPOGENESIS',pvalue_table=F) 
    dev.off()


    pdf("HALLMARK_FATTY_ACID_METABOLISM.pdf",width = 12,height = 6) 
    gseaplot2(egmt,geneSetID = 'HALLMARK_FATTY_ACID_METABOLISM',pvalue_table=F) 
    dev.off()



    # gsea analysis
    mdb_h <- msigdbr(species = "Homo sapiens", category = "H")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    markers$genes = rownames(markers)
    cluster0.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
    ranks<- deframe(cluster0.genes)

    fgseaRes <- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

    png("TP53.png",width = 1080,height = 768)
    plotEnrichment(fgsea_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) +
        labs(title = "HALLMARK_ALLOGRAFT_REJECTION")
    dev.off()
    library(enrichplot)
    png("TP53.png",width = 1080,height = 768)
    gseaplot2(fgseaRes,geneSetID = 'HALLMARK_MTORC1_SIGNALING',pvalue_table=T)
    dev.off()


    gp <- ggplot(fgseaRes %>% 
        as_tibble() %>% 
        arrange(desc(NES)) %>% 
        filter(pval < 0.05) %>% 
        head(n= 20), aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = NES)) +
            coord_flip() +
            labs(x = "Hallmark", y = "Normalized Enrichment Score", 
            title = "HALLMARK gene sets NES from GSEA") +
            theme_bw()

    ggsave(filename = "gsea_barplot.png", width = 8, height = 6, units = c("in"))
    ggsave(filename = "gsea_barplot.pdf", width = 8, height = 6, units = c("in"))



    # Pseudo Volcalno Plot
    log2FC_threshold = 1
    p_value_threshold = 0.01
    markers[which(markers$p_val_adj < 0.01 & markers$avg_log2FC <= -log2FC_threshold),'sig'] <- 'Epithelial'
    markers[which(markers$p_val_adj < 0.01 & markers$avg_log2FC >= log2FC_threshold),'sig'] <- 'Malignant'
    markers[which(markers$p_val_adj >= 0.01 | abs(markers$avg_log2FC) < log2FC_threshold),'sig'] <- 'None'
    # 横轴 log2FC，纵轴 -log10(adj.P.Val)，颜色表示差异
    markers$Gene <- row.names(markers)
    markers$percent_difference <- markers$pct.1 - markers$pct.2
    p <- ggplot(markers, aes(x = percent_difference, y = avg_log2FC, color = sig)) +
        geom_point(alpha = 0.8, size = 0.6) +
        scale_colour_manual(values = c("red2", "blue2", "gray"), limits = c("Malignant", "Epithelial", "None")) +
        theme(panel.grid = element_blank(), 
            panel.background = element_rect(color = "black", fill = "transparent"), 
            plot.title = element_text(hjust = 0.5)) +
        theme(legend.key = element_rect(fill = "transparent"), 
            legend.background = element_rect(fill = "transparent"), 
            legend.position = c(0.9, 0.93)) +
        xlim(-2, 2) +
        ylim(-4, 4) +
        labs(
            x = "\nPercentage Difference",
            y = "Log2Fold Change\n", 
            color = "", 
            title = "Malignant vs Epithelial\n"
        )

    up <- subset(markers, sig == 'Malignant')
    up <- up[order(up$p_val_adj), ][1:40, ]
    down <- subset(markers, sig == 'Epithelial')
    down <- down[order(down$p_val_adj), ][1:40, ]

    p1 <- p + theme(legend.position = 'right') +
            geom_label_repel(data = rbind(up, down), 
                    aes(x = percent_difference, y = avg_log2FC, label = Gene),
                    size = 3,
                    box.padding = unit(0.5, 'lines'), 
                    segment.color = 'black',
                    alpha = 0.7,
                    max.overlaps = Inf, 
                    show.legend = FALSE)

    ggsave('Malignant_Epithelial_VolcanoPlot.png', p1, width = 6 * 2, height = 4 * 2)
    ggsave('Malignant_Epithelial_VolcanoPlot.pdf', p1, width = 6 * 2, height = 4 * 2)

}






# Figure10: Find All Markers among Samples
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure10")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    malignant_sce <- subset(global_sce, clusters %in% c("EPCAM+"))

    # malignant cell num of each sample should be >= 50
    tumor_sample_array <- names(table(malignant_sce@meta.data$orig.ident[grepl("T", malignant_sce@meta.data$orig.ident)])[table(malignant_sce@meta.data$orig.ident[grepl("T", malignant_sce@meta.data$orig.ident)]) > 50])
    malignant_sce <- subset(malignant_sce,orig.ident %in% tumor_sample_array)

    findmarker_list <- list()
    for (i in 1:length(tumor_sample_array)) {
        sample_id <- tumor_sample_array[i]
        markers <- FindMarkers(malignant_sce, group.by = "orig.ident", ident.1 = sample_id)
        markers$gene <- row.names(markers)
        markers$sample_id <- sample_id
        findmarker_list[[i]] <- markers
    }

    allmarkers_frame <- do.call(rbind,findmarker_list)

    AllMakers <- 'all_markers.csv'
    allmarkers_frame <- allmarkers_frame %>% group_by(sample_id)
    write.csv(allmarkers_frame, file=AllMakers, quote=F)

    top20 <- allmarkers_frame %>% group_by(sample_id) %>% top_n(n = 15, wt = avg_log2FC)

    write.csv(top20, file="top_marker.tsv", quote=F)

    #sce <- subset(sce,downsample = 500)
    Idents(malignant_sce) <- malignant_sce$orig.ident
    graph <- 'top_markers_heatmap.pdf'
    p1 <- DoHeatmap(malignant_sce, features = top20$gene, label = T)
    ggsave(graph, p1, width = 14, height = 40)

    graph <- 'top_markers_heatmap.png'
    ggsave(graph, p1, width = 14, height = 40,units = "in")


    # GSEA analysis

    AllMarkers <- read.csv(file = "all_markers.csv",row.names = 1)
    ZP04T_markers <- AllMarkers %>% filter(sample_id %in% c("ZP04T"))

    mdb_h <- msigdbr(species = "Homo sapiens", category = "H")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)
    cluster0.genes<- ZP04T_markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC)
    ranks<- deframe(cluster0.genes)

    fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

    png("TP53.png",width = 1080,height = 768)
    plotEnrichment(fgsea_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) +
        labs(title = "HALLMARK_ALLOGRAFT_REJECTION")
    dev.off()

    gp <- ggplot(fgseaRes %>% 
        as_tibble() %>% 
        arrange(desc(NES)) %>% 
        filter(pval < 0.05) %>% 
        head(n= 20), aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = NES)) +
            coord_flip() +
            labs(x = "Hallmark", y = "Normalized Enrichment Score", 
            title = "HALLMARK gene sets NES from GSEA") +
            theme_bw()

    ggsave(filename = "gsea_barplot.png", width = 8, height = 6, units = c("in"))
    ggsave(filename = "gsea_barplot.pdf", width = 8, height = 6, units = c("in"))



}





# Figure 11 CopyKat by Sample

# step 1 
if (F) {
   
    fig_outdir <- paste0(outdir, "/", "Figure11")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("EPCAM+")]
    
    group_id_array <- paste0(
        substr(malignant_sce@meta.data$orig.ident, 1, 2),
        substr(malignant_sce@meta.data$orig.ident, 5, 5)
    )
    malignant_sce <- AddMetaData(malignant_sce,
        group_id_array,
        col.name = "group_id"
    )

    malignant_sce <- subset(malignant_sce, group_id %in% c("ZPT", "ZRT"))

    counts <- as.matrix(malignant_sce@assays$RNA@counts)

    sample_array <- as.character(names(table(malignant_sce$orig.ident)[table(malignant_sce$orig.ident) >= 100]))


    for (sample_id in sample_array) {
        cat("Processing", sample_id, "\n")
        sub_malignant_sce <- malignant_sce[, malignant_sce$orig.ident %in% c(sample_id)]
        cat("Sample:", sample_id, nrow(sub_malignant_sce@meta.data),"\n")
        if (nrow(sub_malignant_sce@meta.data) > 100) {
            sample_dir <- file.path(fig_outdir,sample_id)
            if (!file.exists(sample_dir)) {
                dir.create(sample_dir)
            }
            setwd(sample_dir)

            merge_sce <- sub_malignant_sce

            counts <- as.matrix(merge_sce@assays$RNA@counts)
            write.table(counts, file = paste0(sample_id, ".merge_express.matrix.txt"), sep = "\t", quote = FALSE)

        }
    }

}


# step 2: CopyKat
if (F) {

    library(copykat)

    fig_outdir <- paste0(outdir, "/", "Figure11")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    for (sample_id in list.files(fig_outdir,"ZP02T")) {
        setwd(file.path(fig_outdir, sample_id))
        cat("Parsing", sample_id, "\n")
        exp.rawdata = read.table(
            file = paste0(sample_id, ".merge_express.matrix.txt"),
            sep = "\t", 
            header = T, 
            stringsAsFactors = F
        )
        copykat.test <- copykat(
            rawmat = exp.rawdata,
            id.type = "S", 
            ngene.chr = 5, 
            win.size = 25, 
            KS.cut = 0.3, 
            sam.name = sample_id, 
            distance = "euclidean", 
            norm.cell.names = "", 
            output.seg = "FLASE", 
            plot.genes = "TRUE", 
            n.cores = 28
        )
    }



}






# Figure 12 Demo NMF program inference


if (F) {

    step1 = function(dir_input = "count_data", 
        dir_output = "res1", 
        k = 3:10, 
        iteration=200, 
        feature=2000){
        
        library(tidyverse)

        files=dir(dir_input)
        for (i in files) {
            
            one.file.path=paste(dir_input,"/",i,sep = "")
            prefix=strsplit(i,"\\.")[[1]][1]
            k.choices=paste(k,collapse = " ")
            junk.file=paste(dir_output,"/",prefix,"/","cnmf_tmp/",prefix,".spectra.k_*.iter_*.df.npz",sep = "")
            
            system(paste("/mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py prepare --output-dir ",
                dir_output, " --name ", prefix, " -c ", one.file.path, " -k ", k.choices, " --n-iter ", iteration, " --total-workers 1 --numgenes ", feature,
                sep = ""
            ))
            system(paste("/mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py factorize --output-dir ",
                dir_output, " --name ", prefix, " --worker-index 0",
                sep = ""
            ))
            system(paste("/mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py combine --output-dir ",
                dir_output, " --name ", prefix,
                sep = ""
            ))
            system(paste("rm -f ",junk.file,sep = ""))
            system(paste("MPLBACKEND='Agg' /mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py k_selection_plot --output-dir ",
                dir_output, " --name ", prefix,
                sep = ""
            ))
        }
        
        dir.create( paste(dir_output,"/k_selection",sep = ""))
        file.create(paste(dir_output,"/k_selection.txt",sep = ""))
        path3=paste(dir_output,"/k_selection.txt",sep = "")
        
        for (i in str_replace(files,"\\..*$","")) {
            path1=paste(dir_output,"/",i,"/",i,".k_selection.png",sep = "")
            path2=paste(dir_output,"/k_selection/",i,".k_selection.png",sep = "")
            system(paste("cp ",path1," ",path2,sep = ""))
            
            cat(i,",\n",sep = "",file=path3,append = T)
        }
    
    }


    step2 = function(dir_input="res1",
        dir_output="res2", 
        dir_count="count_data", 
        usage_filter = 0.03, 
        top_gene = 30, 
        cor_min = 0, 
        cor_max = 0.6, 
        color = NULL,
        cluster_method = "complete",
        scale_min = -2,
        scale_max = 2,
        cluster_method2 = "complete"){
        
        library(tidyverse)
        
        dir.create(dir_output)
        dirs=setdiff(dir(dir_input),c("k_selection","k_selection.txt"))
        ref.file=read.table(paste(dir_input,"/k_selection.txt",sep = ""),header = F,sep = ",",stringsAsFactors = F)
        colnames(ref.file)=c("sample","k")
        rownames(ref.file)=ref.file$sample
        
        for (i in dirs) {
            cat("Parsing consensus ",i,"\n")
            system(paste("MPLBACKEND='Agg' /mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py consensus --output-dir ",dir_input," --name ",i," --components ",ref.file[i,"k"]," --local-density-threshold 0.02 --show-clustering",sep = ""))
            
            path1=paste(dir_input,"/",i,"/",i,".usages.k_",ref.file[i,"k"],".dt_0_02.consensus.txt",sep = "")
            path2=paste(dir_input,"/",i,"/",i,".gene_spectra_score.k_",ref.file[i,"k"],".dt_0_02.txt",sep = "")
            
            if( file.exists(path1) & file.exists(path2) ){
            system(paste("cp ",path1," ",path2," ",dir_output,"/",sep = ""))
            }
        }
        
        ###################################################################################
        for (i in dirs) {
            cat("Quality control ", i, "\n")
            usage.file=dir(dir_output,pattern = paste(i,".usages",sep = ""))
            usage.df=read.table(paste(dir_output,"/",usage.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
            colnames(usage.df)=paste(i,1:dim(usage.df)[2],sep = "_")
            
            #normalize
            usage.df=usage.df / rowSums(usage.df)
            write.table(usage.df,file = paste(dir_output,"/",i,"_program.usage.norm.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
            
            #QC1
            tmpdf1=gather(usage.df,"program","ratio")
            tmpdf1%>%ggplot(aes(x=program,y=ratio))+geom_boxplot(outlier.shape = NA)+geom_jitter(color="red",alpha=0.4,width = 0.2)+
            labs(title = i)+
            theme(
                axis.text.x.bottom = element_text(angle = 45,hjust = 1),
                plot.title = element_text(hjust = 0.5,size=20)
            )
            ggsave(paste(dir_output,"/",i,"_program.usage.norm.QC.png",sep = ""),device = "png",width = 20,height = 16,units = c("cm"))
            
            #score
            score.file=dir(dir_output,pattern = paste(i,".gene_spectra_score",sep = ""))
            score.df=read.table(paste(dir_output,"/",score.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
            score.df=as.data.frame(t(score.df))
            colnames(score.df)=paste(i,1:dim(score.df)[2],sep = "_")
            
            topn.df=as.data.frame(matrix(nrow = top_gene,ncol = ncol(score.df)))
            colnames(topn.df)=colnames(score.df)
            
            for (k in colnames(score.df)) {
            tmpv=score.df[,k]
            names(tmpv)=rownames(score.df)
            topn.df[,k]=names(rev(tail(sort(tmpv),top_gene)))
            }
            
            #save
            write.table(topn.df, file = paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
            score.df$gene=rownames(score.df)
            write.table(score.df,file = paste(dir_output,"/",i,"_program.Zscore.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
        }
        
        ###################################################################################
        check.usage=data.frame()
        
        for (i in dirs) {
            usage.file = paste(dir_output,"/",i,"_program.usage.norm.txt",sep = "")
            usage.df = read.table(usage.file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
            check.usage = rbind(check.usage,as.data.frame(colMeans(usage.df)))
        }
        colnames(check.usage)=c("mean_ratio")
        
        check.usage$sample_programs=rownames(check.usage)
        check.usage=check.usage%>%arrange(mean_ratio)
        check.usage$sample_programs=factor(check.usage$sample_programs,levels = check.usage$sample_programs)
        
        linex=sum(check.usage$mean_ratio < usage_filter)
        check.usage%>%ggplot(aes(x=sample_programs,y=mean_ratio))+geom_point()+
            geom_hline(yintercept = usage_filter,color="red")+
            geom_vline(xintercept = linex+0.5,color="red")+
            theme(
            axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)
            )
        ggsave(paste(dir_output,"/","check.usage.png",sep = ""),width = 30,height = 16,device = "png",units = "cm")
        
        maybe.bg=as.character(check.usage$sample_programs[check.usage$mean_ratio < usage_filter])
        
        ###################################################################################
        library(pheatmap)
        library(RColorBrewer)
        library(scales)
        
        all.score.df = data.frame()
        all.score.topn.df = data.frame()
        
        for (i in dirs) {
            score.file=paste(dir_output,"/",i,"_program.Zscore.txt",sep = "")
            score.df=read.table(score.file,header = T,sep = "\t",stringsAsFactors = F)
            if (i==dirs[1]) {all.score.df=score.df}
            if (i!=dirs[1]) {
            all.score.df=all.score.df%>%inner_join(score.df,by="gene")
            }
            
            score.topn.file=paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = "")
            score.topn.df=read.table(score.topn.file,header = T,sep = "\t",stringsAsFactors = F)
            if (i==dirs[1]) {all.score.topn.df=score.topn.df}
            if (i!=dirs[1]) {
            all.score.topn.df=cbind(all.score.topn.df,score.topn.df)
            }
        }
        
        rownames(all.score.df)=all.score.df$gene
        all.score.df$gene=NULL
        all.score.df=all.score.df[rowSums(is.na(all.score.df)) == 0,] #可能有空值，需要去掉
        all.score.rm.df=all.score.df[,setdiff(colnames(all.score.df),maybe.bg)] #在质控这一步检测出来的噪声
        all.score.rm.df.cor=cor(all.score.rm.df,method = "pearson")
        
        all.score.rm.df.cor[all.score.rm.df.cor < cor_min]=cor_min
        all.score.rm.df.cor[all.score.rm.df.cor > cor_max]=cor_max
        
        colanno=as.data.frame(colnames(all.score.rm.df.cor))
        colnames(colanno)="colnames"
        colanno$sample=str_replace(colanno$colnames,"_.*","")
        rownames(colanno)=colanno$colnames
        colanno$colnames=NULL
        
        rowanno=as.data.frame(rownames(all.score.rm.df.cor))
        colnames(rowanno)="rownames"
        rowanno$sample=str_replace(rowanno$rownames,"_.*","")
        rownames(rowanno)=rowanno$rownames
        rowanno$rownames=NULL
        
        #指定注释条的颜色
        if (is.null(color)){
            color_v=colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(dirs))
            names(color_v)=dirs
        }else{
            color_v=color
        }
        ann_colors = list(sample = color_v)
        
        #画图
        tmpp=pheatmap(all.score.rm.df.cor,cluster_rows = T,cluster_cols = T,
                        clustering_method = cluster_method, 
                        show_colnames = F,
                        treeheight_row=30,treeheight_col=0,
                        border_color=NA,
                        annotation_row = rowanno,annotation_col = colanno,
                        annotation_names_row = F,annotation_names_col = F,
                        annotation_colors = ann_colors,
                        color = colorRampPalette(c("white","yellow", "red","#67001F"))(50),
                        fontsize_row=12,
                        width = 11.5,height = 9,
                        filename = paste(dir_output,"/","program_pearson_cor.",cluster_method,".heatmap.pdf",sep = "")
        )
        
        #保存画图数据
        write.table(all.score.rm.df.cor,file = paste(dir_output,"/","cor_heatmap_data.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
        #保存program和(最相关)gene的对应关系
        all.score.topn.rm.df=all.score.topn.df[,setdiff(colnames(all.score.topn.df),maybe.bg)]#在质控这一步检测出来的噪声
        write.table(all.score.topn.rm.df,file = paste(dir_output,"/","program_topngene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
        
        ###################################################################################
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(xlsx)
        
        hsets <- read.gmt("hallmark_cancersea.gmt")
        enrich.result=data.frame()
        pathway_v=c()
        program_v=c()
        
        program_topn=read.table(paste(dir_output,"/","program_topngene.txt",sep = ""),header = T,sep = "\t",stringsAsFactors = F)
        for (i in 1:dim(program_topn)[2]) {
            tmp <- enricher(program_topn[,i], TERM2GENE = hsets)
            if (is.null(tmp)) {
            next
            }
            
            tmp1=head(tmp@result)
            tmp1$program=colnames(program_topn)[i]
            rownames(tmp1)=NULL
            enrich.result=rbind(enrich.result,tmp1)
            
            program_v=append(program_v,colnames(program_topn)[i])
            pathway_v=append(pathway_v,paste(tmp1$Description,collapse = ","))
        }
        
        write.xlsx(enrich.result,file = paste(dir_output,"/","program_topngene_enrichment.xlsx",sep = ""),row.names = F)
        
        enrich.df=data.frame(program=program_v,pathway=pathway_v)
        enrich.df$program=factor(enrich.df$program,levels = tmpp$tree_row$labels[tmpp$tree_row$order])
        enrich.df=enrich.df%>%arrange(program)
        write.csv(enrich.df,file = paste(dir_output,"/","program_topngene_enrichment_order.csv",sep = ""),row.names = F,quote = F)
        
        ###################################################################################
        for (i in dirs) {
            
            one_matrix=read.table(paste(dir_count,"/",i,".count.txt",sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
            RowSum=rowSums(one_matrix)
            one_matrix=log1p((one_matrix / RowSum) * 10000)
            one_matrix=as.data.frame(t(one_matrix))
            
            plot_gene=read.table(paste(dir_output,"/","program_topngene.txt",sep = ""),header = T,sep = "\t",stringsAsFactors = F)
            plot_gene=plot_gene[,str_detect(colnames(plot_gene),i)]
            plot_gene=gather(plot_gene,key = "program",value = "gene")
            plot_gene=plot_gene%>%arrange(program)
            
            one_matrix=one_matrix[plot_gene$gene,]
            one_matrix=t(scale(t(one_matrix)))
            one_matrix[one_matrix < scale_min] = scale_min
            one_matrix[one_matrix > scale_max] = scale_max
            
            if(length(unique(plot_gene$gene)) < length(plot_gene$gene)){
            plot_gene$gene=rownames(one_matrix)
            }
            
            pn=length(unique(plot_gene$program))
            tmpp2=pheatmap(one_matrix,cluster_rows = F,cluster_cols = T,
                        color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
                        treeheight_col=0,
                        clustering_method = cluster_method2,
                        show_colnames=F,
                        border_color=NA,
                        gaps_row=as.numeric(cumsum(table(plot_gene$program))[- pn])
            )
            write.table(one_matrix,paste(dir_output,"/",i,"_data_heatmap.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
            
            cell_sort=tmpp2$tree_col$labels[tmpp2$tree_col$order]
            gene_sort=plot_gene$gene
            matrix_new=as.data.frame(one_matrix)
            matrix_new$gene=rownames(matrix_new)
            matrix_new=matrix_new%>%reshape2::melt(id="gene")
            colnames(matrix_new)[c(2,3)]=c("cell","exp")
            matrix_new$gene=factor(matrix_new$gene,levels = rev(gene_sort))
            matrix_new$cell=factor(matrix_new$cell,levels = cell_sort)
            
            plot1=matrix_new%>%ggplot(aes(x=cell,y=gene,fill=exp))+geom_tile()+
            geom_hline(yintercept = as.numeric(cumsum(table(plot_gene$program))[- pn])+0.5,color="black",linetype=5)+
            labs(title = paste(i,": ",length(cell_sort)," cells; ",pn," programs",sep = ""))+
            scale_x_discrete("",expand = c(0,0))+
            scale_y_discrete("",expand = c(0,0))+
            scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100))+
            theme_bw()+
            theme(
                panel.grid = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x.bottom = element_blank(),axis.text.y.left = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 20),
                legend.position = "left"
            )
            
            dim1=5*pn
            dim2=(top_gene / 10) *2
            gene.text=t(matrix(plot_gene$gene,nrow = dim2,ncol = dim1))
            rownames(gene.text)= seq(-0.5,-(dim1-0.5),-1)
            colnames(gene.text)=seq(0.5,(dim2-0.5),1)
            gene.text=as.data.frame(gene.text)
            gene.text$dim1=rownames(gene.text)
            gene.text=reshape2::melt(gene.text,id="dim1")
            colnames(gene.text)[2:3]=c("dim2","gene")
            gene.text$dim1=as.numeric(as.character(gene.text$dim1))
            gene.text$dim2=as.numeric(as.character(gene.text$dim2))
            
            plot2=gene.text%>%ggplot(aes(x=dim2,y=dim1))+geom_text(aes(label=gene))+
            geom_hline(yintercept = seq(-dim1,0,5)[-c(1,length(seq(-dim1,0,5)))],color="black",linetype=5)+
            scale_x_continuous("",expand = c(0,0),limits = c(0,10))+
            scale_y_continuous("",expand = c(0,0),limits = c(-dim1,0))+
            labs(title = paste(unique(plot_gene$program),collapse = "; "))+
            theme_bw()+
            theme(
                panel.grid = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 20)
            )
            
            library(patchwork)
            plot3=plot1+plot2+plot_layout(widths = c(1,2))
            ggsave(filename = paste(dir_output,"/",i,"_program_gene.heatmap.pdf",sep = ""),plot = plot3,width = 46,height = 16,units = "cm")
        }
    }

    setwd("/seurat_jerry/hcc_project/P04_MalignantCell/P999_Malignant_Analysis/Figure12")
    # step 3 use funtion step1 and step2
    library(reticulate)

    ### 调用子环境的python
    use_condaenv(condaenv = "cnmf_env", required = T,conda = "/mnt/share01/tools/bin/conda")
    py_config() #如果显示cnmf_env环境里面的python就OK

    ### 第一步 ################
    # 表达矩阵文件命名尽量和我的示例文件一致(HNSCC17.count.txt)，保证第一个点(.)前面是样本名称/ID。这些文件都保存在某个文件夹里面，比如示例中的"count_data"
    # source("1.R")

    step1(
        dir_input = "count_data",
        dir_output = "res1", 
        k = 3:5, 
        iteration = 50
    ) # 这里为了演示方便，取值都比较小

    # ###
    # 到这儿先停一停，step1会在dir_output下面为每一个样本生成一个文件夹。
    # 此外还有一个k_selection文件夹以及一个k_selection.txt文件。
    # k_selection文件夹里面有每一个样本选k值的图片（点击图片可以在浏览器打开），以此为依据，将k_selection.txt补充完整（可以直接在Rstudio打开）。
    # k_selection.txt有两列，以逗号分割，第二列是空的，需要你填上去。
    # ###


    ### 第二步 ################
    # source("2.R")

    step2(
        dir_input = "res1",
        dir_output = "res2", 
        dir_count = "count_data", 
        usage_filter = 0.03, 
        top_gene = 30, 
        cor_min = 0, 
        cor_max = 0.6
    )

    #过滤阈值，10X 0.01 smart-seq 0.03，这是我用过的阈值，不一定适用于所有情况
    #cor_min可以设为0，尽管真实范围会小于0；因为原文用的是0.6，所以我这里也限制最大相关系数为0.6，主要是为了图明显一些。试过10X数据，相关系数比较高(dropout导致的假象)，所以没有限制
    #top_gene按照20 30 50来设计的
    #有一个color参数，可以给样本涂色，格式是命名后的字符串向量；函数里面有一个默认的配色，也还OK。

    ### 实际分析中的一点点经验 ################
    #1. 相关性聚类这一步需要多做几次，所以原始画图的数据我都导出了，方便你自己改图
    #根据每个program的大致功能，如果发现另一种功能的program聚到某一种meta模块里面，这时可以将乱入的program删掉，再做一次相关性聚类。
    #如果我们认定应该属于同一个meta模块的program分散在两个地方，可以试试调整聚类方法(参数clustering_method)，或者像"_1""_2"这样定义
    #2. 对top基因做注释，尽量选择和研究背景有关联的基因集，比如研究肿瘤细胞异质性时，选择和肿瘤相关的基因集


    library(tidyverse)

    mat.files=dir("./res2/",pattern = "dt_0_02.txt$")

    all.mat=data.frame()
    for (fi in mat.files) {
    tmp.mat=read.table(paste0("./res2/",fi),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    tmp.mat=as.data.frame(t(tmp.mat))
    sampleid=str_replace(fi,"\\..*$","")
    colnames(tmp.mat)=paste(sampleid,colnames(tmp.mat),sep = "_")
    tmp.mat$gene=rownames(tmp.mat)
    
    if (sampleid == "HNSCC22") {
        all.mat=tmp.mat
    }else{
        all.mat=all.mat %>% full_join(tmp.mat,by="gene") #元素的并集进行合并
    }
    }

    # 对于某一个模块
    signature.programs=c("HNSCC6_3","HNSCC22_4","HNSCC5_3")
    signature.loading=all.mat[,c("gene",signature.programs)]

    used.gene=c()
    for (pi in signature.programs) {
    tmp.df=signature.loading[,c("gene",pi)]
    tmp.loading=tmp.df[,2]
    names(tmp.loading)=tmp.df[,1]
    
    tmp.loading=tmp.loading[!is.na(tmp.loading)]
    used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
    }
    used.gene=unique(used.gene)

    signature.loading=signature.loading[signature.loading$gene %in% used.gene,]
    rownames(signature.loading)=signature.loading$gene
    signature.loading$gene=NULL
    signature.loading[is.na(signature.loading)]<-0
    signature.loading$total_loading=rowSums(signature.loading)
    signature.loading$average_loading=signature.loading$total_loading / length(signature.programs)

    signature.loading=signature.loading%>%arrange(desc(average_loading))
    head(rownames(signature.loading),30)

    # library(Seurat)
    # 
    # tmp.gene=c("TUBA1B", "HMGB2",  "TOP2A",  "UBE2C",  "NUSAP1", "MKI67",  "PBK",
    #            "BIRC5",  "CDK1",   "AURKB",  "UBE2T",  "CKS1B",  "H2AFZ",  "TPX2",   
    #            "TK1",    "CCNA2",  "GTSE1", "CEP55",  "KIF23",  "CENPF",  "CKS2",
    #            "HMGB1",  "PTTG1",  "CDCA3",  "CDKN3",  "PRC1",   "NUF2",   "CCNB1",
    #            "CCNB2",  "FOXM1")
    # sum(tmp.gene %in% c(Seurat::cc.genes.updated.2019$s.genes,Seurat::cc.genes.updated.2019$g2m.genes))


























}






# Figure 13 HCC NMF program inference
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure13")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)
    malignant_sce <- global_sce[, global_sce@meta.data$clusters %in% c("EPCAM+")]
    
    group_id_array <- paste0(
        substr(malignant_sce@meta.data$orig.ident, 1, 2),
        substr(malignant_sce@meta.data$orig.ident, 5, 5)
    )
    malignant_sce <- AddMetaData(malignant_sce,
        group_id_array,
        col.name = "group_id"
    )

    malignant_sce <- subset(malignant_sce, group_id %in% c("ZPT", "ZRT"))

    counts <- as.matrix(malignant_sce@assays$RNA@counts)


    sample_array <- as.character(names(table(malignant_sce$orig.ident)[table(malignant_sce$orig.ident) >= 100]))

    count_dir <- file.path(fig_outdir,"count_data")
    if (!file.exists(count_dir)) {
        dir.create(count_dir)
    }

    for (sample_id in sample_array) {
        cat("Processing", sample_id, "\n")
        sub_malignant_sce <- malignant_sce[, malignant_sce$orig.ident %in% c(sample_id)]
        cat("Sample:", sample_id, nrow(sub_malignant_sce@meta.data),"\n")

        if (nrow(sub_malignant_sce@meta.data) > 30) {
            
            sce_merge <- sub_malignant_sce

            count <- as.data.frame(sce_merge@assays$RNA@counts)
            count=count[rowSums(count) > 0,]
            count = count[!str_detect(row.names(count), "^MT-"),]

            # transpose the matrix
            cat("Transposing matrix...\n")
            count = data.frame(t(count))
            cat("Transposing complete.\n")
            
            cat("Output matrix to file...\n")
            write.table(count, 
                        file = file.path(count_dir,paste0(sample_id,".count.txt")),
                        quote = F,
                        sep = "\t",
                        row.names = T,
                        col.names = T) 
        }
    }


    step1 = function(dir_input = "count_data", 
        dir_output = "res1", 
        k = 3:10, 
        iteration=200, 
        feature=2000){
        
        library(tidyverse)

        files=dir(dir_input)
        for (i in files) {
            
            one.file.path=paste(dir_input,"/",i,sep = "")
            prefix=strsplit(i,"\\.")[[1]][1]
            k.choices=paste(k,collapse = " ")
            junk.file=paste(dir_output,"/",prefix,"/","cnmf_tmp/",prefix,".spectra.k_*.iter_*.df.npz",sep = "")
            
            system(paste("/mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py prepare --output-dir ",
                dir_output, " --name ", prefix, " -c ", one.file.path, " -k ", k.choices, " --n-iter ", iteration, " --total-workers 1 --numgenes ", feature,
                sep = ""
            ))
            system(paste("/mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py factorize --output-dir ",
                dir_output, " --name ", prefix, " --worker-index 0",
                sep = ""
            ))
            system(paste("/mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py combine --output-dir ",
                dir_output, " --name ", prefix,
                sep = ""
            ))
            system(paste("rm -f ",junk.file,sep = ""))
            system(paste("MPLBACKEND='Agg' /mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py k_selection_plot --output-dir ",
                dir_output, " --name ", prefix,
                sep = ""
            ))
        }
        
        dir.create( paste(dir_output,"/k_selection",sep = ""))
        file.create(paste(dir_output,"/k_selection.txt",sep = ""))
        path3=paste(dir_output,"/k_selection.txt",sep = "")
        
        for (i in str_replace(files,"\\..*$","")) {
            path1=paste(dir_output,"/",i,"/",i,".k_selection.png",sep = "")
            path2=paste(dir_output,"/k_selection/",i,".k_selection.png",sep = "")
            system(paste("cp ",path1," ",path2,sep = ""))
            
            cat(i,",\n",sep = "",file=path3,append = T)
        }
    
    }


    step2 = function(dir_input="res1",
        dir_output="res2", 
        dir_count="count_data", 
        usage_filter = 0.03, 
        top_gene = 30, 
        cor_min = 0, 
        cor_max = 0.6, 
        color = NULL,
        cluster_method = "complete",
        scale_min = -2,
        scale_max = 2,
        cluster_method2 = "complete"){
        
        library(tidyverse)
        
        dir.create(dir_output)
        dirs=setdiff(dir(dir_input),c("k_selection","k_selection.txt"))
        ref.file=read.table(paste(dir_input,"/k_selection.txt",sep = ""),header = F,sep = ",",stringsAsFactors = F)
        colnames(ref.file)=c("sample","k")
        rownames(ref.file)=ref.file$sample
        
        for (i in dirs) {
            cat("Parsing consensus ",i,"\n")
            system(paste("MPLBACKEND='Agg' /mnt/share01/tools/miniconda/envs/cnmf_env/bin/python -W ignore cnmf.py consensus --output-dir ",dir_input," --name ",i," --components ",ref.file[i,"k"]," --local-density-threshold 0.02 --show-clustering",sep = ""))
            
            path1=paste(dir_input,"/",i,"/",i,".usages.k_",ref.file[i,"k"],".dt_0_02.consensus.txt",sep = "")
            path2=paste(dir_input,"/",i,"/",i,".gene_spectra_score.k_",ref.file[i,"k"],".dt_0_02.txt",sep = "")
            
            if( file.exists(path1) & file.exists(path2) ){
            system(paste("cp ",path1," ",path2," ",dir_output,"/",sep = ""))
            }
        }
        
        ###################################################################################
        for (i in dirs) {
            cat("Quality control ", i, "\n")
            usage.file=dir(dir_output,pattern = paste(i,".usages",sep = ""))
            usage.df=read.table(paste(dir_output,"/",usage.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
            colnames(usage.df)=paste(i,1:dim(usage.df)[2],sep = "_")
            
            #normalize
            usage.df=usage.df / rowSums(usage.df)
            write.table(usage.df,file = paste(dir_output,"/",i,"_program.usage.norm.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
            
            #QC1
            tmpdf1=gather(usage.df,"program","ratio")
            tmpdf1%>%ggplot(aes(x=program,y=ratio))+geom_boxplot(outlier.shape = NA)+geom_jitter(color="red",alpha=0.4,width = 0.2)+
            labs(title = i)+
            theme(
                axis.text.x.bottom = element_text(angle = 45,hjust = 1),
                plot.title = element_text(hjust = 0.5,size=20)
            )
            ggsave(paste(dir_output,"/",i,"_program.usage.norm.QC.png",sep = ""),device = "png",width = 20,height = 16,units = c("cm"))
            
            #score
            score.file=dir(dir_output,pattern = paste(i,".gene_spectra_score",sep = ""))
            score.df=read.table(paste(dir_output,"/",score.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
            score.df=as.data.frame(t(score.df))
            colnames(score.df)=paste(i,1:dim(score.df)[2],sep = "_")
            
            topn.df=as.data.frame(matrix(nrow = top_gene,ncol = ncol(score.df)))
            colnames(topn.df)=colnames(score.df)
            
            for (k in colnames(score.df)) {
            tmpv=score.df[,k]
            names(tmpv)=rownames(score.df)
            topn.df[,k]=names(rev(tail(sort(tmpv),top_gene)))
            }
            
            #save
            write.table(topn.df, file = paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
            score.df$gene=rownames(score.df)
            write.table(score.df,file = paste(dir_output,"/",i,"_program.Zscore.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
        }
        
        ###################################################################################
        check.usage=data.frame()
        
        for (i in dirs) {
            usage.file = paste(dir_output,"/",i,"_program.usage.norm.txt",sep = "")
            usage.df = read.table(usage.file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
            check.usage = rbind(check.usage,as.data.frame(colMeans(usage.df)))
        }
        colnames(check.usage)=c("mean_ratio")
        
        check.usage$sample_programs=rownames(check.usage)
        check.usage=check.usage%>%arrange(mean_ratio)
        check.usage$sample_programs=factor(check.usage$sample_programs,levels = check.usage$sample_programs)
        
        linex=sum(check.usage$mean_ratio < usage_filter)
        check.usage%>%ggplot(aes(x=sample_programs,y=mean_ratio))+geom_point()+
            geom_hline(yintercept = usage_filter,color="red")+
            geom_vline(xintercept = linex+0.5,color="red")+
            theme(
            axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)
            )
        ggsave(paste(dir_output,"/","check.usage.png",sep = ""),width = 30,height = 16,device = "png",units = "cm")
        
        maybe.bg=as.character(check.usage$sample_programs[check.usage$mean_ratio < usage_filter])
        
        ###################################################################################
        library(pheatmap)
        library(RColorBrewer)
        library(scales)
        
        all.score.df = data.frame()
        all.score.topn.df = data.frame()
        
        for (i in dirs) {
            score.file=paste(dir_output,"/",i,"_program.Zscore.txt",sep = "")
            score.df=read.table(score.file,header = T,sep = "\t",stringsAsFactors = F)
            if (i==dirs[1]) {all.score.df=score.df}
            if (i!=dirs[1]) {
            all.score.df=all.score.df%>%inner_join(score.df,by="gene")
            }
            
            score.topn.file=paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = "")
            score.topn.df=read.table(score.topn.file,header = T,sep = "\t",stringsAsFactors = F)
            if (i==dirs[1]) {all.score.topn.df=score.topn.df}
            if (i!=dirs[1]) {
            all.score.topn.df=cbind(all.score.topn.df,score.topn.df)
            }
        }
        
        rownames(all.score.df)=all.score.df$gene
        all.score.df$gene=NULL
        all.score.df=all.score.df[rowSums(is.na(all.score.df)) == 0,] #可能有空值，需要去掉
        all.score.rm.df=all.score.df[,setdiff(colnames(all.score.df),maybe.bg)] #在质控这一步检测出来的噪声
        all.score.rm.df.cor=cor(all.score.rm.df,method = "pearson")
        
        all.score.rm.df.cor[all.score.rm.df.cor < cor_min]=cor_min
        all.score.rm.df.cor[all.score.rm.df.cor > cor_max]=cor_max
        
        colanno=as.data.frame(colnames(all.score.rm.df.cor))
        colnames(colanno)="colnames"
        colanno$sample=str_replace(colanno$colnames,"_.*","")
        rownames(colanno)=colanno$colnames
        colanno$colnames=NULL
        
        rowanno=as.data.frame(rownames(all.score.rm.df.cor))
        colnames(rowanno)="rownames"
        rowanno$sample=str_replace(rowanno$rownames,"_.*","")
        rownames(rowanno)=rowanno$rownames
        rowanno$rownames=NULL
        
        #指定注释条的颜色
        if (is.null(color)){
            color_v=colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(dirs))
            names(color_v)=dirs
        }else{
            color_v=color
        }
        ann_colors = list(sample = color_v)
        
        #画图
        tmpp=pheatmap(all.score.rm.df.cor,cluster_rows = T,cluster_cols = T,
                        clustering_method = cluster_method, 
                        show_colnames = F,
                        treeheight_row=30,treeheight_col=0,
                        border_color=NA,
                        annotation_row = rowanno,annotation_col = colanno,
                        annotation_names_row = F,annotation_names_col = F,
                        annotation_colors = ann_colors,
                        color = colorRampPalette(c("white","yellow", "red","#67001F"))(50),
                        fontsize_row=12,
                        width = 11.5,height = 9,
                        filename = paste(dir_output,"/","program_pearson_cor.",cluster_method,".heatmap.pdf",sep = "")
        )
        
        #保存画图数据
        write.table(all.score.rm.df.cor,file = paste(dir_output,"/","cor_heatmap_data.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
        #保存program和(最相关)gene的对应关系
        all.score.topn.rm.df=all.score.topn.df[,setdiff(colnames(all.score.topn.df),maybe.bg)]#在质控这一步检测出来的噪声
        write.table(all.score.topn.rm.df,file = paste(dir_output,"/","program_topngene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
        
        ###################################################################################
        library(clusterProfiler)
        library(org.Hs.eg.db)
        library(xlsx)
        
        hsets <- read.gmt("hallmark_cancersea.gmt")
        enrich.result=data.frame()
        pathway_v=c()
        program_v=c()
        
        program_topn=read.table(paste(dir_output,"/","program_topngene.txt",sep = ""),header = T,sep = "\t",stringsAsFactors = F)
        for (i in 1:dim(program_topn)[2]) {
            tmp <- enricher(program_topn[,i], TERM2GENE = hsets)
            if (is.null(tmp)) {
            next
            }
            
            tmp1=head(tmp@result)
            tmp1$program=colnames(program_topn)[i]
            rownames(tmp1)=NULL
            enrich.result=rbind(enrich.result,tmp1)
            
            program_v=append(program_v,colnames(program_topn)[i])
            pathway_v=append(pathway_v,paste(tmp1$Description,collapse = ","))
        }
        
        write.xlsx(enrich.result,file = paste(dir_output,"/","program_topngene_enrichment.xlsx",sep = ""),row.names = F)
        
        enrich.df=data.frame(program=program_v,pathway=pathway_v)
        enrich.df$program=factor(enrich.df$program,levels = tmpp$tree_row$labels[tmpp$tree_row$order])
        enrich.df=enrich.df%>%arrange(program)
        write.csv(enrich.df,file = paste(dir_output,"/","program_topngene_enrichment_order.csv",sep = ""),row.names = F,quote = F)
        
        ###################################################################################
        for (i in dirs) {
            
            one_matrix=read.table(paste(dir_count,"/",i,".count.txt",sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
            RowSum=rowSums(one_matrix)
            one_matrix=log1p((one_matrix / RowSum) * 10000)
            one_matrix=as.data.frame(t(one_matrix))
            
            plot_gene=read.table(paste(dir_output,"/","program_topngene.txt",sep = ""),header = T,sep = "\t",stringsAsFactors = F)
            plot_gene=plot_gene[,str_detect(colnames(plot_gene),i)]
            plot_gene=gather(plot_gene,key = "program",value = "gene")
            plot_gene=plot_gene%>%arrange(program)
            
            one_matrix=one_matrix[plot_gene$gene,]
            one_matrix=t(scale(t(one_matrix)))
            one_matrix[one_matrix < scale_min] = scale_min
            one_matrix[one_matrix > scale_max] = scale_max
            
            if(length(unique(plot_gene$gene)) < length(plot_gene$gene)){
            plot_gene$gene=rownames(one_matrix)
            }
            
            pn=length(unique(plot_gene$program))
            tmpp2=pheatmap(one_matrix,cluster_rows = F,cluster_cols = T,
                        color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
                        treeheight_col=0,
                        clustering_method = cluster_method2,
                        show_colnames=F,
                        border_color=NA,
                        gaps_row=as.numeric(cumsum(table(plot_gene$program))[- pn])
            )
            write.table(one_matrix,paste(dir_output,"/",i,"_data_heatmap.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
            
            cell_sort=tmpp2$tree_col$labels[tmpp2$tree_col$order]
            gene_sort=plot_gene$gene
            matrix_new=as.data.frame(one_matrix)
            matrix_new$gene=rownames(matrix_new)
            matrix_new=matrix_new%>%reshape2::melt(id="gene")
            colnames(matrix_new)[c(2,3)]=c("cell","exp")
            matrix_new$gene=factor(matrix_new$gene,levels = rev(gene_sort))
            matrix_new$cell=factor(matrix_new$cell,levels = cell_sort)
            
            plot1=matrix_new%>%ggplot(aes(x=cell,y=gene,fill=exp))+geom_tile()+
            geom_hline(yintercept = as.numeric(cumsum(table(plot_gene$program))[- pn])+0.5,color="black",linetype=5)+
            labs(title = paste(i,": ",length(cell_sort)," cells; ",pn," programs",sep = ""))+
            scale_x_discrete("",expand = c(0,0))+
            scale_y_discrete("",expand = c(0,0))+
            scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100))+
            theme_bw()+
            theme(
                panel.grid = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x.bottom = element_blank(),axis.text.y.left = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 20),
                legend.position = "left"
            )
            
            dim1=5*pn
            dim2=(top_gene / 10) *2
            gene.text=t(matrix(plot_gene$gene,nrow = dim2,ncol = dim1))
            rownames(gene.text)= seq(-0.5,-(dim1-0.5),-1)
            colnames(gene.text)=seq(0.5,(dim2-0.5),1)
            gene.text=as.data.frame(gene.text)
            gene.text$dim1=rownames(gene.text)
            gene.text=reshape2::melt(gene.text,id="dim1")
            colnames(gene.text)[2:3]=c("dim2","gene")
            gene.text$dim1=as.numeric(as.character(gene.text$dim1))
            gene.text$dim2=as.numeric(as.character(gene.text$dim2))
            
            plot2=gene.text%>%ggplot(aes(x=dim2,y=dim1))+geom_text(aes(label=gene))+
            geom_hline(yintercept = seq(-dim1,0,5)[-c(1,length(seq(-dim1,0,5)))],color="black",linetype=5)+
            scale_x_continuous("",expand = c(0,0),limits = c(0,10))+
            scale_y_continuous("",expand = c(0,0),limits = c(-dim1,0))+
            labs(title = paste(unique(plot_gene$program),collapse = "; "))+
            theme_bw()+
            theme(
                panel.grid = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(),
                plot.title = element_text(hjust = 0.5,size = 20)
            )
            
            library(patchwork)
            plot3=plot1+plot2+plot_layout(widths = c(1,2))
            ggsave(filename = paste(dir_output,"/",i,"_program_gene.heatmap.pdf",sep = ""),plot = plot3,width = 46,height = 16,units = "cm")
        }
    }

    setwd(fig_outdir)
    # step 3 use funtion step1 and step2
    library(reticulate)

    ### 调用子环境的python
    use_condaenv(condaenv = "cnmf_env", required = T,conda = "/mnt/share01/tools/bin/conda")
    py_config() #如果显示cnmf_env环境里面的python就OK

    ### 第一步 ################
    # 表达矩阵文件命名尽量和我的示例文件一致(HNSCC17.count.txt)，保证第一个点(.)前面是样本名称/ID。这些文件都保存在某个文件夹里面，比如示例中的"count_data"
    # source("1.R")

    step1(
        dir_input = "count_data",
        dir_output = "res1", 
        k = 3:10, 
        iteration = 100
    ) # 这里为了演示方便，取值都比较小

    # ###
    # 到这儿先停一停，step1会在dir_output下面为每一个样本生成一个文件夹。
    # 此外还有一个k_selection文件夹以及一个k_selection.txt文件。
    # k_selection文件夹里面有每一个样本选k值的图片（点击图片可以在浏览器打开），以此为依据，将k_selection.txt补充完整（可以直接在Rstudio打开）。
    # k_selection.txt有两列，以逗号分割，第二列是空的，需要你填上去。
    # ###


    ### 第二步 ################
    # source("2.R")

    step2(
        dir_input = "res1",
        dir_output = "res2", 
        dir_count = "count_data", 
        usage_filter = 0.02, 
        top_gene = 30, 
        cor_min = 0, 
        cor_max = 0.8
    )

    #过滤阈值，10X 0.01 smart-seq 0.03，这是我用过的阈值，不一定适用于所有情况
    #cor_min可以设为0，尽管真实范围会小于0；因为原文用的是0.6，所以我这里也限制最大相关系数为0.6，主要是为了图明显一些。试过10X数据，相关系数比较高(dropout导致的假象)，所以没有限制
    #top_gene按照20 30 50来设计的
    #有一个color参数，可以给样本涂色，格式是命名后的字符串向量；函数里面有一个默认的配色，也还OK。

    ### 实际分析中的一点点经验 ################
    #1. 相关性聚类这一步需要多做几次，所以原始画图的数据我都导出了，方便你自己改图
    #根据每个program的大致功能，如果发现另一种功能的program聚到某一种meta模块里面，这时可以将乱入的program删掉，再做一次相关性聚类。
    #如果我们认定应该属于同一个meta模块的program分散在两个地方，可以试试调整聚类方法(参数clustering_method)，或者像"_1""_2"这样定义
    #2. 对top基因做注释，尽量选择和研究背景有关联的基因集，比如研究肿瘤细胞异质性时，选择和肿瘤相关的基因集


    library(tidyverse)

    mat.files=dir("res2/",pattern = "dt_0_02.txt$")

    all.mat=data.frame()
    for (fi in mat.files) {
        cat("Parsing ",fi,"\n")
        tmp.mat=read.table(paste0("res2/",fi),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
        tmp.mat=as.data.frame(t(tmp.mat))
        sampleid=str_replace(fi,"\\..*$","")
        colnames(tmp.mat)=paste(sampleid,colnames(tmp.mat),sep = "_")
        tmp.mat$gene=rownames(tmp.mat)
        
        if (sampleid == "ZP02T") {
            all.mat=tmp.mat
        }else{
            all.mat = all.mat %>% full_join(tmp.mat,by="gene") #元素的并集进行合并
        }
    }


    # 对于某一个模块
    signature.programs = c(
        "ZR04T_1",
        "ZP09T_3",
        "ZR05T_1",
        "ZP06T_4",
        "ZR01T_4"
    )

    signature.programs = c(
        "ZP06T_2",
        "ZR01T_2",
        "ZP02T_1",
        "ZR03T_3",
        "ZP05T_3",
        "ZP10T_1",
        "ZR07T_2",
        "ZP07T_3",
        "ZP04T_2",
        "ZP08T_1",
        "ZP03T_2",
        "ZR01T_3"
    )
    
    signature.programs = c(
        "ZP02T_8",
        "ZP06T_5",
        "ZP04T_4",
        "ZR04T_3",
        "ZP05T_1",
        "ZR03T_1",
        "ZR05T_3",
        "ZP08T_3",
        "ZP09T_2",
        "ZR07T_1"
    )

    signature.programs = c(
        "ZP09T_1",
        "ZR04T_2",
        "ZP05T_2",
        "ZR05T_2",
        "ZP04T_1",
        "ZR07T_3",
        "ZP08T_4",
        "ZR01T_5",
        "ZP07T_1",
        "ZP03T_1"
    )


    signature.programs = c(
        "ZP02T_3",
        "ZR01T_1",
        "ZR07T_2",
        "ZP08T_4",
        "ZP04T_2",
        "ZP07T_2",
        "ZP10T_2",
        "ZR03T_3",
        "ZP04T_6",
        "ZP05T_8",
        "ZP06T_1",
        "ZP03T_6",
        "ZP05T_6"
    )
    signature.programs = c(
        "ZP06T_1",
        "ZP08T_2",
        "ZR07T_4",
        "ZP02T_7",
        "ZR03T_4",
        "ZP03T_3",
        "ZP04T_5",
        "ZR01T_1"
    )
    signature.programs = c(
        "ZP02T_3",
        "ZP03T_4"
    )

    signature.loading = all.mat[,c("gene",signature.programs)]

    used.gene=c()
    for (pi in signature.programs) {
    tmp.df=signature.loading[,c("gene",pi)]
    tmp.loading=tmp.df[,2]
    names(tmp.loading)=tmp.df[,1]
    
    tmp.loading=tmp.loading[!is.na(tmp.loading)]
    used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
    }
    used.gene=unique(used.gene)

    signature.loading=signature.loading[signature.loading$gene %in% used.gene,]
    rownames(signature.loading)=signature.loading$gene
    signature.loading$gene=NULL
    signature.loading[is.na(signature.loading)]<-0
    signature.loading$total_loading=rowSums(signature.loading)
    signature.loading$average_loading=signature.loading$total_loading / length(signature.programs)

    signature.loading=signature.loading%>%arrange(desc(average_loading))
    gene_set <- head(rownames(signature.loading),50)




    library(Seurat)
    # 
    # tmp.gene=c("TUBA1B", "HMGB2",  "TOP2A",  "UBE2C",  "NUSAP1", "MKI67",  "PBK",
    #            "BIRC5",  "CDK1",   "AURKB",  "UBE2T",  "CKS1B",  "H2AFZ",  "TPX2",   
    #            "TK1",    "CCNA2",  "GTSE1", "CEP55",  "KIF23",  "CENPF",  "CKS2",
    #            "HMGB1",  "PTTG1",  "CDCA3",  "CDKN3",  "PRC1",   "NUF2",   "CCNB1",
    #            "CCNB2",  "FOXM1")
    sum(gene_set %in% c(Seurat::cc.genes.updated.2019$s.genes,Seurat::cc.genes.updated.2019$g2m.genes))







    # GeneSet Analysis
    entrezid_frame <- bitr(gene_set,
        fromType = "SYMBOL", # 输入为SYMBOL格式
        toType = "ENTREZID", # 转为ENTERZID格式
        OrgDb = "org.Hs.eg.db"
    )

    data(geneList, package = "DOSE")

    ego_ALL <- enrichGO(
        gene = entrezid_frame$ENTREZID,
        universe = names(geneList), # 背景基因集
        OrgDb = org.Hs.eg.db, # 没有organism="human"，改为OrgDb=org.Hs.eg.db
        # keytype = 'ENSEMBL',
        ont = "ALL", # 也可以是 CC  BP  MF中的一种
        pAdjustMethod = "BH", # 矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
        pvalueCutoff = 1, # P值会过滤掉很多，可以全部输出
        qvalueCutoff = 1,
        readable = TRUE
    )

    pdf("GO_GSEA.pdf", width = 8, height = 9)
    dotplot(ego_ALL, title = "EnrichmentGO_MF_dot")
    dev.off()

    m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
        dplyr::select(gs_name, entrez_gene)
    em <- enricher(entrezid_frame$ENTREZID, TERM2GENE = m_t2g)
    
    pdf("HALLMARK.pdf", width = 8, height = 9)
    dotplot(em, title = "EnrichmentHallMark_dot")
    dev.off()










}






# Figure 14 SCEVAN demo
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure14")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)
    load("/P04_MalignantCell/P999_Malignant_Analysis/Figure14/MGH106_data.RData")

    devtools::reload(pkgload::inst("SCEVAN"))
    library(SCEVAN)
    suppressPackageStartupMessages(library(ggtree))
    results <- SCEVAN::pipelineCNA(count_mtx,
        sample = "MGH106", 
        par_cores = 30, 
        SUBCLONES = T, 
        plotTree = F
    )

}


# Figure 15 SCEVAN ZPT and ZRT
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure15")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)

    
    devtools::reload(pkgload::inst("SCEVAN"))
    suppressPackageStartupMessages(library(ggtree))

    for (tumor_sample_id in names(tumor_sample_colors)) {
        sample_outdir <- file.path(fig_outdir,tumor_sample_id)
        if (!file.exists(sample_outdir)) {
            dir.create(sample_outdir)
        }
        setwd(sample_outdir)
        message("Parsing ", tumor_sample_id)
        tumor_sce <- subset(global_sce, orig.ident %in% c(tumor_sample_id) & clusters %in% c("EPCAM+"))
        cat("Total cells:", length(Cells(tumor_sce)))
        if (length(Cells(tumor_sce)) < 100) {
            next
        }
        results <- SCEVAN::pipelineCNA(tumor_sce@assays$RNA@counts,
            sample = paste0(tumor_sample_id, "_"),
            par_cores = 2,
            SUBCLONES = TRUE,
            plotTree = TRUE
        )
        saveRDS(results, file = paste0(tumor_sample_id,"_scevan.rds"),compress = F)
    }



 




        # 保存之后生成热图
        saveRDS(gsea_result_frame,file = "gsea_result.rds",compress = F)

        # generate heatmap figures
        all_gsea_frame <- gsea_result_frame %>%
            select(c("pathway", "sample_id", "subclone_id", "NES")) %>%
            mutate(sample_clone_id = paste0(sample_id, "_", subclone_id)) %>%
            select(-c("sample_id", "subclone_id")) %>%
            spread(sample_clone_id, NES) %>%
            mutate_if(is.numeric, ~replace(., is.na(.), 0))
        
        all_gsea_frame <- as.data.frame(all_gsea_frame)
        row.names(all_gsea_frame) <- all_gsea_frame$pathway
        all_gsea_frame$pathway <- NULL
        all_gsea_frame <- all_gsea_frame[order(rowSums(all_gsea_frame != 0), decreasing = T),]
        ht <- Heatmap(
            as.matrix(head(all_gsea_frame,n = 20)),
            col = colorRamp2(c(0, 2), c("#2B74C8","#FAFA0E")),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 36),
            column_names_gp = gpar(fontsize = 36),
            rect_gp = gpar(col = "white", lwd = 3),
            heatmap_legend_param = list(title = "",
                    grid_width = unit(1, "cm"),
                    legend_height = unit(8, "cm"),
                    title_gp = gpar(fontsize = 20),
                    labels_gp = gpar(fontsize = 30)),
            cluster_rows = FALSE,
            cluster_columns = FALSE
        )
        png("gsea_heatmap.png",width = 2080 * 2,height = 1080 * 2)
        draw(ht, padding = unit(c(260, 360, 20, 20), "mm"))
        dev.off()

        pdf("gsea_heatmap.pdf",width = 28 * 2,height = 20 * 2)
        draw(ht, padding = unit(c(260, 360, 20, 20), "mm"))
        dev.off()





    
}




# Figure 18 SCEVAN ZPT and ZRT
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure18")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)
    library(dior)
    
    for (tumor_sample_id in names(tumor_sample_colors)) {
        tumor_sample_id <- "PH02T"
        message(sprintf("Parsing %s", (tumor_sample_id)))
        sample_outdir <- file.path(fig_outdir,tumor_sample_id)
        if (!file.exists(sample_outdir)) {
            dir.create(sample_outdir)
        }
        setwd(sample_outdir)
        # tumor_sce <- subset(global_sce, sample_id %in% c(tumor_sample_id) & clusters %in% c("EPCAM+"))
        tumor_sce <- subset(global_sce, clusters %in% c("EPCAM+"))
        cat("Total cells:", length(Cells(tumor_sce)))
        # if (length(Cells(tumor_sce)) < 100) {
        #     next
        # }
        tumor_sce <- RunPCA(object = tumor_sce, resolution = 1)
        tumor_sce <- FindNeighbors(object = tumor_sce)
        tumor_sce <- FindClusters(object = tumor_sce)

        tumor_sce <- RunUMAP(tumor_sce, dims = 1:10)
        gp <- DimPlot(tumor_sce, reduction = "umap", 
                    pt.size = 3, 
                    label = TRUE, 
                    label.size = 4)
        ggsave(gp, file = "PH02T_UMAP.png",w = 8,h = 8)


        write_h5(tumor_sce, 
                file = sprintf("%s_sce.h5",tumor_sample_id), 
                object.type = 'seurat', 
                assay.name = 'RNA', 
                save.graphs = TRUE, 
                save.scale=FALSE)

        saveRDS(results, file = paste0(tumor_sample_id,"_scevan.rds"),compress = F)
    }



 




        # 保存之后生成热图
        saveRDS(gsea_result_frame,file = "gsea_result.rds",compress = F)

        # generate heatmap figures
        all_gsea_frame <- gsea_result_frame %>%
            select(c("pathway", "sample_id", "subclone_id", "NES")) %>%
            mutate(sample_clone_id = paste0(sample_id, "_", subclone_id)) %>%
            select(-c("sample_id", "subclone_id")) %>%
            spread(sample_clone_id, NES) %>%
            mutate_if(is.numeric, ~replace(., is.na(.), 0))
        
        all_gsea_frame <- as.data.frame(all_gsea_frame)
        row.names(all_gsea_frame) <- all_gsea_frame$pathway
        all_gsea_frame$pathway <- NULL
        all_gsea_frame <- all_gsea_frame[order(rowSums(all_gsea_frame != 0), decreasing = T),]
        ht <- Heatmap(
            as.matrix(head(all_gsea_frame,n = 20)),
            col = colorRamp2(c(0, 2), c("#2B74C8","#FAFA0E")),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 36),
            column_names_gp = gpar(fontsize = 36),
            rect_gp = gpar(col = "white", lwd = 3),
            heatmap_legend_param = list(title = "",
                    grid_width = unit(1, "cm"),
                    legend_height = unit(8, "cm"),
                    title_gp = gpar(fontsize = 20),
                    labels_gp = gpar(fontsize = 30)),
            cluster_rows = FALSE,
            cluster_columns = FALSE
        )
        png("gsea_heatmap.png",width = 2080 * 2,height = 1080 * 2)
        draw(ht, padding = unit(c(260, 360, 20, 20), "mm"))
        dev.off()

        pdf("gsea_heatmap.pdf",width = 28 * 2,height = 20 * 2)
        draw(ht, padding = unit(c(260, 360, 20, 20), "mm"))
        dev.off()

    # 画一个热图展示GSEA的通路，选细胞多的样本

}


# plot significant GSEA pathway genes using Seurat
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure16")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)
    tumor_sample_id <- "ZP03T"
    tumor_sce <- subset(global_sce, orig.ident %in% c("ZP03T") &
        clusters %in% c("EPCAM+"))
    scevan_metadata <- readRDS(file = file.path(outdir,"Figure15",tumor_sample_id, paste0(tumor_sample_id, "_scevan.rds")))
    tumor_sce <- AddMetaData(tumor_sce, metadata = scevan_metadata)


    stat_table <- table(tumor_sce$subclone)
    tumor_sce <- subset(malignant_sce, subclone %in% names(stat_table))
    Idents(tumor_sce) <- "subclone"


    mdb_h <- msigdbr(species = "Homo sapiens", category = "H")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    p <- DoHeatmap(tumor_sce,features = fgsea_sets$HALLMARK_HYPOXIA,slot = "scale.data")
    ggsave("pathway_heatmap.png",width = 8,height = 20)

}




# Figure 18 SCEVAN ZPT and ZRT
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure19")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)
    library(dior)
    
    scevan_dir <- "/P04_MalignantCell/P999_Malignant_Analysis/Figure15"
    scevan_frame <- data.frame()
    for (tumor_sample_id in names(tumor_sample_colors)) {
        message(sprintf("Parsing %s", (tumor_sample_id)))
        scevan_file <- file.path(scevan_dir,tumor_sample_id,paste0(tumor_sample_id,"_scevan.rds"))
        if (!file.exists(scevan_file)) {
            next
        }
        scevan <- readRDS(scevan_file)
        scevan_frame <- rbind(scevan_frame,scevan)
    }

    epithelial_sce <- subset(global_sce,cells = row.names(scevan_frame))
    epithelial_sce <- AddMetaData(epithelial_sce, metadata = scevan_frame)



    write_h5(epithelial_sce, 
            file = "epithelial_sce.h5", 
            object.type = 'seurat', 
            assay.name = 'RNA', 
            save.graphs = TRUE, 
            save.scale=FALSE)

    malignant_sce <- subset(epithelial_sce,class %in% c("tumor"))
    malignant_sce$subclone <- paste0(malignant_sce$sample_id,"_",malignant_sce$subclone)

    library(msigdbr)
    library(GSVA)
    library(limma)

    genesets <- msigdbr(species ="Homo sapiens",category ="H")
    genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
    genesets <- split(genesets$gene_symbol,genesets$gs_name)

    Idents(malignant_sce) <-"subclone"
    expr <- AverageExpression(malignant_sce,assays ="RNA",slot = "data")
    expr <- expr$RNA
    expr <- expr[rowSums(expr)>0,]#选取非零基因
    expr <- as.matrix(expr)

    gsva.res <- gsva(expr,genesets,method = "zscore")

    gsva.df <- data.frame(Genesets=rownames(gsva.res),gsva.res,check.names = F)
    # write.csv(gsva.df,"gsva_res.csv",row.names = F)
    pdf(file = "subclone.pdf",width = 18,height = 8)
    pheatmap::pheatmap(gsva.res,
        show_colnames = T,
        scale = "row",
        angle_col = "45",
        color = colorRampPalette(c("navy","white","firebrick3"))(50))
    dev.off()



    # metabolism analysis
    gmtfile = "/P04_MalignantCell/P999_Malignant_Analysis/data/KEGG_metabolism.gmt"
    genesets <- clusterProfiler::read.gmt(gmtfile)
    colnames(genesets) <- c("gs_name","gene_symbol")
    genesets <- split(genesets$gene_symbol,genesets$gs_name)


    Idents(malignant_sce) <-"subclone"
    expr <- AverageExpression(malignant_sce,assays ="RNA",slot = "data")
    expr <- expr$RNA
    expr <- expr[rowSums(expr)>0,]#选取非零基因
    expr <- as.matrix(expr)

    gsva.res <- gsva(expr,genesets,method = "zscore")

    gsva.df <- data.frame(Genesets=rownames(gsva.res),gsva.res,check.names = F)
    # write.csv(gsva.df,"gsva_res.csv",row.names = F)
    pdf(file = "metabolism_subclone.pdf",width = 18,height = 12)
    pheatmap::pheatmap(gsva.res,
        show_colnames = T,
        scale = "row",
        angle_col = "45",
        color = colorRampPalette(c("#2F398D","white","#9C1F28"))(50))
    dev.off()




}






























