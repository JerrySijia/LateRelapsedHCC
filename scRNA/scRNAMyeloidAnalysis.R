# load library
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)
library(ggsignif)
library(circlize)
library(future)
library(ggrepel)
library(dplyr)
library(destiny)
library(msigdbr)
library(fgsea)

#myeloid_seurat_rds <- "/P03_MyeloidCell/P03_Myeloid_UMAPByCellType/myeloid.rds"
myeloid_seurat_rds <- "/P03_MyeloidCell/P999_Myeloid_Analysis/Figure1/myeloid_sce.rds"
cell_type <- "/P03_MyeloidCell/P999_Myeloid_Analysis/data/cell.type.txt"
sample_list <- "/P03_MyeloidCell/P999_Myeloid_Analysis/data/sample_list.txt"

celltype_group <- "/P01_GlobalCell/P999_GlobalCell_Analysis/data/celltype_group.txt"
inflammatory_signature_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/inflammatory_signature.txt"
dc_signature_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/dc_signature.txt"
dc_marker_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/dc_marker.txt"
outdir <- "/P03_MyeloidCell/P999_Myeloid_Analysis"


# Read myeloid Seurat RDS Data
cat("Reading Myeloid Seurat RDS data...\n")
myeloid_sce <- readRDS(myeloid_seurat_rds)
cat("Complete.\n")

plan("multicore", workers = 8)
options(future.globals.maxSize = 150 * 1024^3)


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
    ZP04N = "#DCC1DD",
    ZP05T = "#CCE0F5",
    ZP05N = "#CCC9E6",
    ZP06T = "#625D9E",
    ZP06N = "#68A180",
    ZP07T = "#3A6963",
    ZP07N = "#968175"
)

group_cols <- c(
    PHN = "#ef9020",
    PHT = "#00af3e", 
    LRN = "#0081b4", 
    LRT = "#CD1E24"
)

celltype_colors <- c(
    "Monocyte-C1" = "#262C68",
    "Monocyte-C2" = "#CD1E24",
    "Monocyte-C3" = "#1E843F",
    "Monocyte-C4" = "#84278B",
    "Monocyte-C5" = "#0B6E78",
    "Macrophage-C1" = "#EF7A2A",
    "Macrophage-C2" = "#B968A5",
    "Macrophage-C3" = "#BDD240",
    "Macrophage-C4" = "#E63863",
    "Macrophage-C5" = "#E4C755",
    "DC-C1" = "#D6E7A3",
    "DC-C2" = "#DCC1DD",
    "DC-C3" = "#58A4C3",
    "Cycling" = "#AB3282"
)

#VIP Figure 1 UMAP Plot by Celltype
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure1")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)




    meta_frame <- myeloid_sce@meta.data %>% 
        dplyr::select(c("orig.ident"))
    meta_frame$sample_id <- gsub("ZP","PH",meta_frame$orig.ident)
    meta_frame$sample_id <- gsub("ZR","LR",meta_frame$sample_id)
    meta_frame$group_id <- paste0(substr(meta_frame$sample_id,1,2),substr(meta_frame$sample_id,5,5))
    myeloid_sce = AddMetaData(myeloid_sce, meta_frame)

    saveRDS(myeloid_sce,file = "myeloid_sce.rds",compress = F)

    data = myeloid_sce@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(clusters = myeloid_sce@meta.data$clusters)

    # step2 获取要添加标签的位置
    class_avg <- data %>% 
        group_by(clusters) %>% 
        summarise(
        UMAP_1 = median(UMAP_1),
        UMAP_2 = median(UMAP_2)
    )

    # step2 绘图
    umap <- ggplot(data ,aes(x=UMAP_1,y=UMAP_2))+
        geom_point(aes(color=clusters),size = 0.5, alpha = 0.2)+ 
        scale_color_manual(values = celltype_colors)+
        geom_text(aes(label = clusters), data = class_avg, size = 3)+
        xlab('UMAP dimension 1')+
        ylab('UMAP dimension 2')+
        theme_bw()+
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                legend.title = element_blank(),
                axis.title.x = element_text(size=18), 
                axis.title.y = element_text(size=18),
                axis.text = element_text(size=18),
                plot.margin=unit(rep(1,4),'cm'),
                legend.text = element_text(size=18),
                legend.key.size = unit(0.4, "inches")
        ) + guides(colour = guide_legend(override.aes = list(size = 6)))

    ggsave("MyeloidUMAPByCellType.pdf",plot = umap,width = 12,height = 8)
    ggsave("MyeloidUMAPByCellType.png",plot = umap,width = 12,height = 8)

}


# Figure 2 Marker VlnPlot by Cell type
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure2")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    myeloid_marker_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/TAMsM1M2_markers.txt"
    myeloid_marker_frame = read.table(file = myeloid_marker_file,sep = "\t",header = F,stringsAsFactors = F)

    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
        p <- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = "clusters", log = TRUE, ...) +
            xlab("") +
            ylab(feature) +
            ggtitle("") +
            scale_x_discrete(limits = names(celltype_colors)) +
            theme(
                legend.position = "none",
                # axis.text.x = element_text(size = rel(1), angle = 30),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(),
                axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
                plot.margin = plot.margin
            )
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, ...))
        plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
            theme(
                axis.text.x = element_text(size = rel(0.6), angle = 60),
                axis.ticks.x = element_line()
            )
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
        return(p)
    }

    p <- StackedVlnPlot(myeloid_sce, myeloid_marker_frame$V2,
        slot = "scale.data", 
        pt.size = 0.6, 
        cols = celltype_colors
    )

    ggsave(filename = "Figure2-TAMM1M2.png", plot = p, width = 8, height = 16, units = c("in"))
    ggsave(filename = "Figure2-TAMM1M2.pdf", plot = p, width = 8, height = 16, units = c("in"))

}

# Figure 3: Plot myeloid Marker Gene Heatmap By Celltype
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure3")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    cell_marker <- "/P03_MyeloidCell/P999_Myeloid_Analysis/data/selected_markers.txt"
    Idents(myeloid_sce) <- "seurat_clusters"
    cell_marker_frame = read.table(file = cell_marker,
                sep = "\t", 
                header = T,
                stringsAsFactors = F)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$gene),]

    AverageExp <- AverageExpression(myeloid_sce,
        features = cell_marker_frame$gene, 
        group.by = "clusters"
    )

    expr <- AverageExp$RNA

    heatmapColorRamp = colorRampPalette(c("#0658D6","#FDF8EE","#ED7A03"))(64)

    expr <- ScaleData(expr)
    color_array <- c("#949483",
                    "#F47B7B",
                    "#9F1F5C",
                    "#EF9020",
                    "#00AF3E",
                    "#85B7E2",
                    "#29245C",
                    "#FFD616",
                    "#E5352B",
                    "#FFD616",
                    "#E990AB",
                    "#0081B4",
                    "#96CBB3",
                    "#91BE3E",
                    "#39A6DD",
                    "#EB0973",
                    "#DDE2E0",
                    "#333C41")

    ordered_celltype_array <- sort(colnames(expr))

    celltype_color_structure <- structure(color_array[1:(length(unique(cell_marker_frame$celltype)))],
                                        names = unique(cell_marker_frame$celltype)) 

    ha = rowAnnotation(labels = cell_marker_frame$celltype,
                    col = list(bar = celltype_color_structure))
    expr <- expr[, ordered_celltype_array]
    ht <- Heatmap(as.matrix(expr),
            col = heatmapColorRamp,
            right_annotation = ha,
            rect_gp = gpar(col = "white", lwd = 1),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 16),
            row_split = cell_marker_frame$celltype,
            #row_title_side = "right",
            #row_title_rot = 0,
            row_title = NULL,
            cluster_columns = F,
            cluster_rows = F)

    png("MyeloidMarkerHeatmap.png",width = 1080,height = 1080 * 1.5)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

    pdf("MyeloidMarkerHeatmap.pdf",width = 8,height = 16)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

}


#VIP Figure5C Macrophage cell Monocyte, M1, M2 and MDSC signature VlnPlot and DotPlot for each Macrophage type
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure4")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    M1M2_signature_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/M1M2_signature.txt"

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macrophage", myeloid_sce@meta.data$clusters)]))

    M1M2_signature_frame <- read.table(file = M1M2_signature_file, sep = "\t",header = F,stringsAsFactors = F)
    M1M2_signature_frame <- M1M2_signature_frame %>% filter(!(V1 %in% c("Monocyte")))
    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)

    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(M1M2_signature_frame, M1M2_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        macrophage_sce <- AddModuleScore(
                        object = macrophage_sce,
                        features = list(gene_array),
                        name = signature)
        colnames(macrophage_sce@meta.data) <- gsub(colnames(macrophage_sce@meta.data),
            pattern = paste0(signature, 1),
            replacement = signature
        )

    }


    p <- DotPlot(
        object = macrophage_sce,
        features = names(split_signature_list),
        dot.scale = 12,
        cols = c("#0E01FA", "#FB010F"),
        group.by = "clusters",
    ) + theme(
        axis.text.x = element_text(size = 30, angle = 90),
        axis.text.y = element_text(size = 30)
    ) + theme_bw()

    ggsave("Figure4-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure4-DotPlot.pdf",p,width = 8,height = 8,units = "in")

}


# Figure 5 Macrophage cell inflammatory signature VlnPlot and DotPlot for each Macrophage type
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure5")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macrophage", myeloid_sce@meta.data$clusters)]))

    inflammatory_signature_frame <- read.table(file = inflammatory_signature_file, sep = "\t",header = F,stringsAsFactors = F)
    
    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)

    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(inflammatory_signature_frame, inflammatory_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        macrophage_sce <- AddModuleScore(
                        object = macrophage_sce,
                        features = list(gene_array),
                        name = signature)
    }
    
    plot <- VlnPlot(macrophage_sce,
        features = paste0(names(split_signature_list),"1"),
        group.by = "clusters",
        col = celltype_colors,
        ncol = 2
    )
    
    ggsave("Figure5-VlnPlot.png",plot = plot,width = 8,height = 4.5,units = "in")
    ggsave("Figure5-VlnPlot.pdf",plot = plot,width = 8,height = 4.5,units = "in")


    p <- DotPlot(
        object = macrophage_sce,
        features = paste0(names(split_signature_list), "1"),
        dot.scale = 12,
        cols = c("#0E01FA", "#FB010F"),
        group.by = "clusters",
    ) + theme(axis.text.x = element_text(angle = 90)) + theme_bw()

    ggsave("Figure5-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure5-DotPlot.pdf",p,width = 8,height = 8,units = "in")

}


#VIP Figure 6 DC cell signature VlnPlot and DotPlot for each DC type
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure6")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    dc_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("DC", myeloid_sce@meta.data$clusters)]))

    dc_signature_frame <- read.table(file = dc_signature_file, sep = "\t",header = F,stringsAsFactors = F)
    
    dc_sce <- subset(myeloid_sce, clusters %in% dc_celltype_array)

    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(dc_signature_frame, dc_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        dc_sce <- AddModuleScore(
                        object = dc_sce,
                        features = list(gene_array),
                        name = signature)
    }
    
    plot <- VlnPlot(dc_sce,
        features = paste0(names(split_signature_list),"1"),
        group.by = "clusters",
        col = celltype_colors,
        ncol = 2
    )
    
    ggsave("Figure6-VlnPlot.png",plot = plot,width = 8,height = 4.5,units = "in")
    ggsave("Figure6-VlnPlot.pdf",plot = plot,width = 8,height = 4.5,units = "in")


    p <- DotPlot(
        object = dc_sce,
        features = paste0(names(split_signature_list), "1"),
        dot.scale = 12,
        cols = c("#0E01FA", "#FB010F"),
        group.by = "clusters",
    ) + theme(axis.text.x = element_text(angle = 90)) + theme_bw()

    ggsave("Figure6-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure6-DotPlot.pdf",p,width = 8,height = 8,units = "in")

}




# Figure 7 TAMs, M1 and M2 Marker VlnPlot by Macrophage Cell type
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure7")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macrophage", myeloid_sce@meta.data$clusters)]))

    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)

    TAMsM1M2_marker_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/TAMsM1M2_markers.txt"

    tams_marker_frame = read.table(file = TAMsM1M2_marker_file,sep = "\t",header = F,stringsAsFactors = F)

    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
        p <- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = "clusters", log = TRUE, ...) +
            xlab("") + ylab(feature) + ggtitle("") +
            theme(
                legend.position = "none",
                # axis.text.x = element_text(size = rel(1), angle = 30),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(),
                axis.title.y = element_text(size = rel(1.5), angle = 0, vjust = 0.5),
                plot.margin = plot.margin
            )
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, ...))
        plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
            theme(
                axis.text.x = element_text(size = rel(1.0), angle = 30),
                axis.ticks.x = element_line()
            )
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
        return(p)
    }

    p <- StackedVlnPlot(macrophage_sce, tams_marker_frame$V2,
        slot = "scale.data", 
        pt.size = 0.3, 
        cols = celltype_colors
    )

    ggsave(filename = "Figure7.png", plot = p, width = 8, height = 25, units = c("in"))
    ggsave(filename = "Figure7.pdf", plot = p, width = 8, height = 25, units = c("in"))


}



# Figure 8 Marker of Macrophage cell type VlnPlot by group (ZPT, ZRT,ZPN, and ZRN)
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure8")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    samples_list <- read.table(sample_list, header = F, sep = "\t")
    group_array <- samples_list[, "V2"]
    names(group_array) <- samples_list[, "V1"]
    myeloid_sce$group <- group_array[myeloid_sce$orig.ident]

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macrophage", myeloid_sce@meta.data$clusters)]))

    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)
    markers <- c(
        "PDCD1LG2",
        "CD80",
        "CD86",
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "PVR",
        "LGALS9",
        "TNFRSF14",
        "SPP1"
    )

    #orig.ident
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
        compaired <- list(
            c("ZPT", "ZRT"),
            c("ZPN", "ZRN"),
            c("ZPT", "ZPN"),
            c("ZRT", "ZRN")
        )

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'group', log=TRUE,... ) +
            xlab("") + 
            ylab(feature) + 
            ggtitle("") +
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
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 3)
        return(p)
    }

    p <- StackedVlnPlot(macrophage_sce,
        markers,
        pt.size = 0, 
        cols = group_cols
    )


    ggsave(filename = "Figure8.png", plot = p, width = 12, height = 12, units = c("in"))
    ggsave(filename = "Figure8.pdf", plot = p, width = 12, height = 12, units = c("in"))

}



# Figure 9 DC marker DotPlot for DC celltype
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure9")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    Idents(myeloid_sce) <- "seurat_clusters"


    dc_celltype_array <- unique(as.character(dc_sce@meta.data$clusters[grepl("DC", myeloid_sce@meta.data$clusters)]))

    dc_sce <- subset(myeloid_sce, clusters %in% dc_celltype_array)

    cell_marker_frame = read.table(file = dc_marker_file,
                sep = "\t", 
                header = F,
                stringsAsFactors = F)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$V2),]

    dp <- DotPlot(dc_sce,
        features = cell_marker_frame$V2,
        group.by = "clusters"
    ) +
        coord_flip() +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        ) + scale_color_gradientn(values = seq(0, 1, 0.2), colours = c("#330066", "#336699", "#66CC66", "#FFCC33")) +
        labs(x = NULL) +
        guides(size = guide_legend(order = 3))

    ggsave(
        filename = "Figure9.png",
        plot = dp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "Figure9.pdf", 
       plot = dp, 
       width = 8, 
       height = 12, 
       units = c("in"))

}




#VIP Figure 10 DC marker Heatmap for DC celltype
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure10")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    Idents(myeloid_sce) <- "seurat_clusters"


    dc_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("DC", myeloid_sce@meta.data$clusters)]))

    dc_sce <- subset(myeloid_sce, clusters %in% dc_celltype_array)

    cell_marker_frame = read.table(file = dc_marker_file,
                sep = "\t", 
                header = F,
                stringsAsFactors = F)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$V2),]

    AverageExp <- AverageExpression(dc_sce,slot = 'data',group.by = "clusters",features=cell_marker_frame$V2)

    expr <- AverageExp$RNA
    colorRamp = colorRampPalette(c("#0658D6","#FDF8EE","#ED7A03"))(64)
    scale_row <- function (x){
        m = apply(x, 1, mean, na.rm = T)
        s = apply(x, 1, sd, na.rm = T)
        return((x - m)/s)
    }

    expr <- scale_row(as.matrix(expr))

    # Complex Heatmap
    ht <- Heatmap(expr,
            col = colorRamp,
            rect_gp = gpar(col = "white", lwd = 1), 
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 15),
            column_names_gp = gpar(fontsize = 20),
            row_title_side = "right",
            row_title_rot = 0,
            cluster_columns = F,
            cluster_rows = F)

    png("Figure10.png",width = 768,height = 1080)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()


    pdf("Figure10.pdf",width = 3,height = 6) 
    draw(ht,padding = unit(c(2, 2, 2, 2), "mm"))
    dev.off()

}



# Figure 11 Myeloid Cell type Bar Plot
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure11")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    immune_cell_stat_frame <- as.data.frame.array(table(
        as.character(myeloid_sce$orig.ident),
        as.character(myeloid_sce$clusters)
    ))

    write.table(
        x = immune_cell_stat_frame,
        file = "myeloid_immune_stat.xls",
        sep = "\t",
        col.name = TRUE,
        row.name = TRUE,
        quote=FALSE
    )

    immune_cell_stat_frame$sample_group <- paste0(
        substr(row.names(immune_cell_stat_frame), 1, 2),
        substr(row.names(immune_cell_stat_frame), 5, 6)
    )
    immune_cell_stat_frame$sample <- row.names(immune_cell_stat_frame)
    row.names(immune_cell_stat_frame) <- NULL

    bp <- immune_cell_stat_frame %>% gather(celltype,count,-sample_group)  %>% 
        ggplot(aes(
                x = sample_group,
                y = count,
                fill = factor(celltype)
            )
            ) +
            geom_bar(
                position = "fill",
                stat = "identity", width = 0.5
            ) +
            scale_fill_manual(values = celltype_colors) +
            xlab("") +
            ylab("Percent") +
            theme_bw() +
            guides(fill = guide_legend(title = NULL)) +
            theme(
                axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 25),
                legend.text = element_text(size = 25),
                axis.text.y = element_text(size = 20),
                axis.title.y = element_text(size = 30)
            )

    pdf("Figure11.pdf", w = 6, h = 8)
    print(bp)
    dev.off()
    png("Figure11.png", width = 1200, height = 1600)
    print(bp)
    dev.off()
}




#VIP Figure 12: Macrophage cell type frequency BoxPlot by PHT, LRT, PHN and LRN
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure12")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)



    cell_stat_frame <- as.data.frame.array(table(
        as.character(myeloid_sce$sample_id),
        as.character(myeloid_sce$clusters)
    ))

    celltype_array = colnames(cell_stat_frame)

    cell_stat_frame$Total_cell <- rowSums(cell_stat_frame)

    cell_stat_frame$Sample <- row.names(cell_stat_frame)
    row.names(cell_stat_frame) <- NULL

    cell_stat_frame$sample_group <- paste0(substr(cell_stat_frame$Sample,1,2),
                                        substr(cell_stat_frame$Sample,5,6))

    for(celltype in celltype_array){
        cell_stat_frame[celltype] = cell_stat_frame[celltype] / cell_stat_frame["Total_cell"]
    }
    macrophage_celltype_array = colnames(cell_stat_frame %>% dplyr::select(starts_with("Macrophage")))
    cell_stat_frame <- cell_stat_frame[, c(macrophage_celltype_array, "Sample", "sample_group")]

    cell_stat_frame["Macrophage"] = rowSums(cell_stat_frame[macrophage_celltype_array])

    cell_stat_frame[c("Total_cell")] <- NULL

    melt_cell_stat_frame = melt(cell_stat_frame,id = c("Sample","sample_group"))
    compaired <- list(c("PHT","PHN"),c("LRT","LRN"),c("PHN","LRN"),c("PHT","LRT"))

    p1 <- ggplot(melt_cell_stat_frame,aes(x = sample_group,y = value,fill = sample_group)) + 
        geom_boxplot(outlier.colour = NA) +
        scale_fill_manual(values = group_cols) +
        theme_bw() +
        theme(strip.background = element_rect(fill="white", colour="black", size=0),
        # remove legend title
                legend.title=element_blank()) +
        xlab("") + 
        ylab("") +
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    test.args = list(var.equal = TRUE,alternative = "greater"), # one side or two side?
                    test = wilcox.test) +
        facet_grid(~variable)

    ggsave("Figure12.png",p1,width = 10,height = 4.5,units = "in")
    ggsave("Figure12.pdf",p1,width = 10,height = 4.5,units = "in")

}



#VIP Figure 13: DC cell type frequency BoxPlot by PHT, LRT, PHN and LRN
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure13")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)



    cell_stat_frame <- as.data.frame.array(table(
        as.character(myeloid_sce$orig.ident),
        as.character(myeloid_sce$clusters)
    ))

    celltype_array = colnames(cell_stat_frame)

    cell_stat_frame$Monocyte <- rowSums(cell_stat_frame)

    cell_stat_frame$Sample <- row.names(cell_stat_frame)
    row.names(cell_stat_frame) <- NULL

    cell_stat_frame$sample_group <- paste0(substr(cell_stat_frame$Sample,1,2),
                                        substr(cell_stat_frame$Sample,5,6))

    for(celltype in celltype_array){
        cell_stat_frame[celltype] = cell_stat_frame[celltype] / cell_stat_frame["Monocyte"]
    }
    macrophage_celltype_array = colnames(cell_stat_frame %>% dplyr::select(starts_with("DC")))
    cell_stat_frame <- cell_stat_frame[, c(macrophage_celltype_array, "Sample", "sample_group")]

    melt_cell_stat_frame = melt(cell_stat_frame,id = c("Sample","sample_group"))
    compaired <- list(c("ZPT","ZPN"),c("ZRT","ZRN"),c("ZPN","ZRN"),c("ZPT","ZRT"))

    p1 <- ggplot(melt_cell_stat_frame,aes(x = sample_group,y = value,fill = sample_group)) + 
    geom_boxplot(outlier.colour = NA) +
    scale_fill_manual(values = group_cols) +
    theme_bw() +
    theme(strip.background = element_rect(fill="white", colour="black", size=0),
    # remove legend title
            legend.title=element_blank()) +
    xlab("") + 
    ylab("") +
    geom_signif(comparisons = compaired,
                step_increase = 0.1,
                map_signif_level = F,
                #test.args = list(var.equal = TRUE,alternative = "greater"), # one side or two side?
                test = wilcox.test) +
    facet_grid(~variable)

    ggsave("Figure13.png",p1,width = 8,height = 4.5,units = "in")
    ggsave("Figure13.pdf",p1,width = 8,height = 4.5,units = "in")

}


#VIP Figure 14: Monocyte cell type frequency BoxPlot by PHT, LRT, PHN and LRN
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure14")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cell_stat_frame <- as.data.frame.array(table(
        as.character(myeloid_sce$sample_id),
        as.character(myeloid_sce$clusters)
    ))

    celltype_array = colnames(cell_stat_frame)

    cell_stat_frame$Total_cell <- rowSums(cell_stat_frame)

    cell_stat_frame$Sample <- row.names(cell_stat_frame)
    row.names(cell_stat_frame) <- NULL

    cell_stat_frame$sample_group <- paste0(substr(cell_stat_frame$Sample,1,2),
                                        substr(cell_stat_frame$Sample,5,6))

    for(celltype in celltype_array){
        cell_stat_frame[celltype] = cell_stat_frame[celltype] / cell_stat_frame["Total_cell"]
    }
    monocyte_celltype_array = colnames(cell_stat_frame %>% dplyr::select(starts_with("Monocyte")))
    cell_stat_frame <- cell_stat_frame[, c(monocyte_celltype_array, "Sample", "sample_group")]

    cell_stat_frame["Monocyte"] = rowSums(cell_stat_frame[monocyte_celltype_array])

    melt_cell_stat_frame = melt(cell_stat_frame,id = c("Sample","sample_group"))
    compaired <- list(c("PHT","PHN"),c("LRT","LRN"),c("PHN","LRN"),c("PHT","LRT"))

    p1 <- ggplot(melt_cell_stat_frame,aes(x = sample_group,y = value,fill = sample_group)) + 
        geom_boxplot(outlier.colour = NA) +
        scale_fill_manual(values = group_cols) +
        theme_bw() +
        theme(strip.background = element_rect(fill="white", colour="black", size=0),
        # remove legend title
                legend.title=element_blank()) +
        xlab("") + 
        ylab("") +
        geom_signif(comparisons = compaired,
                    step_increase = 0.1,
                    map_signif_level = T,
                    test.args = list(var.equal = TRUE,alternative = "less"), # one side or two side?
                    test = wilcox.test) +
        facet_grid(~variable)

    ggsave("Figure14.png",p1,width = 8,height = 4.5,units = "in")
    ggsave("Figure14.pdf",p1,width = 8,height = 4.5,units = "in")

}


# Figure 15: Plot Myeloid UMAP by marker
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure15")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    markers <- c("AIF1","APOE","FCN1","CLEC9A","LAMP3","CD207","LILRA4","CCL17","PCNA")
        
    p1 <- FeaturePlot(myeloid_sce,
        features = markers, 
        cols = c("lightgrey", "orange"), 
        pt.size = 0.5,
        raster = F,
        ncol = 3
    )

    ggsave(plot = p1, filename = paste0("Figure15.pdf"), w = 12, h = 8, units = c("in"))
    ggsave(plot = p1, filename = paste0("Figure15.png"), w = 12, h = 8, units = c("in"))
}



# Figure 16 Macrophage Vocalno Plot by ZPN and ZRN
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure16")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macrophage", myeloid_sce@meta.data$clusters)]))

    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)
    macrophage_normal_sce <- subset(macrophage_sce, group %in% c("ZPN", "ZRN"))

    if (T) {
        allCells = colnames(macrophage_normal_sce)
        allType = names(table(macrophage_normal_sce@meta.data$orig.ident))

        choose_Cells = unlist(lapply(allType, function(x) {
            cgCells = allCells[macrophage_normal_sce@meta.data$orig.ident == x]
            num = ceiling(dim(as.data.frame(cgCells))[1] * 0.3)
            cg = sample(cgCells, num)
            cg
        }))

        macrophage_normal_sce = macrophage_normal_sce[, allCells %in% choose_Cells]
    }
    
    diff_marker_frame <- FindMarkers(macrophage_normal_sce,
        ident.1 = "ZPN", 
        group.by = "group", 
        assay = "RNA", 
        slot = "counts", 
        logfc.threshold = 0, 
        min.pct = 0.01,
        test.use="DESeq2"
    )

    write.table(
        x = diff_marker_frame,
        file = "macrophage_normal_diffexp.gene.txt",
        sep = "\t",
        col.name = TRUE,
        row.name = F,
        quote=FALSE
    )

    log2FC_threshold = 0.8 
    p_value_threshold = 0.01
    diff_marker_frame[which(diff_marker_frame$p_val_adj < 0.01 & diff_marker_frame$avg_log2FC <= -log2FC_threshold),'sig'] <- 'ZRN'
    diff_marker_frame[which(diff_marker_frame$p_val_adj < 0.01 & diff_marker_frame$avg_log2FC >= log2FC_threshold),'sig'] <- 'ZPN'
    diff_marker_frame[which(diff_marker_frame$p_val_adj >= 0.01 | abs(diff_marker_frame$avg_log2FC) < log2FC_threshold),'sig'] <- 'None'
    # 横轴 log2FC，纵轴 -log10(adj.P.Val)，颜色表示差异

    diff_marker_frame$Gene <- row.names(diff_marker_frame)

    p <- ggplot(diff_marker_frame, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
        geom_point(alpha = 0.8, size = 0.6) +
        scale_colour_manual(values = c("red2", "blue2", "gray"), limits = c("ZPN", "ZRN", "None")) +
        theme(panel.grid = element_blank(), 
            panel.background = element_rect(color = "black", fill = "transparent"), 
            plot.title = element_text(hjust = 0.5)) +
        theme(legend.key = element_rect(fill = "transparent"), 
            legend.background = element_rect(fill = "transparent"), 
            legend.position = c(0.9, 0.93)) +
        geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), color = "gray", size = 0.3) +
        geom_hline(yintercept = -log(0.01, 10), color = "gray", size = 0.3) +
        xlim(-2, 2) +
        ylim(0, 10) +
        labs(x = "\nLog2 Fold Change", y = "-log10(pvalue)\n", color = "", title = "ZPN vs ZRN\n") +
        ggtitle("ZPN vs ZRN")

    up <- subset(diff_marker_frame, sig == 'ZPN')
    up <- up[order(up$p_val_adj), ][1:40, ]
    down <- subset(diff_marker_frame, sig == 'ZRN')
    down <- down[order(down$p_val_adj), ][1:40, ]

    p1 <- p + theme(legend.position = 'right') +
            geom_label_repel(data = rbind(up, down), 
                    aes(x = avg_log2FC, y = -log10(p_val_adj), label = Gene),
                    size = 3,
                    box.padding = unit(0.5, 'lines'), 
                    segment.color = 'black',
                    alpha = 0.7,
                    max.overlaps = Inf, 
                    show.legend = FALSE)

    ggsave('Figure16.png', p1, width = 5, height = 6)
    ggsave('Figure16.pdf', p1, width = 5, height = 6)

}



# Figure17: Find All Markers in Myeloid Cell types
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure17")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    marker_frame <- FindAllMarkers(myeloid_sce,group.by = "clusters", only.pos = TRUE)
    AllMakers <- 'all_markers.csv'
    marker_frame <- marker_frame %>% group_by(cluster)
    write.csv(marker_frame, file=AllMakers, quote=F)

    marker_frame <- marker_frame %>%
        group_by(cluster) %>%
        top_n(3, avg_log2FC)
    marker_frame <- marker_frame[!duplicated(marker_frame$gene),]
    #DoHeatmap(global_sce, features = marker_frame$gene, group.by = "clusters", label = TRUE)
    dp <- DotPlot(myeloid_sce,
        features = marker_frame$gene,
        group.by = "clusters"
    ) +
        coord_flip() +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        ) + scale_color_gradientn(values = seq(0, 1, 0.2), 
        colours = c("#330066", "#336699", "#66CC66", "#FFCC33")) +
        labs(x = NULL) +
        guides(size = guide_legend(order = 3))

    ggsave(
        filename = "Figure17.png",
        plot = dp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "Figure17.pdf", 
       plot = dp, 
       width = 8, 
       height = 12, 
       units = c("in"))

}



# Figure 18: Macrophage Diffusion map
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure18")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    myeloid_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macro", myeloid_sce@meta.data$clusters)]))
    myeloid_sce <- subset(myeloid_sce, clusters %in% myeloid_celltype_array)

    # Running diffusion map
    cat("Running DiffusionMap...\n")
    dm <- DiffusionMap(as.SingleCellExperiment(myeloid_sce), verbose = TRUE)
    cat("Complete\n")
    saveRDS(dm, file = "DM.rds", compress = F)

    dm <- readRDS("/P03_MyeloidCell/P999_Myeloid_Analysis/Figure18/DM.rds")
    # Ploting the diffusion map
    dc_frame <- data.frame(DC1 = eigenvectors(dm)[, 1],
                    DC2 = eigenvectors(dm)[, 2],
                    DC3 = eigenvectors(dm)[, 3],
                    DC4 = eigenvectors(dm)[, 4],
                    cluster = myeloid_sce$clusters)
    gp <- ggplot(dc_frame, aes(x = DC1, y = DC2, colour = cluster)) +
        geom_point(size = 0.5)  + 
        scale_color_manual(values = celltype_colors[myeloid_celltype_array]) +
        xlab("Diffusion component 1") + 
        ylab("Diffusion component 2") +
        theme_classic()
    ggsave(filename = "DC1_DC2.pdf", plot = gp, width = 8, height = 8, units = c("in"))
    ggsave(filename = "DC1_DC2.png", plot = gp, width = 8, height = 8, units = c("in"))

    # Plotting cell progression along the diffusion map components
    myeloid_sce$pseud_dm1 <- rank(eigenvectors(dm)[,1])      # ramyeloid cells by their dpt dm1
    myeloid_sce$pseud_dm2 <- rank(eigenvectors(dm)[,2])      # ramyeloid cells by their dpt dm2
    myeloid_sce$pseud_dm1R <- rank(-eigenvectors(dm)[,1])    # ramyeloid cells by their dpt dm1 reverse order
    myeloid_sce$pseud_dm2R <- rank(-eigenvectors(dm)[,2])    # ramyeloid cells by their dpt dm2 reverse order
    
    group_id_array <- paste0(substr(myeloid_sce@meta.data$orig.ident, 1, 2), substr(myeloid_sce@meta.data$orig.ident, 5, 5))
    myeloid_sce <- AddMetaData(myeloid_sce,
        group_id_array,
        col.name = "group_id"
    )

    myeloid_sce_exp <- as.SingleCellExperiment(myeloid_sce)

    SortedDM1 <- data.frame(
        DM1Sort = as.data.frame(colData(myeloid_sce_exp))$pseud_dm1,
        Cluster = as.data.frame(colData(myeloid_sce_exp))$clusters,
        Groups = as.data.frame(colData(myeloid_sce_exp))$group_id
    )
    SortedDM2 <- data.frame(
        DM2Sort = as.data.frame(colData(myeloid_sce_exp))$pseud_dm2,
        Cluster = as.data.frame(colData(myeloid_sce_exp))$clusters,
        Groups = as.data.frame(colData(myeloid_sce_exp))$group_id
    )
    SortedDM1R <- data.frame(
        DM1SortR = as.data.frame(colData(myeloid_sce_exp))$pseud_dm1R,
        Cluster = as.data.frame(colData(myeloid_sce_exp))$clusters,
        Groups = as.data.frame(colData(myeloid_sce_exp))$group_id
    )
    SortedDM2R <- data.frame(
        DM2SortR = as.data.frame(colData(myeloid_sce_exp))$pseud_dm2R,
        Cluster = as.data.frame(colData(myeloid_sce_exp))$clusters,
        Groups = as.data.frame(colData(myeloid_sce_exp))$group_id
    )
    # By clusters
    gp1 <- ggplot(SortedDM1, aes(x = SortedDM1[, 1], y = Cluster, fill = Cluster)) +
        geom_density_ridges(alpha = 0.5, height = 0.5) +
        geom_boxplot(width = 0.1,alpha = 0.5,outlier.alpha = 0) +
        scale_fill_manual(values = celltype_colors[myeloid_celltype_array]) +
        xlab("Diffusion component 1 (DC1)") +
        ylab("Cluster") +
        ggtitle("Cells ordered by DC1") +
        theme_bw()
    ggsave(filename = "SortedDM1.pdf", plot = gp1, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedDM1.png", plot = gp1, width = 8, height = 8, units = c("in"))

    gp2 <- ggplot(SortedDM2, aes(x = SortedDM2[, 1], y = Cluster, fill = Cluster)) +
        geom_density_ridges(alpha = 0.5, height = 0.5) +
        geom_boxplot(width = 0.1,alpha = 0.5,outlier.alpha = 0) +
        scale_fill_manual(values = celltype_colors[myeloid_celltype_array]) +
        xlab("Diffusion component 2 (DC2)") +
        ylab("Cluster") +
        ggtitle("Cells ordered by DC2") +
        theme_bw()
    ggsave(filename = "SortedDM2.pdf", plot = gp2, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedDM2.png", plot = gp2, width = 8, height = 8, units = c("in"))


    gp3 <- ggplot(SortedDM1R, aes(x = SortedDM1R[,1], y = Cluster, color = Cluster)) +
        geom_jitter() + xlab("Minus Diffusion component 1 (DC1)") + ylab("Cluster") +
        ggtitle("Cells ordered by reversed DC1")
    ggsave(filename = "SortedSampleDM1.pdf", plot = gp3, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedSampleDM1.png", plot = gp3, width = 8, height = 8, units = c("in"))


    gp4 <- ggplot(SortedDM2R, aes(x=SortedDM2R[,1], y=Cluster,color=Cluster)) +
        geom_jitter() + xlab("Minus Diffusion component 2 (DC2)") + ylab("Cluster") +
        ggtitle("Cells ordered by reversed DC2")
    ggsave(filename = "SortedSampleDM2.pdf", plot = gp4, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedSampleDM2.png", plot = gp4, width = 8, height = 8, units = c("in"))

    #By Groups

    gp1 <- ggplot(SortedDM1, aes(x = SortedDM1[, 1], y = Groups, color = Groups)) +
        geom_boxplot() +
        scale_color_manual(values = group_cols) +
        xlab("Diffusion component 1 (DC1)") +
        ylab("Groups") +
        ggtitle("Cells ordered by DC1") +
        theme_bw()
    ggsave(filename = "SortedGroupsDM1.pdf", plot = gp1, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedGroupsDM1.png", plot = gp1, width = 8, height = 8, units = c("in"))

    gp2 <- ggplot(SortedDM2, aes(x = SortedDM2[, 1], y = Groups, color = Groups)) +
        geom_boxplot() +
        scale_color_manual(values = group_cols) +
        xlab("Diffusion component 2 (DC2)") +
        ylab("Groups") +
        ggtitle("Cells ordered by DC2") +
        theme_bw()
    ggsave(filename = "SortedGroupsDM2.pdf", plot = gp2, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedGroupsDM2.png", plot = gp2, width = 8, height = 8, units = c("in"))


    gp3 <- ggplot(SortedDM1R, aes(x = SortedDM1R[,1], y = Groups, color = Groups)) +
        geom_jitter() + xlab("Minus Diffusion component 1 (DC1)") + ylab("Groups") +
        ggtitle("Cells ordered by reversed DC1")
    ggsave(filename = "SortedGroupsDM1-1.pdf", plot = gp3, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedGroupsDM1-1.png", plot = gp3, width = 8, height = 8, units = c("in"))


    gp4 <- ggplot(SortedDM2R, aes(x=SortedDM2R[,1], y=Groups,color=Groups)) +
        geom_jitter() + xlab("Minus Diffusion component 2 (DC2)") + ylab("Groups") +
        ggtitle("Cells ordered by reversed DC2")
    ggsave(filename = "SortedGroupsDM2-1.pdf", plot = gp4, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedGroupsDM2-1.png", plot = gp4, width = 8, height = 8, units = c("in"))
 



    if (F) {
        # Make and interactive 2D and 3D diffusion map figure
        library(plotly)

        # interactive 2D
        p2 = plot_ly(
            x = tmp$DC1,
            y = tmp$DC2, 
            type = "scatter", 
            mode = "markers", 
            color = tmp$Samples, 
            marker.size = 0.5
        )
        htmlwidgets::saveWidget(as_widget(p2), "Interactive2D_DiffM.html", title = "Diffusion map")

        #interactive 3D
        p = plot_ly(x=tmp$DC1, y=tmp$DC2, z=tmp$DC3, type="scatter3d", mode="markers", color=tmp$Samples, marker = list(size = 2 ))
        htmlwidgets::saveWidget(as_widget(p), "Interactive3D.html", title = "Diffusion map")
    }
}




# Figure 19: Monocytes Macrophage Diffusion map
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure19")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    myeloid_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Mono|Macro", myeloid_sce@meta.data$clusters)]))
    myeloid_sce <- subset(myeloid_sce, clusters %in% myeloid_celltype_array)

    # Running diffusion map
    cat("Running DiffusionMap...\n")
    myeloid_sce_exp <- as.SingleCellExperiment(myeloid_sce)
    dm <- DiffusionMap(as.SingleCellExperiment(myeloid_sce), verbose = TRUE)
    cat("Complete\n")
    saveRDS(dm, file = "Monocyte_Macrophage_DM.rds", compress = F)

    dm <- readRDS("/P02_myeloidCell/P999_myeloid_Analysis/Figure20/DM.rds")
    # Ploting the diffusion map
    dc_frame <- data.frame(DC1 = eigenvectors(dm)[, 1],
                    DC2 = eigenvectors(dm)[, 2],
                    DC3 = eigenvectors(dm)[, 3],
                    DC4 = eigenvectors(dm)[, 4],
                    cluster = myeloid_sce$clusters)
    gp <- ggplot(dc_frame, aes(x = DC1, y = DC2, colour = cluster)) +
        geom_point(size = 0.5)  + 
        scale_color_manual(values = celltype_colors[myeloid_celltype_array]) +
        xlab("Diffusion component 1") + 
        ylab("Diffusion component 2") +
        theme_classic()
    ggsave(filename = "DC1_DC2.pdf", plot = gp, width = 8, height = 8, units = c("in"))
    ggsave(filename = "DC1_DC2.png", plot = gp, width = 8, height = 8, units = c("in"))

    # Plotting cell progression along the diffusion map components
    myeloid_sce$pseud_dm1 <- ramyeloid(eigenvectors(dm)[,1])      # ramyeloid cells by their dpt dm1
    myeloid_sce$pseud_dm2 <- ramyeloid(eigenvectors(dm)[,2])      # ramyeloid cells by their dpt dm2
    myeloid_sce$pseud_dm1R <- ramyeloid(-eigenvectors(dm)[,1])    # ramyeloid cells by their dpt dm1 reverse order
    myeloid_sce$pseud_dm2R <- ramyeloid(-eigenvectors(dm)[,2])    # ramyeloid cells by their dpt dm2 reverse order

    SortedDM1 <- data.frame(DM1Sort = as.data.frame(colData(myeloid_sce_exp))$pseud_dm1,
                            Samples = as.data.frame(colData(myeloid_sce_exp))$ident)
    SortedDM2 <- data.frame(DM2Sort = as.data.frame(colData(myeloid_sce_exp))$pseud_dm2,
                            Samples = as.data.frame(colData(myeloid_sce_exp))$ident)
    SortedDM1R <- data.frame(DM1SortR = as.data.frame(colData(myeloid_sce_exp))$pseud_dm1R,
                            Samples = as.data.frame(colData(myeloid_sce_exp))$ident)
    SortedDM2R <- data.frame(DM2SortR = as.data.frame(colData(myeloid_sce_exp))$pseud_dm2R,
                            Samples = as.data.frame(colData(myeloid_sce_exp))$ident)

    gp1 <- ggplot(SortedDM1, aes(x = SortedDM1[, 1], y = Samples, color = Samples)) +
        geom_boxplot() +
        scale_color_manual(values = celltype_colors[myeloid_celltype_array]) +
        xlab("Diffusion component 1 (DC1)") +
        ylab("Samples") +
        ggtitle("Cells ordered by DC1") +
        theme_bw()
    ggsave(filename = "SortedDM1.pdf", plot = gp1, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedDM1.png", plot = gp1, width = 8, height = 8, units = c("in"))

    gp2 <- ggplot(SortedDM2, aes(x = SortedDM2[, 1], y = Samples, color = Samples)) +
        geom_boxplot() +
        scale_color_manual(values = celltype_colors[myeloid_celltype_array]) +
        xlab("Diffusion component 2 (DC2)") +
        ylab("Samples") +
        ggtitle("Cells ordered by DC2") +
        theme_bw()
    ggsave(filename = "SortedDM2.pdf", plot = gp2, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedDM2.png", plot = gp2, width = 8, height = 8, units = c("in"))

    
    ggplot(SortedDM1R, aes(x=SortedDM1R[,1], y=Samples,color=Samples)) +
        geom_jitter() + xlab("Minus Diffusion component 1 (DC1)") + ylab("Samples") +
        ggtitle("Cells ordered by reversed DC1")
        ggplot(SortedDM2R, aes(x=SortedDM2R[,1], y=Samples,color=Samples)) +
        geom_jitter() + xlab("Minus Diffusion component 2 (DC2)") + ylab("Samples") +
        ggtitle("Cells ordered by reversed DC2")

    if (F) {
        # Make and interactive 2D and 3D diffusion map figure
        library(plotly)

        # interactive 2D
        p2 = plot_ly(
            x = tmp$DC1,
            y = tmp$DC2, 
            type = "scatter", 
            mode = "markers", 
            color = tmp$Samples, 
            marker.size = 0.5
        )
        htmlwidgets::saveWidget(as_widget(p2), "Interactive2D_DiffM.html", title = "Diffusion map")

        #interactive 3D
        p = plot_ly(x=tmp$DC1, y=tmp$DC2, z=tmp$DC3, type="scatter3d", mode="markers", color=tmp$Samples, marker = list(size = 2 ))
        htmlwidgets::saveWidget(as_widget(p), "Interactive3D.html", title = "Diffusion map")
    }
}




# Figure 20 Marker of Monocyte cell type VlnPlot by group
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure20")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Mono", myeloid_sce@meta.data$clusters)]))

    cd8_sce <- subset(myeloid_sce, clusters %in% cd8_celltype_array)
    markers <- read.table(
        file = "/P02_myeloidCell/P999_myeloid_Analysis/Figure16/markers.txt",
        header = F, sep = "\t"
    )[, "V1"]

    #orig.ident
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
        compaired <- list(
            c("ZPT", "ZRT"),
            c("ZPN", "ZRN"),
            c("ZPT", "ZPN"),
            c("ZRT", "ZRN")
        )

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'group', log=TRUE,... ) +
            xlab("") + 
            ylab(feature) + 
            ggtitle("") +
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
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 3)
        return(p)
    }

    p <- StackedVlnPlot(cd8_sce,
        markers,
        pt.size = 0, 
        cols = group_cols
    )

    ggsave(filename = "Figure20.png", plot = p, width = 8, height = 18, units = c("in"))
    ggsave(filename = "Figure20.pdf", plot = p, width = 8, height = 18, units = c("in"))

}




# Figure 21 Macrophage Vocalno Plot by ZPT and ZRT

if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure21")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macrophage", myeloid_sce@meta.data$clusters)]))

    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)
    group_id_array <- paste0(substr(macrophage_sce@meta.data$orig.ident, 1, 2), substr(macrophage_sce@meta.data$orig.ident, 5, 5))
    macrophage_sce <- AddMetaData(macrophage_sce,
        group_id_array,
        col.name = "group_id"
    )

    macrophage_tumor_sce <- subset(macrophage_sce, group_id %in% c("ZPT", "ZRT"))

    if (F) {
        allCells = colnames(macrophage_tumor_sce)
        allType = names(table(macrophage_tumor_sce@meta.data$orig.ident))

        choose_Cells = unlist(lapply(allType, function(x) {
            cgCells = allCells[macrophage_tumor_sce@meta.data$orig.ident == x]
            num = ceiling(dim(as.data.frame(cgCells))[1] * 0.6)
            cg = sample(cgCells, num)
            cg
        }))

        macrophage_tumor_sce = macrophage_tumor_sce[, allCells %in% choose_Cells]
    }
    options(future.globals.maxSize = 891289600)
    deg <- FindMarkers(macrophage_tumor_sce,
        ident.1 = "ZPT", 
        group.by = "group_id", 
        assay = "RNA", 
        slot = "counts", 
        logfc.threshold = 0, 
        min.pct = 0
    )
    
    saveRDS(deg, file = "ZPT_ZRT_macrophage_deg.rds", compress = F)
    deg <- readRDS("ZPT_ZRT_macrophage_deg.rds")
    if (T) {
        log2FC_threshold = 0.8
        p_value_threshold = 0.01
        deg[which(deg$p_val_adj < 0.05 & deg$avg_log2FC <= -log2FC_threshold),'sig'] <- 'LRT'
        deg[which(deg$p_val_adj < 0.05 & deg$avg_log2FC >= log2FC_threshold),'sig'] <- 'PHT'
        deg[which(deg$p_val_adj >= 0.05 | abs(deg$avg_log2FC) < log2FC_threshold),'sig'] <- 'None'


        # 横轴 log2FC，纵轴 -log10(adj.P.Val)，颜色表示差异
        deg$Gene <- row.names(deg)
        deg$percent_difference <- deg$pct.1 - deg$pct.2
        p <- ggplot(deg, aes(x = percent_difference, y = avg_log2FC, color = sig)) +
            geom_point(alpha = 0.8, size = 0.6) +
            scale_colour_manual(values = c("#00af3e", "#CD1E24", "gray"), limits = c("PHT", "LRT", "None")) +
            theme(panel.grid = element_blank(), 
                panel.background = element_rect(color = "black", fill = "transparent"), 
                plot.title = element_text(hjust = 0.5)) +
            theme(legend.key = element_rect(fill = "transparent"), 
                legend.background = element_rect(fill = "transparent"), 
                legend.position = c(0.9, 0.93)) +
            xlim(-1, 1) +
            ylim(-2, 2) +
            labs(
                x = "\nPercentage Difference",
                y = "Log2Fold Change\n", 
                color = "", 
                title = "Macrophage ZPT vs ZRT\n"
            )

        up <- subset(deg, sig == 'PHT')
        up <- up[order(up$p_val_adj), ][1:40, ]
        down <- subset(deg, sig == 'LRT')
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

        ggsave('PHT_LRT_VolcanoPlot.png', p1, width = 6 * 2, height = 4 * 2)
        ggsave('PHT_LRT_VolcanoPlot.pdf', p1, width = 6 * 2, height = 4 * 2)

    }

    mdb_h <- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    # GSEA BarPlot
    deg$genes = rownames(deg)
    cluster0.genes<- deg %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
    ranks<- deframe(cluster0.genes)

    fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

    #png("TP53.png",width = 1080,height = 768)
    #plotEnrichment(fgsea_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) +
    #    labs(title = "HALLMARK_ALLOGRAFT_REJECTION")
    #dev.off()

    gp <- ggplot(fgseaRes %>% 
        as_tibble() %>% 
        arrange(desc(NES)) %>% 
        filter(pval < 0.05) %>% 
        head(n= 60), aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = NES)) +
            coord_flip() +
            labs(x = "Hallmark", y = "Normalized Enrichment Score", 
            title = "HALLMARK gene sets NES from GSEA") +
            theme_bw()

    ggsave(filename = "macrophage_C5_gsea_barplot.png", width = 12, height = 12, units = c("in"))
    ggsave(filename = "macrophage_C5_gsea_barplot.pdf", width = 12, height = 12, units = c("in"))


    mdb_h <- msigdbr(species = "Homo sapiens", category = "C2")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    # GSEA BarPlot
    deg$genes = rownames(deg)
    cluster0.genes<- deg %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
    ranks<- deframe(cluster0.genes)

    fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)

    #png("TP53.png",width = 1080,height = 768)
    #plotEnrichment(fgsea_sets[[head(fgseaRes[order(pval), ], 1)$pathway]], ranks) +
    #    labs(title = "HALLMARK_ALLOGRAFT_REJECTION")
    #dev.off()

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

    ggsave(filename = "macrophage_C5_gsea_barplot.png", width = 12, height = 6, units = c("in"))
    ggsave(filename = "macrophage_C5_gsea_barplot.pdf", width = 12, height = 6, units = c("in"))



}



# Figure 22 Monocytes Vocalno Plot by ZPT and ZRT
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure22")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    monocyte_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Monocyte", myeloid_sce@meta.data$clusters)]))

    monocyte_sce <- subset(myeloid_sce, clusters %in% monocyte_celltype_array)
    group_id_array <- paste0(substr(monocyte_sce@meta.data$orig.ident, 1, 2), substr(monocyte_sce@meta.data$orig.ident, 5, 5))
    monocyte_sce <- AddMetaData(monocyte_sce,
        group_id_array,
        col.name = "group_id"
    )

    monocyte_tumor_sce <- subset(monocyte_sce, group_id %in% c("ZPT", "ZRT"))

    if (F) {
        allCells = colnames(monocyte_tumor_sce)
        allType = names(table(monocyte_tumor_sce@meta.data$orig.ident))

        choose_Cells = unlist(lapply(allType, function(x) {
            cgCells = allCells[monocyte_tumor_sce@meta.data$orig.ident == x]
            num = ceiling(dim(as.data.frame(cgCells))[1] * 0.6)
            cg = sample(cgCells, num)
            cg
        }))

        monocyte_tumor_sce = monocyte_tumor_sce[, allCells %in% choose_Cells]
    }

    options(future.globals.maxSize = 891289600)
    deg <- FindMarkers(monocyte_tumor_sce,
        ident.1 = "ZPT", 
        group.by = "group_id", 
        assay = "RNA", 
        slot = "counts", 
        logfc.threshold = 0, 
        min.pct = 0
    )

    if (F) {
        log2FC_threshold = 0.8
        p_value_threshold = 0.01
        deg[which(deg$p_val_adj < 0.05 & deg$avg_log2FC <= -log2FC_threshold),'sig'] <- 'LRT'
        deg[which(deg$p_val_adj < 0.05 & deg$avg_log2FC >= log2FC_threshold),'sig'] <- 'PHT'
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
            xlim(-1, 1) +
            ylim(-2, 2) +
            labs(
                x = "\nPercentage Difference",
                y = "Log2Fold Change\n", 
                color = "", 
                title = "Monocyte ZPT vs ZRT\n"
            )

        up <- subset(deg, sig == 'ZPT')
        up <- up[order(up$p_val_adj), ][1:40, ]
        down <- subset(deg, sig == 'ZRT')
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

        ggsave('ZPT_ZRT_VolcanoPlot.png', p1, width = 6 * 2, height = 4 * 2)
        ggsave('ZPT_ZRT_VolcanoPlot.pdf', p1, width = 6 * 2, height = 4 * 2)


    }


    mdb_h <- msigdbr(species = "Homo sapiens", category = "C7")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    # GSEA BarPlot
    deg$genes = rownames(deg)
    cluster0.genes<- deg %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
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

    ggsave(filename = "monocyte_C7_gsea_barplot.png", width = 12, height = 6, units = c("in"))
    ggsave(filename = "monocyte_C7_gsea_barplot.pdf", width = 12, height = 6, units = c("in"))


    mdb_h <- msigdbr(species = "Homo sapiens", category = "C5")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    # GSEA BarPlot
    deg$genes = rownames(deg)
    cluster0.genes<- deg %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
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

    ggsave(filename = "monocyte_C5_gsea_barplot.png", width = 12, height = 6, units = c("in"))
    ggsave(filename = "monocyte_C5_gsea_barplot.pdf", width = 12, height = 6, units = c("in"))


}



# Figure 23 DC pyscenic analysis by cell types

if (F) {
    library(SCENIC)
    library(SCopeLoomR)
    library(pheatmap) 

    fig_outdir <- paste0(outdir, "/", "Figure23")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    dc_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("DC", myeloid_sce@meta.data$clusters)]))
    dc_sce <- subset(myeloid_sce, clusters %in% dc_celltype_array)

    #dc_matrix <-dc_sce@assays$RNA@counts[(rowSums(dc_sce@assays$RNA@counts == 0) / dim(dc_sce@assays$RNA@counts)[2] < 0.99),]

    write.csv(t(as.matrix(dc_sce@assays$RNA@counts)), file = "dc.csv")

    # run00-03.sh

    scenicLoomPath='dc_SCENIC.loom'
    loom <- open_loom(scenicLoomPath)
    # Read information from loom file: 
    regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
    regulons <- regulonsToGeneLists(regulons_incidMat)
    regulonAUC <- get_regulons_AUC(loom, column.attr.name="RegulonsAUC")

    regulonAucThresholds <- get_regulon_thresholds(loom)
    tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

    embeddings <- get_embeddings(loom)  
    close_loom(loom)

    rownames(regulonAUC)
    names(regulons)

    
    sub_regulonAUC <- regulonAUC[, match(colnames(dc_sce), colnames(regulonAUC))]
    identical(colnames(sub_regulonAUC), colnames(dc_sce))

    cellClusters <- data.frame(row.names = colnames(dc_sce), 
                           seurat_clusters = as.character(dc_sce$seurat_clusters))
    cellTypes <- data.frame(row.names = colnames(dc_sce), 
                            celltype = dc_sce$clusters)
    head(cellTypes)
    head(cellClusters)
    sub_regulonAUC[1:4,1:4] 



    rss <- calcRSS(
        AUC = getAUC(regulonAUC),
        cellAnnotation = cellTypes[
            colnames(sub_regulonAUC),
            selectedResolution
        ]
    )




































    cellTypes <- as.data.frame(dc_sce$clusters)
    colnames(cellTypes) <- "celltype"
    selectedResolution <- "celltype"
    cellTypes$celltype <- as.character(cellTypes$celltype)

    cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution]) 

    sub_regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ] # remove extended regulons
    dim(sub_regulonAUC)

    regulonActivity_byGroup <- sapply(cellsPerGroup, function(cells) {
        rowMeans(getAUC(sub_regulonAUC[,cells]))
    })

    regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale = T))

    rss <- regulonActivity_byGroup_Scaled

    rss <- rbind(
        rss[order(rss[, 1], decreasing = T), ][1:10, ],
        rss[order(rss[, 2], decreasing = T), ][1:10, ],
        rss[order(rss[, 3], decreasing = T), ][1:10, ])



    df = do.call(rbind, lapply(1:ncol(rss), function(i) {
        dat = data.frame(
            path = rownames(rss),
            cluster = colnames(rss)[i],
            sd.1 = rss[,i],
            sd.2 = apply(rss[,-i],1,median)
        )
    }))

    df$fc = df$sd.1 - df$sd.2
    top5 <- df %>%
        group_by(cluster) %>%
        top_n(10, fc)
    rowcn = data.frame(path = top5$cluster)
    n = rss[top5$path, ]

    n1 <- regulonActivity_byGroup[top5$path,]

    ht <- Heatmap(n,
        rect_gp = gpar(col = "white", lwd = 2),
        heatmap_legend_param = list(title = "", 
                        title_gp = gpar(fontsize = 8), 
                        labels_gp = gpar(fontsize = 8))
    )
    png("RSSHeatmap.png",width = 768,height = 1080)
    draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

    pdf("RSSHeatmap.pdf",width = 6,height = 10)
    draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()



    rss <- calcRSS(
        AUC = getAUC(regulonAUC),
        cellAnnotation = cellTypes[
            colnames(regulonAUC),
            selectedResolution
        ]
    )
    png("RegulonSpecificityScoreDotPlot.png",width = 1080,height = 1080)
    rss=na.omit(rss) 
    rssPlot <- plotRSS(rss)
    rssPlot$plot
    dev.off()

    png("RegulonSpecificityScoreHeatmap.png", width = 1080, height = 1080)
    scale_rss <- scale(rss)[order(scale(rss)[, "DC-C1"]), ]
    Heatmap(scale(rss))
    dev.off()

    library(dplyr)
    rss 


}




# Figure 24 Marker of DC cell type VlnPlot by group
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure24")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    samples_list <- read.table(sample_list, header = F, sep = "\t")
    group_array <- samples_list[, "V2"]
    names(group_array) <- samples_list[, "V1"]
    myeloid_sce$group <- group_array[myeloid_sce$orig.ident]

    dc_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("DC-C1", myeloid_sce@meta.data$clusters)]))

    dc_sce <- subset(myeloid_sce, clusters %in% dc_celltype_array)
    markers_array <- c("HLA-A","HLA-B","HLA-C","BIRC3","HLA-DPB1","HLA-DQB2")

    #orig.ident
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
        compaired <- list(
            c("ZPT", "ZRT"),
            c("ZPN", "ZRN"),
            c("ZPT", "ZPN"),
            c("ZRT", "ZRN")
        )

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'group', log=TRUE,... ) +
            xlab("") + 
            ylab(feature) + 
            ggtitle("") +
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
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 3)
        return(p)
    }

    p <- StackedVlnPlot(dc_sce,
        markers_array,
        pt.size = 0, 
        cols = group_cols
    )

    ggsave(filename = "Figure24.png", plot = p, width = 8, height = 8, units = c("in"))
    ggsave(filename = "Figure24.pdf", plot = p, width = 8, height = 8, units = c("in"))

}





# Figure 25: Plot myeloid Marker Gene Heatmap By Celltype
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure25")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    cell_marker <- "
		MonoIV	CKB
		MonoIV	LILRA5
		MonoIV	FAM65B
		MonoIV	LIMD2
		MonoIV	RHOC
		MonoIV	HES4
		MonoIV	LILRB2
		MonoIV	POU2F2
		MonoIV	MTSS1
		MonoIV	LILRA1
		MonoIV	TESC
		MonoIV	CD79B
		MonoIII	IFITM2
		MonoIII	CDKN1C
		MonoIII	LYPD2
		MonoIII	PTX3
		MonoIII	THBS1
		MonoIII	APOBEC3A
		MonoIII	PTGS2
		MonoIII	NLRP3
		MonoIII	IL1B
		MonoIII	CD300E
		MonoIII	AREG
		MonoIII	G0S2
		MonoIII	FCN1
		MonoIII	EREG
		MonoIII	VCAN
		MonoIII	SERPINB2
		MonoIII	S100A12
		MacrophageII	RARRES1
		MacrophageII	LGMN
		MacrophageII	SLC40A1
		MacrophageII	CXCL9
		MacrophageII	PLA2G7
		MacrophageII	CCL13
		MacrophageII	MMP9
		MacrophageII	RNASE1
		MacrophageII	CCL2
		MacrophageII	PTGDS
		MacrophageII	CHIT1
		MacrophageII	CHI3L1
		MacrophageII	SEPP1
		MacrophageII	SPP1
		MacrophageII	PPBP
		MacrophageI	AKR1C3
		MacrophageI	PCOLCE2
		MacrophageI	RND3
		MacrophageI	CAMP
		MacrophageI	DEFB1
		MacrophageI	IGFBP2
		MacrophageI	FHL1
		MacrophageI	MME
		MacrophageI	FOLR3
		MacrophageI	INHBA
		MacrophageI	GPD1
		MacrophageI	HP
		MacrophageI	MARCO
		MacrophageI	MRC1
		MacrophageI	PLBD1
		MacrophageI	FABP4"

    Idents(myeloid_sce) <- "seurat_clusters"
    cell_marker_frame = marker_frame <- read.table(
			text = gsub("\t\t","",cell_marker),
			sep = "\t", 
			col.names = c("celltype", "gene")
		)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$gene),]

    AverageExp <- AverageExpression(myeloid_sce,
        features = cell_marker_frame$gene, 
        group.by = "clusters"
    )

    expr <- AverageExp$RNA

    heatmapColorRamp = colorRampPalette(c("#0658D6","#FDF8EE","#ED7A03"))(64)

    expr <- ScaleData(expr)
    color_array <- c("#949483",
                    "#F47B7B",
                    "#9F1F5C",
                    "#EF9020",
                    "#00AF3E",
                    "#85B7E2",
                    "#29245C",
                    "#FFD616",
                    "#E5352B",
                    "#FFD616",
                    "#E990AB",
                    "#0081B4",
                    "#96CBB3",
                    "#91BE3E",
                    "#39A6DD",
                    "#EB0973",
                    "#DDE2E0",
                    "#333C41")

    ordered_celltype_array <- sort(colnames(expr))

    celltype_color_structure <- structure(color_array[1:(length(unique(cell_marker_frame$celltype)))],
                                        names = unique(cell_marker_frame$celltype)) 

    ha = rowAnnotation(labels = cell_marker_frame$celltype,
                    col = list(bar = celltype_color_structure))
    expr <- expr[, ordered_celltype_array]
    ht <- Heatmap(as.matrix(expr),
            col = heatmapColorRamp,
            right_annotation = ha,
            rect_gp = gpar(col = "white", lwd = 1),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 16),
            row_split = cell_marker_frame$celltype,
            row_title_side = "right",
            row_title_rot = 0, 
            row_title_gp = gpar(fontsize = 25),
            cluster_columns = F,
            cluster_rows = F)

    png("Figure25.png",width = 1080,height = 1080 * 1.5)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

    pdf("Figure25.pdf",width = 12,height = 16)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

}




#VIP Figure 6C and Figure 6D Macrophage signature VlnPlot by group
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure26")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    sample_list <- "
        ZP01N	PHN
        ZP01T	PHT
        ZP02N	PHN
        ZP02T	PHT
        ZP03N	PHN
        ZP03T	PHT
        ZP04N	PHN
        ZP04T	PHT
        ZP05N	PHN
        ZP05T	PHT
        ZP06N	PHN
        ZP06T	PHT
        ZP07N	PHN
        ZP07T	PHT
        ZP08N	PHN
        ZP08T	PHT
        ZP09N	PHN
        ZP09T	PHT
        ZP10N	PHN
        ZP10T	PHT
        ZP11N	PHN
        ZP11T	PHT
        ZR01N	LRN
        ZR01T	LRT
        ZR02N	LRN
        ZR02T	LRT
        ZR03N	LRN
        ZR03T	LRT
        ZR04N	LRN
        ZR04T	LRT
        ZR05N	LRN
        ZR05T	LRT
        ZR06N	LRN
        ZR06T	LRT
        ZR07N	LRN
        ZR07T	LRT"
    
    samples_list_frame <- read.table(
        text = gsub("        ","",sample_list),
        sep = "\t", 
        col.names = c("V1", "V2")
    )
    group_array <- samples_list_frame$V2
    names(group_array) <- samples_list_frame$V1
    myeloid_sce$group <- group_array[myeloid_sce$orig.ident]
   
    M1M2_signature_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/M1M2_signature.txt"

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macrophage", myeloid_sce@meta.data$clusters)]))

    M1M2_signature_frame <- read.table(file = M1M2_signature_file, sep = "\t",header = F,stringsAsFactors = F)
    M1M2_signature_frame <- M1M2_signature_frame %>% filter(V1 %in% c("M2","Monocyte"))
    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)
    macrophage_sce <- subset(macrophage_sce,group %in% c("PHT","LRT"))
    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(M1M2_signature_frame, M1M2_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        macrophage_sce <- AddModuleScore(
                        object = macrophage_sce,
                        features = list(gene_array),
                        name = signature)
		colnames(macrophage_sce@meta.data) <- gsub(colnames(macrophage_sce@meta.data),
      		pattern = paste0(signature, 1),
			replacement = signature
		)
	}
    
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
        compaired <- list(
            c("LRT", "PHT")
        )

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'group', log=TRUE,... ) +
            xlab("") + 
            ylab(feature) + 
            ggtitle("") +
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
                #test.args = "less", 
                test = wilcox.test
            ) + geom_boxplot(width = 0.3, fill = "white")
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.25, 0, -0.25, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
                plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
                theme(axis.text.x=element_text(), axis.ticks.x = element_line())
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 2)
        return(p)
    }

    p <- StackedVlnPlot(macrophage_sce,
        names(split_signature_list),
        pt.size = 0, 
        cols = group_cols
    )


    ggsave(filename = "MacrophageSignatureVlnPlot.png", plot = p, width = 10, height = 6, units = c("in"))
    ggsave(filename = "MacrophageSignatureVlnPlot.pdf", plot = p, width = 10, height = 6, units = c("in"))



    M1M2_signature_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/M1M2_signature.txt"

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Monocyte", myeloid_sce@meta.data$clusters)]))

    M1M2_signature_frame <- read.table(file = M1M2_signature_file, sep = "\t",header = F,stringsAsFactors = F)
    M1M2_signature_frame <- M1M2_signature_frame %>% filter(V1 %in% c("M2","Monocyte"))
    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)
    macrophage_sce <- subset(macrophage_sce,group %in% c("PHT","LRT"))
    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(M1M2_signature_frame, M1M2_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        macrophage_sce <- AddModuleScore(
                        object = macrophage_sce,
                        features = list(gene_array),
                        name = signature)
		colnames(macrophage_sce@meta.data) <- gsub(colnames(macrophage_sce@meta.data),
      		pattern = paste0(signature, 1),
			replacement = signature
		)
	}
    
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
        compaired <- list(
            c("LRT", "PHT")
        )

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'group', log=TRUE,... ) +
            xlab("") + 
            ylab(feature) + 
            ggtitle("") +
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
                #test.args = "less", 
                test = wilcox.test
            ) + geom_boxplot(width = 0.3, fill = "white")
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.25, 0, -0.25, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
                plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
                theme(axis.text.x=element_text(), axis.ticks.x = element_line())
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 2)
        return(p)
    }

    p <- StackedVlnPlot(macrophage_sce,
        names(split_signature_list),
        pt.size = 0, 
        cols = group_cols
    )


    ggsave(filename = "MonocyteSignatureVlnPlot.png", plot = p, width = 10, height = 6, units = c("in"))
    ggsave(filename = "MonocyteSignatureVlnPlot.pdf", plot = p, width = 10, height = 6, units = c("in"))



}



# Figure 27 Myeloid subpopulation cell type pysenic analysis
if (F) {
    library(SCENIC)
    library(SCopeLoomR)
    library(pheatmap)

    fig_outdir <- paste0(outdir, "/", "Figure27")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    # myeloid_matrix <- myeloid_sce@assays$RNA@counts[(rowSums(myeloid_sce@assays$RNA@counts == 0) / dim(myeloid_sce@assays$RNA@counts)[2] < 0.99),]

    write.csv(t(as.matrix(myeloid_sce@assays$RNA@counts)), file = "myeloid.csv")

    # run00-03.sh

    scenicLoomPath = "myeloid_SCENIC.loom"
    loom <- open_loom(scenicLoomPath)
    # Read information from loom file:
    regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
    regulons <- regulonsToGeneLists(regulons_incidMat)
    regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

    regulonAucThresholds <- get_regulon_thresholds(loom)
    tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

    embeddings <- get_embeddings(loom)
    close_loom(loom)

    rownames(regulonAUC)
    names(regulons)


    sub_regulonAUC <- regulonAUC[, match(colnames(dc_sce), colnames(regulonAUC))]
    identical(colnames(sub_regulonAUC), colnames(dc_sce))

    cellClusters <- data.frame(
        row.names = colnames(dc_sce),
        seurat_clusters = as.character(dc_sce$seurat_clusters)
    )
    cellTypes <- data.frame(
        row.names = colnames(dc_sce),
        celltype = dc_sce$clusters
    )
    head(cellTypes)
    head(cellClusters)
    sub_regulonAUC[1:4, 1:4]



    rss <- calcRSS(
        AUC = getAUC(regulonAUC),
        cellAnnotation = cellTypes[
            colnames(sub_regulonAUC),
            selectedResolution
        ]
    )
}






# Figure 28 Macrophage signature correlated with CD8 exhausted
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure28")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    sample_list <- "
        ZP01N	PHN
        ZP01T	PHT
        ZP02N	PHN
        ZP02T	PHT
        ZP03N	PHN
        ZP03T	PHT
        ZP04N	PHN
        ZP04T	PHT
        ZP05N	PHN
        ZP05T	PHT
        ZP06N	PHN
        ZP06T	PHT
        ZP07N	PHN
        ZP07T	PHT
        ZP08N	PHN
        ZP08T	PHT
        ZP09N	PHN
        ZP09T	PHT
        ZP10N	PHN
        ZP10T	PHT
        ZP11N	PHN
        ZP11T	PHT
        ZR01N	LRN
        ZR01T	LRT
        ZR02N	LRN
        ZR02T	LRT
        ZR03N	LRN
        ZR03T	LRT
        ZR04N	LRN
        ZR04T	LRT
        ZR05N	LRN
        ZR05T	LRT
        ZR06N	LRN
        ZR06T	LRT
        ZR07N	LRN
        ZR07T	LRT"
    
    samples_list_frame <- read.table(
        text = gsub("        ","",sample_list),
        sep = "\t", 
        col.names = c("V1", "V2")
    )
    group_array <- samples_list_frame$V2
    names(group_array) <- samples_list_frame$V1
    myeloid_sce$group <- group_array[myeloid_sce$orig.ident]
   
    M1M2_signature_file = "/P03_MyeloidCell/P999_Myeloid_Analysis/data/M1M2_signature.txt"

    macrophage_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macrophage", myeloid_sce@meta.data$clusters)]))

    M1M2_signature_frame <- read.table(file = M1M2_signature_file, sep = "\t",header = F,stringsAsFactors = F)
    
    macrophage_sce <- subset(myeloid_sce, clusters %in% macrophage_celltype_array)

    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(M1M2_signature_frame, M1M2_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        macrophage_sce <- AddModuleScore(
                        object = macrophage_sce,
                        features = list(gene_array),
                        name = signature)
		colnames(macrophage_sce@meta.data) <- gsub(colnames(macrophage_sce@meta.data),
      		pattern = paste0(signature, 1),
			replacement = signature
		)
	}
    
    sample_M2_frame <- macrophage_sce@meta.data %>%
        group_by(orig.ident) %>%
        summarise(mean(M2)) %>%
        filter(grepl("T", orig.ident)) %>%
        rename("sample_id" = "orig.ident", "M2" = "mean(M2)")

    
    # CD8 exhausted cell proprotion

    cell_stat_frame <- as.data.frame.array(table(
        as.character(tnk_sce$orig.ident),
        as.character(tnk_sce$clusters)
    ))

    celltype_array = colnames(cell_stat_frame)

    cell_stat_frame$Total_cell <- rowSums(cell_stat_frame)

    cell_stat_frame$Sample <- row.names(cell_stat_frame)
    row.names(cell_stat_frame) <- NULL

    cell_stat_frame$sample_group <- paste0(substr(cell_stat_frame$Sample,1,2),
                                        substr(cell_stat_frame$Sample,5,6))

    for(celltype in celltype_array){
        cell_stat_frame[celltype] = cell_stat_frame[celltype] / cell_stat_frame["Total_cell"]
    }
    cd8_celltype_array = colnames(cell_stat_frame %>% dplyr::select(starts_with("CD8")))
    cell_stat_frame <- cell_stat_frame[, c(cd8_celltype_array, "Sample", "sample_group")]

    exhausted_stat_frame <- select(cell_stat_frame, c("Sample", "CD8-C4")) %>% rename(exhausted = "CD8-C4")



    # CSF1 expression vs M2

    csf1_frame <- tnk_sce[["RNA"]]@scale.data["CSF1", ]
    


    M2_exh_frame <- merge(sample_M2_frame, exhausted_stat_frame,
        by.x = c("sample_id"), 
        by.y = "Sample"
    )
    cor <- cor.test(M2_exh_frame$M2,M2_exh_frame$exhausted)
    gp <- ggplot(M2_exh_frame, aes(x = M2, y = exhausted)) +
        geom_point(col = "blue", size = 2, alpha = .05, shape = 19) +
        geom_smooth(method = lm, colour = "#D25565", fill = "#fe5f55", size = 1) +
        theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            axis.text = element_blank()
        ) +
        ggtitle(paste0("estimate=",cor$estimate," p=",cor$p.value))

    ggsave("Figure28.png",plot = gp,width = 8, height = 8)
    ggsave("Figure28.pdf",plot = gp,width = 8, height = 8)

}




# Figure29: velocyto analysis
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure29")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    library(SeuratWrappers)
    library(stringr)
    library(ggplot2)
    library(Seurat)
    library(dior)

    mono_celltype_array <- unique(as.character(myeloid_sce@meta.data$clusters[grepl("Macro|Mono", myeloid_sce@meta.data$clusters)]))

    myeloid_sce <- subset(myeloid_sce, clusters %in% mono_celltype_array)
    write_h5(myeloid_sce, 
        file = "mono_macro.h5", 
        object.type = 'seurat', 
        assay.name = 'RNA', 
        save.graphs = TRUE, 
        save.scale=FALSE)

}




# Figure30: monocle3 analysis
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure30")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    ##创建CDS对象并预处理数据
    data <- GetAssayData(myeloid_sce, assay = 'RNA', slot = 'counts')
    cell_metadata <- myeloid_sce@meta.data
    gene_annotation <- data.frame(gene_short_name = rownames(data))
    rownames(gene_annotation) <- rownames(data)
    cds <- new_cell_data_set(data,
                            cell_metadata = cell_metadata,
                            gene_metadata = gene_annotation)
    #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
    cds <- preprocess_cds(cds, num_dim = 50)
    #umap降维
    cds <- reduce_dimension(cds, preprocess_method = "PCA")
    p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
    ggsave(file = "step1.png",plot = p1, width = 8,height = 8)
    ##从seurat导入整合过的umap坐标
    cds.embed <- cds@int_colData$reducedDims$UMAP
    int.embed <- Embeddings(myeloid_sce, reduction = "umap")
    int.embed <- int.embed[rownames(cds.embed),]
    cds@int_colData$reducedDims$UMAP <- int.embed
    p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')



}






# Figure 31 Myeloid subpopulation survival analysis
if (F) {

    library(data.table)
    library(magrittr)
    library(biomaRt)
    library(limma)
    library(survival)
    library(survminer)
    library(dplyr)
    library(Seurat)

    fig_outdir <- paste0(outdir, "/", "Figure31")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    # input parameter
    global_seurat_rds = "/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds"
    subcluster_seurat_rds = myeloid_seurat_rds
    

    # load TCGA data: LIHC_expr: LIHC TCGA expression matrix;
    # LIHC_survival: LIHC TCGA survival information, sample/OS/OS.time;

    LIHC_expr <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.expr.rds")
    LIHC_survival <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.survival.rds")

    # load global Seurat object
    global_sce <- readRDS(global_seurat_rds)

    # load subtype Seurat object
    subcluster_sce <- readRDS(subcluster_seurat_rds)

    celltype_array <- as.character(unique(subcluster_sce$clusters))
    
    for (celltype in celltype_array) {
        # get cell id in the specific subpopulations
        specific_cells <- Cells(subcluster_sce[,subcluster_sce@meta.data$clusters == celltype])
        global_sce@meta.data$yesORno <- ifelse(Cells(global_sce) %in% specific_cells, "yes", "no")

        # find markers (specific subpopulations vs global cells)
        markers <- FindMarkers(global_sce, ident.1 = "yes", 
                            group.by = 'yesORno',
                            only.pos = TRUE)
        markers$Gene.name.uniq <- rownames(markers)

        write.csv(markers,file = paste0(celltype,"_allmarkers.csv"),quote = F)

        top20 <- markers %>% top_n(n = 30, wt = avg_log2FC)

        gene_array <- unique(top20$Gene.name.uniq)



        ### input gene signature array to evaluate the survival 
        gene <- gene_array
        gene <- intersect(gene,row.names(LIHC_expr))
        data <- LIHC_expr[gene,]
        data <- t(data)
        mean_expr <- apply(data[,gene],1,mean)
        attr(mean_expr,"names") <- NULL
        data <- cbind(data,mean_expr)
        data <- data.frame(data)
        data$sample <- row.names(data)
        row.names(data) <- NULL
        surdata <- merge(LIHC_survival,data,by = c("sample"))
        surdata$level <- ifelse(surdata[,"mean_expr"] > median(surdata[,"mean_expr"]),'High','Low')
        #surdata$level <- ifelse(surdata[,"mean_expr"] > as.numeric(quantile(surdata[,"mean_expr"],0.75)),'High','Low')
        fit <- survfit(Surv(OS.time,OS)~level,data = surdata)
        print(surv_pvalue(fit)$pval)
        
        p <- ggsurvplot(fit,
            title = celltype,
            pval = T, 
            pval.method = T, 
            risk.table = T,
            size = 2, # line size
            linetype = "strata", # line type by groups
            palette = c("#E7B800", "#2E9FDF"),
            conf.int = T # add confidence interval
        )
        
        png(file = paste0(celltype,".survival.png"),width = 1024,height = 768)
        print(p)
        dev.off()

        pdf(file = paste0(celltype,".survival.pdf"),width = 12,height = 8)
        print(p)
        dev.off()
    }



}
















# Figure 36: Euclidean distance between ZPT&ZPN vs ZRTvsZRN

if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure36")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    cell_stat_frame <- as.data.frame.array(table(
        as.character(myeloid_sce$orig.ident),
        as.character(myeloid_sce$clusters)
    ))

    immune_celltype_array <- as.character(unique(myeloid_sce$clusters))

    #cell_stat_frame <- cell_stat_frame[immune_celltype_array]

    cell_stat_frame$Total_cell <- rowSums(cell_stat_frame)

    cell_stat_frame$Sample <- row.names(cell_stat_frame)
    row.names(cell_stat_frame) <- NULL

    cell_stat_frame$sample_group <- paste0(substr(cell_stat_frame$Sample,1,2),
                                        substr(cell_stat_frame$Sample,5,6))

    for(celltype in immune_celltype_array){
        cell_stat_frame[celltype] = cell_stat_frame[celltype] / cell_stat_frame["Total_cell"]
    }

    cell_stat_frame[c("Total_cell")] <- NULL

    write.table(cell_stat_frame,
        file = "cell_stat.tsv", 
        row.names = F, 
        col.names = T, 
        sep = "\t", 
        quote = F
    )
    row.names(cell_stat_frame) <- cell_stat_frame$Sample
    cell_stat_frame$sample_group <- NULL
    cell_stat_frame$Sample <- NULL

    # 比列先取log10
    #cellstat_log10_frame <- as.data.frame(apply(cell_stat_frame, 2, log10))
    cellstat_log10_frame <- as.data.frame(apply((cell_stat_frame + 0.0001), 2, log10))
    # 计算欧式距离
    dist <- dist(cellstat_log10_frame,method='euclidean') %>% as.matrix()

    dist <- dist[upper.tri(dist)]

    
    dist_frame <- data.frame(
        t(combn(row.names(cellstat_log10_frame), 2)),
        as.numeric(dist)
    )
    colnames(dist_frame) <- c("sample1","sample2","dist")

    write.table(
        x = dist_frame,
        file = "distance_euclidean.xls",
        sep = "\t",
        col.name = TRUE,
        row.name = FALSE,
        quote=FALSE
    )
    if (F) {
    consistent_sample_dist_frame <- dist_frame %>%
        filter(substr(sample1, 1, 4) == substr(sample2, 1, 4)) %>%
        mutate(group = substr(sample1,1,2))
    }

    sample_dist_frame <- dist_frame %>%
        mutate(type = case_when(
            ((substr(sample1, 1, 4) == substr(sample2, 1, 4)) & (substr(sample1, 1, 2) == "ZP")) ~ "ZP_tumor_matched",
            ((substr(sample1, 1, 4) == substr(sample2, 1, 4)) & (substr(sample1, 1, 2) == "ZR")) ~ "ZR_tumor_matched",
            ((substr(sample1, 1, 4) != substr(sample2, 1, 4)) &
                (substr(sample1, 5, 5) == "T") & (substr(sample2,5,5) == "T")) ~ "tumor_unmatched",
            ((substr(sample1, 1, 4) != substr(sample2, 1, 4)) &
                (substr(sample1, 5, 5) == "N") & (substr(sample2,5,5) == "N")) ~ "normal_unmatched",
            ((substr(sample1, 1, 4) != substr(sample2, 1, 4)) &
                ((substr(sample1, 5, 5) == "T") & (substr(sample2,5,5) == "N") |
                (substr(sample1, 5, 5) == "N") & (substr(sample2,5,5) == "T"))) ~ "tumor_normal_unmatched",
        ))

    compaired <- list(
        c("ZP_tumor_matched", "ZR_tumor_matched"),
        c("tumor_normal_unmatched","ZP_tumor_matched"),
        c("tumor_normal_unmatched","ZR_tumor_matched")
    )

    gp <- ggplot(sample_dist_frame, aes(x = type, dist, fill = type)) +
        geom_boxplot() +
        #scale_fill_manual(values = main_group_cols) +
        theme(
            axis.text.x = element_text(angle=30, hjust=1, vjust=1),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA)
        ) +
        xlab("") +
        ylab("Euclidean distances") +
        geom_signif(
            comparisons = compaired,
            step_increase = 0.3, 
            map_signif_level = F,
            #test.args = "greater",
            test = t.test)

    ggsave(filename = "PHTvsLRT_Myeloid_DistanceBoxPlot.png",plot = gp,width = 8,height = 6,units = c("in"))
    ggsave(filename = "PHTvsLRT_Myeloid_DistanceBoxPlot.pdf",plot = gp,width = 8,height = 6,units = c("in"))

}







#VIP Figure37: OR of tissue preference
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure37")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    library("sscVis")
    library("data.table")
    library("grid")
    library("cowplot")
    library("ggrepel")
    library("readr")
    library("plyr")
    library("ggpubr")
    library("ggplot2")

    cellInfo.tb <- myeloid_sce@meta.data
    cellInfo.tb = data.table(cellInfo.tb)
    meta.cluster <- cellInfo.tb$clusters

    cellInfo.tb$loc <- as.factor(cellInfo.tb$group_id)
    loc.avai.vec <- levels(cellInfo.tb[["loc"]])
    count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
    freq.dist <- sweep(count.dist, 1, rowSums(count.dist), "/")
    # 计算每个细胞类型在各个组织类型中的百分比
	freq.dist.bin <- floor(freq.dist * 100 / 10)
	print(freq.dist.bin)

    # test.dist.table function
    min.rowSum = 0
    # 过滤掉数量小于min.rowSum的细胞类型
    count.dist <- count.dist[rowSums(count.dist) >= min.rowSum, , drop = F]
    # 计算组织类型细胞数量
    sum.col <- colSums(count.dist)
    # 计算细胞类型总数
    sum.row <- rowSums(count.dist)
    # table 转 data.frame (使用as.data.frame)
    count.dist.tb <- as.data.frame(count.dist)

    # 宽数据转长数据
    setDT(count.dist.tb,keep.rownames=T)
    count.dist.melt.tb <- melt(count.dist.tb, id.vars="rn")
    colnames(count.dist.melt.tb) <- c("rid", "cid", "count")
    
    # 循环遍历data.frame
    count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
        # this.row: 细胞类型
        this.row <- count.dist.melt.tb$rid[i]
        # this.col: 组织类型
        this.col <- count.dist.melt.tb$cid[i]
        # this.c: 细胞类型x组织类型-细胞数量
        this.c <- count.dist.melt.tb$count[i]
        # 计算这个组织类型中其他细胞类型数量总和
        other.col.c <- sum.col[this.col] - this.c
        # 生成列联表
        this.m <- matrix(c(
            this.c,
            # 计算这个细胞类型中其他组织类型的数量
            sum.row[this.row] - this.c,
            other.col.c,
            sum(sum.col) - sum.row[this.row] - other.col.c
        ),
        ncol = 2
        )
        # 使用fisher.test计算p.value和OR值
        res.test <- fisher.test(this.m)
        data.frame(rid = this.row,
                cid = this.col,
                p.value = res.test$p.value,
                OR = res.test$estimate)
    }))
    count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
    count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]

    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]


    pdf.width = 3
    pdf.height = 5
    out.prefix = "myeloid"
    sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)


}







