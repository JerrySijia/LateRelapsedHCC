###################################################################
# Copyright(c) Jerry'gene, all rights reserved
# @file        : scRNAGlobalCellAnalysis.R
# @author      : Jerry
# @mail        : cuisijia@qq.com
# @revision    :
# @description :
###################################################################

# load library
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)
library(ggsignif)
# library(MySeuratWrappers)
library(RColorBrewer)
library(ggrepel)

global_seurat_rds <- "/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds"
cell_type <- "/P01_GlobalCell/P999_GlobalCell_Analysis/data/cell.type.txt"
sample_list <- "/P01_GlobalCell/P999_GlobalCell_Analysis/data/sample_list.txt"
cell_marker <- "/P01_GlobalCell/P999_GlobalCell_Analysis/data/cell_marker.txt"
celltype_group <- "/P01_GlobalCell/P999_GlobalCell_Analysis/data/celltype_group.txt"
hypoxia_siganture_file = "/signature/hypoxia_signature.txt"
outdir <- "/P01_GlobalCell/P999_GlobalCell_Analysis"

# Read Global Seurat RDS Data
cat("Reading Global Seurat RDS data...\n")
global_sce <- readRDS(global_seurat_rds)
cat("Complete.\n")



plan("multicore", workers = 16)
options(future.globals.maxSize = 150 * 1024^3)

# Global Color Settings

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


sample_colors <- c(
    PH01T = "#E5D2DD",
    PH01N = "#53A85F",
    PH02T = "#F1BB72",
    PH02N = "#F3B1A0",
    PH03T = "#D6E7A3",
    PH03N = "#57C3F3",
    PH04T = "#476D87",
    PH04N = "#E95C59",
    PH05T = "#E59CC4",
    PH05N = "#AB3282",
    PH06T = "#23452F",
    PH06N = "#BD956A",
    PH07T = "#8C549C",
    PH07N = "#585658",
    PH08T = "#9FA3A8",
    PH08N = "#E0D4CA",
    PH09T = "#5F3D69",
    PH09N = "#C5DEBA",
    PH10T = "#58A4C3",
    PH10N = "#E4C755",
    PH11T = "#F7F398",
    PH11N = "#AA9A59",
    LR01T = "#E63863",
    LR01N = "#E39A35",
    LR02T = "#C1E6F3",
    LR02N = "#6778AE",
    LR03T = "#91D0BE",
    LR03N = "#B53E2B",
    LR04T = "#712820",
    LR04N = "#DCC1DD",
    LR05T = "#CCE0F5",
    LR05N = "#CCC9E6",
    LR06T = "#625D9E",
    LR06N = "#68A180",
    LR07T = "#3A6963",
    LR07N = "#968175"
)

group_cols <- c(
    PHT = "#00af3e", 
    PHN = "#ef9020",
    LRT = "#CD1E24",
    LRN = "#0081b4" 
)

main_group_cols <- c(
    PH = "#00af3e", 
    LR = "#CD1E24"
)

celltype_colors <- c(
    "Myeloid" = "#262C68",
    "Tcell" = "#CD1E24",
    "NK" = "#1E843F",
    "Plasma/B" = "#84278B",
    "Endothelial" = "#0B6E78",
    "Malignant" = "#EF7A2A",
    "HSC" = "#B968A5",
    "Epithelial" = "#BDD240",
    "EPCAM+" = "#EF7A2A",
    "Mastcell" = "#0081b4"
)


immune_celltype_colors <- c(
    "Myeloid" = "#262C68",
    "Tcell" = "#CD1E24",
    "NK" = "#1E843F",
    "Plasma/B" = "#84278B",
    "Mastcell" = "#0081b4"
)

#VIP Figure 0: Add Celltype information to Seurat Object
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure0")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)
    
    global_sce_rdata <- "/P01_GlobalCell/P01_seurat_umap/final_seurat_result/HCC_dim_20_gene_3000_res_1_dist_2_harmony_TRUE_merge_for_tsne-umap.Rdata"

    global_sce <- eval(parse(text = load(global_sce_rdata)))

    cell_type <- "
        0	Myeloid
        1	Tcell
        2	Tcell
        3	Myeloid
        4	Tcell
        5	EPCAM+
        6	NK
        7	Myeloid
        8	NK
        9	NK
        10	Tcell
        11	Myeloid
        13	Endothelial
        14	Myeloid
        15	Myeloid
        16	Plasma/B
        18	EPCAM+
        19	Tcell
        20	Myeloid
        21	EPCAM+
        22	Mastcell
        23	Tcell
        24	Myeloid
        25	EPCAM+
        26	EPCAM+
        27	HSC
        28	EPCAM+
        29	Endothelial
        30	Myeloid
        32	Tcell"

    
    # the first column is cluster index, the second is cell type
    cell_type <- read.table(
        text = gsub("        ","",cell_type),
        sep = "\t",
        col.names = c("index","cell_type")
        stringsAsFactors = F
    )

    global_sce <- global_sce[,(global_sce@meta.data$seurat_clusters %in% cell_type$index)]
    # 取第二列转成character类型的向量
    new.cluster.ids <- as.character(cell_type$cell_type)
    names(new.cluster.ids) <- levels(global_sce)
    global_sce <- RenameIdents(global_sce, new.cluster.ids)
    global_sce@meta.data$clusters = global_sce@active.ident

    global_sce <- readRDS("/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds")

    # step1 生成细胞坐标矩阵
    data = global_sce@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        cbind(clusters = global_sce@meta.data$clusters)

    # step2 获取要添加标签的位置
    class_avg <- data %>% 
        group_by(clusters) %>% 
        summarise(
            UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2)
        )

    # step3 绘图
    umap <- ggplot(data ,aes(x=UMAP_1,y=UMAP_2))+
        geom_point(aes(color = clusters),size = 0.5, alpha = 0.2) + 
        scale_color_manual(values = celltype_colors)+
        geom_text(aes(label = clusters), color = "#002B36", data = class_avg, size = 10) +
        # geom_label(data = class_avg, 
        #    aes(x = UMAP_1, y = UMAP_2, label = clusters), 
        #    color = "#EB7D29",
        #    label.padding = unit(0.5, "lines"),
        #    fontface = "bold",
        #    size = 6) +
        xlab('UMAP dimension 1') +
        ylab('UMAP dimension 2') +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.title.x = element_text(size=18),
                axis.title.y = element_text(size=18),
                axis.text = element_text(size=18),
                legend.position = 'none') + 
        guides(colour = guide_legend(override.aes = list(size = 6)))

    ggsave(file = "GlobalUMAPByCelltype.pdf",plot = umap,width = 12,height = 12)
    ggsave(file = "GlobalUMAPByCelltype.png",plot = umap,width = 12,height = 12)

    # save
    saveRDS(sce, file = paste0('HCC.rds'), compress = FALSE)
}



#VIP Figure 1: Plot Global UMAP by group
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure1")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)

    ## add group information
    samples_list <- read.table(sample_list, header = F, sep = "\t")
    group_array <- samples_list[, "V2"]
    names(group_array) <- samples_list[, "V1"]
    global_sce$group <- group_array[global_sce$orig.ident]

    ### do not use DimPlot Function, change to ggplot2
    ### step1 get umap coordinates
    data = global_sce@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        cbind(group = global_sce@meta.data$group_id)

    ### step2 Plot Figures
    umap <- ggplot(data, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(color = group), size = 0.2, alpha = 1) +
        scale_color_manual(values = group_cols) +
        xlab("UMAP dimension 1") +
        ylab("UMAP dimension 2") +
        theme_bw() +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank(),
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text = element_text(size = 18),
            legend.text = element_text(size = 18),
            legend.key.size = unit(0.4, "inches")
        ) +
        guides(colour = guide_legend(override.aes = list(size = 6)))

    pdf("GlobalUMAPByGroup.pdf", w = 12, h = 8)
    print(umap)
    dev.off()
    png("GlobalUMAPByGroup.png", width = 1200, height = 800)
    print(umap)
    dev.off()

}



# Figure 2: Plot Global UMAP by sample
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure2")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    ## Plot Figures
    data = global_sce@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        cbind(clusters = global_sce@meta.data$orig.ident)

    umap <- ggplot(data ,aes(x=UMAP_1,y=UMAP_2))+
    geom_point(aes(color=clusters),size = 0.2, alpha = 1)+ 
    scale_color_manual(values = sample_colors)+
    xlab('UMAP dimension 1')+
    ylab('UMAP dimension 2')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            legend.title=element_blank(),
            axis.title.x =element_text(size=18), 
            axis.title.y=element_text(size=18),
            axis.text=element_text(size=18),
            legend.text=element_text(size=18),
            legend.key.size = unit(0.4, "inches")
    )+guides(colour = guide_legend(override.aes = list(size=6)))

    pdf("Figure2.pdf",w=12,h=8)
    print(umap)
    dev.off()
    png("Figure2.png",width=1200,height=800)
    print(umap)
    dev.off()
}


# Figure 3: Plot Global Marker Gene Heatmap By Cluster
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure3")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    Idents(global_sce) <- "seurat_clusters"
    cell_marker_frame = read.table(file = cell_marker,
                sep = "\t", 
                header = T,
                stringsAsFactors = F)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$geneSymbol),]

    AverageExp <- AverageExpression(global_sce,
        features = cell_marker_frame$geneSymbol, 
        group.by = "clusters"
    )

    expr <- AverageExp$RNA

    heatmapColorRamp = colorRampPalette(c("#0658D6","#FDF8EE","#ED7A03"))(64)

    expr <- ScaleData(expr)

    ha = rowAnnotation(
        labels = cell_marker_frame$cluster,
        col = list(bar = celltype_colors),
        annotation_legend_param = list(title = "", title_gp = gpar(fontsize = 20), 
                                               labels_gp = gpar(fontsize = 20)))
    expr <- expr[, unique(cell_marker_frame$cluster)]
    expr <- expr[,c("Endothelial","EPCAM+","HSC","Mastcell","Myeloid","NK","Plasma/B","Tcell")]
    ht <- Heatmap(as.matrix(expr),
            col = heatmapColorRamp,
            right_annotation = ha,
            rect_gp = gpar(col = "white", lwd = 1),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 22),
            row_split = cell_marker_frame$cluster,
            row_title_side = "right",
            row_title_rot = 0, 
            row_title_gp = gpar(fontsize = 25),
            heatmap_legend_param = list(title = "", 
                title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 30)),
            cluster_columns = F,
            cluster_rows = F)

    png("Figure3.png",width = 1080,height = 1080 * 1.5)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

    pdf("Figure3.pdf",width = 12,height = 16)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

}


#VIP Figure 4: Plot Global UMAP by marker
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure4")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    marker_frame <- read.table(file = cell_marker, sep = "\t", header = T, stringsAsFactors = F)

    marker_list <- split(marker_frame,marker_frame$cluster)
    for (marker in marker_list){
        cell <- marker[, 1][1]
        cat("Parsing ",cell,"\n")
        markers <- marker[, 2]
        
        plot_array <- c()
        
        p1 <- FeaturePlot(global_sce,
            features = markers, 
            cols = c("lightgrey", "#e32119"), 
            pt.size = 0.5,
            raster = F
        )

        ggsave(plot = p1, filename = paste0(cell, "_Figure4.pdf"), w = 12, h = 8, units = c("in"))
        ggsave(plot = p1, filename = paste0(cell, "_Figure4.png"), w = 12, h = 8, units = c("in"))
    }

}


#VIP Figure 5: Plot Cell type stat BarPlot
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure5")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    celltype_group <- "Myeloid	immune
        Tcell	immune
        NK	immune
        Plasma/B	immune
        Mastcell	immune
        Endothelial	nonimmune
        Malignant	nonimmune
        HSC	nonimmune
        Epithelial	nonimmune"
    # 读取celltype group
    celltype_group_frame <- read.table(
        text = gsub("        ","",celltype_group),
        sep = "\t",
        header = F,
    )
    # 获取immune的细胞类型
    immune_celltype_array <- (celltype_group_frame %>% filter(V2 %in% c("immune")))$V1

    # 生成细胞数量列联表(样本vs细胞类型)
    immune_cell_stat_frame <- as.data.frame.array(table(
        as.character(global_sce$sample_id),
        as.character(global_sce$clusters)
    ))

    immune_cell_stat_frame <- immune_cell_stat_frame[immune_celltype_array]
    write.table(
        x = immune_cell_stat_frame,
        file = "global_immune_stat.xls",
        sep = "\t",
        col.name = TRUE,
        row.name = TRUE,
        quote=FALSE
    )

    immune_cell_stat_frame$sample <- row.names(immune_cell_stat_frame)
    # 添加分组信息
    immune_cell_stat_frame$sample_group <- paste0(
        substr(immune_cell_stat_frame$sample, 1, 2),
        substr(immune_cell_stat_frame$sample, 5, 6)
    )
    immune_cell_stat_frame$sample <- NULL
    immune_cell_stat_frame <- immune_cell_stat_frame %>% gather(celltype,count,-sample_group)
    immune_cell_stat_frame$sample_group <- factor(immune_cell_stat_frame$sample_group,levels = names(group_cols))
    bp <- ggplot(immune_cell_stat_frame,aes(
                x = sample_group,
                y = count,
                fill = factor(celltype)
            )) +
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

    pdf("CelltypeProportionByGroupBarPlot.pdf", w = 6, h = 8)
    print(bp)
    dev.off()
    png("CelltypeProportionByGroupBarPlot.png", width = 1200, height = 1600)
    print(bp)
    dev.off()

}

#VIP Figure 6: Immune cell frequency BarPlot by group
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure6")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

   
    celltype_group <- "Myeloid	immune
        Tcell	immune
        NK	immune
        Plasma/B	immune
        Mastcell	immune
        Endothelial	nonimmune
        Malignant	nonimmune
        HSC	nonimmune
        Epithelial	nonimmune"
    # 读取celltype group
    celltype_group_frame <- read.table(
        text = gsub("        ","",celltype_group),
        sep = "\t",
        header = F,
    )
    immune_celltype_array <- (celltype_group_frame %>% filter(V2 %in% c("immune")))$V1

    cell_stat_frame <- as.data.frame.array(table(
        as.character(global_sce$sample_id),
        as.character(global_sce$clusters)
    ))

    cell_stat_frame <- cell_stat_frame[immune_celltype_array]

    cell_stat_frame$Total_cell <- rowSums(cell_stat_frame)

    cell_stat_frame$Sample <- row.names(cell_stat_frame)
    row.names(cell_stat_frame) <- NULL

    cell_stat_frame$sample_group <- paste0(substr(cell_stat_frame$Sample,1,2),
                                        substr(cell_stat_frame$Sample,5,6))

    for(celltype in immune_celltype_array){
        cell_stat_frame[celltype] = cell_stat_frame[celltype] / cell_stat_frame["Total_cell"]
    }

    cell_stat_frame[c("Total_cell")] <- NULL

    melt_cell_stat_frame = melt(cell_stat_frame, id = c("Sample", "sample_group"))

    melt_cell_stat_frame$sample_group <- factor(melt_cell_stat_frame$sample_group,levels = names(group_cols))

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
                #test.args = list(var.equal = TRUE,alternative = "greater"), # one side or two side?
                test = wilcox.test) +
    facet_grid(~variable)

    ggsave("Figure6.png",p1,width = 8,height = 4.5,units = "in")
    ggsave("Figure6.pdf",p1,width = 8,height = 4.5,units = "in")

}



# Figure 7: Celltype hypoxia signature DotPlot for each celltype
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure7")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    hypoxia_signature_frame <- read.table(file = hypoxia_siganture_file, sep = "\t",header = F,stringsAsFactors = F)
    global_sce <- AddModuleScore(object = global_sce, list(hypoxia_signature_frame$V2),name = "hypoxia_score")

    plot <- VlnPlot(global_sce,
        features = "hypoxia_score1",
        group.by = "clusters",
        col = celltype_colors
    )

    ggsave("Figure7-VlnPlot.png",plot = plot,width = 8,height = 4.5,units = "in")
    ggsave("Figure7-VlnPlot.pdf",plot = plot,width = 8,height = 4.5,units = "in")


    p <- DotPlot(
        object = global_sce,
        features = "hypoxia_score1",
        dot.scale = 12,
        cols = c("lightgrey", "red"),
        group.by = "clusters",
    ) + theme(axis.text.x = element_text(angle = 90)) +
        theme_bw()

    ggsave("Figure7-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure7-DotPlot.pdf",p,width = 8,height = 8,units = "in")

}




# Figure 8: Plot Global Marker DotPlot CellType
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure8")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    Idents(global_sce) <- "seurat_clusters"
    cell_marker_frame = read.table(file = cell_marker,
                sep = "\t", 
                header = T,
                stringsAsFactors = F)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$geneSymbol),]

    dp <- DotPlot(global_sce,
        features = cell_marker_frame$geneSymbol,
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
        filename = "Figure8.png",
        plot = dp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "Figure8.pdf", 
       plot = dp, 
       width = 8, 
       height = 12, 
       units = c("in"))

}



# Figure 9 Global Stack VlnPlot By CellType
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure9")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    Idents(global_sce) <- "seurat_clusters"
    cell_marker_frame = read.table(file = cell_marker,
                sep = "\t", 
                header = T,
                stringsAsFactors = F)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$geneSymbol),]

    vp <- VlnPlot(global_sce,
        features = cell_marker_frame$geneSymbol,
        stack = T,
        pt.size = 0,
        group.by = "clusters",
        cols = celltype_colors,
        #direction = "horizontal",
        x.lab = "",
        y.lab = ""
    ) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    ggsave(
        filename = "Figure9.png",
        plot = vp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "Figure9.pdf", 
       plot = vp, 
       width = 8, 
       height = 12, 
       units = c("in"))

}


# Figure10: Find All Markers in Global Cell types

if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure10")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    marker_frame <- FindAllMarkers(global_sce,group.by = "clusters", only.pos = TRUE)
    AllMakers <- 'all_markers.csv'
    marker_frame <- marker_frame %>% group_by(cluster)
    write.csv(marker_frame, file=AllMakers, quote=F)

    marker_frame <- marker_frame %>%
        group_by(cluster) %>%
        top_n(3, avg_log2FC)
    marker_frame <- marker_frame[!duplicated(marker_frame$gene),]
    #DoHeatmap(global_sce, features = marker_frame$gene, group.by = "clusters", label = TRUE)
    dp <- DotPlot(global_sce,
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
        filename = "Figure10.png",
        plot = dp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "Figure10.pdf", 
       plot = dp, 
       width = 8, 
       height = 12, 
       units = c("in"))

}



# Figure 11: CellPhoneDB interactive bubble plot
if (T) {
    prefix <- "LRT"
    fig_outdir <- paste0(outdir, "/", "Figure11","/",prefix)
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    indir <- "/P01_GlobalCell/P999_GlobalCell_Analysis/Figure18/out"



    means_file <- file.path(indir, "means.txt")
    pval_file <- file.path(indir, "pvalues.txt")
    sig_means_file <- file.path(indir, "significant_means.txt")

    means_frame <- read.delim(
        file = means_file,
        sep = "\t", 
        check.names = F
    )
    pvals_frame <- read.delim(
        file = pval_file,
        sep = "\t", 
        check.names = F
    )
    sig_means_frame <- read.delim(
        file = sig_means_file,
        sep = "\t", 
        check.names = F
    )

    selected_celltypes_file <- "/P01_GlobalCell/P999_GlobalCell_Analysis/data/cellphonedb/columns.txt"
    selected_interactor_file <- "/P01_GlobalCell/P999_GlobalCell_Analysis/data/cellphonedb/rows.txt"

    selected_celltypes_frame <- read.table(file = selected_celltypes_file,sep = "\t",header = F,stringsAsFactors = F)
    selected_interactor_frame <- read.table(file = selected_interactor_file,sep = "\t",header = F,stringsAsFactors = F)

    if (F) {
    order_sequence <- function(df) {
        da <- data.frame()
        for (i in 1:length(df$gene_a)) {
            sub_data <- df[i, ]
            if (sub_data$receptor_b == "False") {
                if (sub_data$receptor_a == "True") {
                    old_names <- colnames(sub_data)
                    my_list <- strsplit(old_names[-c(1:11)], split = "\\|")
                    my_character <- paste(sapply(my_list, "[[", 2L), sapply(my_list, "[[", 1L), sep = "|")
                    new_names <- c(names(sub_data)[1:4], "gene_b", "gene_a", "secreted", "receptor_b", "receptor_a", "annotation_strategy", "is_integrin", my_character)
                    sub_data = dplyr::select(sub_data, new_names)
                    # print('Change sequence!!!')
                    names(sub_data) <- old_names
                    da = rbind(da, sub_data)
                }
            } else {
                da = rbind(da, sub_data)
            }
        }
        return(da)
    }
    }
    # ligand-receptor, one is ligand and the other is receptor
    df <- subset(
        means_frame,
        receptor_a == "True" & receptor_b == "False" | receptor_a == "False" & receptor_b == "True"
    )

    df <- df %>%
        dplyr::mutate(na_count = rowSums(is.na(df) | df == "")) %>%
        subset(na_count == 0) %>%
        dplyr::select(-na_count)

    if (F) {
        means_order <- order_sequence(df) %>%
            tidyr::unite(Pairs, gene_a, gene_b)
        pvals_order <- order_sequence(pvals_frame) %>%
            tidyr::unite(Pairs, gene_a, gene_b)
    }

    means_order <- df %>% tidyr::unite(Pairs, gene_a, gene_b)
    pvals_order <- pvals_frame %>% tidyr::unite(Pairs, gene_a, gene_b)

    means_sub <- means_order[, c('Pairs', colnames(means_frame)[-c(1:11)])]
    pvals_sub <- pvals_order[, c('Pairs', colnames(means_frame)[-c(1:11)])]
    means_gather <- tidyr::gather(means_sub, celltype, mean_expression, names(means_sub)[-1])
    pvals_gather <- tidyr::gather(pvals_sub, celltype, pval, names(pvals_sub)[-1])
    mean_pval <- dplyr::left_join(means_gather, pvals_gather, by = c('Pairs', 'celltype'))
    #create_dt(mean_pval)
    a <- mean_pval %>% dplyr::select(Pairs, celltype, pval) %>% tidyr::spread(key=celltype, value=pval)
    sig_pairs <- a[which(rowSums(a<=0.05)!=0), ]
    dim(sig_pairs)
    ## [1] 72 17

    # 保存显著性表达的受体配体对
    mean_pval_sub <- subset(mean_pval, Pairs %in% sig_pairs$Pairs)
    library(ggpubr)
    myPalette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")), space="Lab") ##color scheme

    mean_pval_sub <- mean_pval_sub %>%
        filter(celltype %in% selected_celltypes_frame$V1) %>%
        filter(Pairs %in% selected_interactor_frame$V1) %>%
        dplyr::arrange(Pairs)

    gp <- ggplot(mean_pval_sub,aes(x=Pairs, y=celltype)) +
        geom_tile(fill='white',color='black') +
        geom_point(aes(x=Pairs, 
                    y=celltype,
                    fill = mean_expression,
                    size = -log10(pval + 0.01)),
                    shape = 21,
                    color = "black",
                    stroke = 0.1) +
        scale_fill_gradientn(colors = c("yellow","#B00E26")) +
            theme(
                axis.line = element_blank(),
                axis.text.x = element_text(angle = 60, hjust = 1),
                panel.background = element_blank(),
                panel.grid.minor = element_blank()
            ) +
            xlab("") + ylab("")



    ggsave(file = "ligands_receptors_dotplot.png",gp,height = 12,width = 8,units = c("in"))
    ggsave(file = "ligands_receptors_dotplot.pdf",gp,height = 12,width = 8,units = c("in"))



}




#VIP Figure 12: Euclidean distance between PHT&PHN vs LRTvsLRN

if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure12")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    celltype_group = "
        Myeloid	immune
        Tcell	immune
        NK	immune
        Plasma/B	immune
        Endothelial	nonimmune
        Malignant	nonimmune
        HSC	nonimmune
        Epithelial	nonimmune"

    celltype_group_frame <- read.table(
        text = gsub("        ","",celltype_group),
        sep = "\t", 
        col.names = c("V1", "V2")
    )

    immune_celltype_array <- (celltype_group_frame %>% filter(V2 %in% c("immune")))$V1

    cell_stat_frame <- as.data.frame.array(table(
        as.character(global_sce$sample_id),
        as.character(global_sce$clusters)
    ))

    #cell_stat_frame <- cell_stat_frame[names(sample_colors),]

    cell_stat_frame <- cell_stat_frame[immune_celltype_array]

    cell_stat_frame$Total_cell <- rowSums(cell_stat_frame)

    cell_stat_frame$Sample <- row.names(cell_stat_frame)
    row.names(cell_stat_frame) <- NULL

    cell_stat_frame$patient_id <- paste0(substr(cell_stat_frame$Sample,1,4))

    for(celltype in immune_celltype_array){
        cell_stat_frame[celltype] = log10(cell_stat_frame[celltype] / cell_stat_frame["Total_cell"] + 0.001)
    }

    cell_stat_frame[c("Total_cell")] <- NULL

    all_dist_frame <- data.frame()
    combine_sample_matrix <- combn(cell_stat_frame$Sample,2)
    for(col_index in 1:ncol(combine_sample_matrix)){
        sample_array <- combine_sample_matrix[,col_index]
        sample_1 = sample_array[1]
        sample_2 = sample_array[2]
        distance <- dist(cell_stat_frame %>% filter(Sample %in% sample_array) %>% select(immune_celltype_array),method='euclidean')
        if ((substr(sample_1,1,2) == substr(sample_2,1,2)) & 
            (substr(sample_1,1,2) == "PH") & 
            (substr(sample_1,5,5) == substr(sample_2,5,5)) &
            (substr(sample_1,5,5) == "T")) {
            group_id <- "PHT-T"
        } else if ((substr(sample_1,1,2) == substr(sample_2,1,2)) & 
            (substr(sample_1,1,2) == "PH") & 
            (substr(sample_1,5,5) == substr(sample_2,5,5)) &
            (substr(sample_1,5,5) == "N")) {
            group_id <- "PHN-N"
        } else if ((substr(sample_1,1,2) == substr(sample_2,1,2)) & 
            (substr(sample_1,1,2) == "PH") & 
            (substr(sample_1,5,5) != substr(sample_2,5,5))) {
            group_id <- "PHT-N"
        } else if ((substr(sample_1,1,2) == substr(sample_2,1,2)) & 
            (substr(sample_1,1,2) == "LR") & 
            (substr(sample_1,5,5) == substr(sample_2,5,5)) &
            (substr(sample_1,5,5) == "T")) {
            group_id <- "LRT-T"
        } else if ((substr(sample_1,1,2) == substr(sample_2,1,2)) & 
            (substr(sample_1,1,2) == "LR") & 
            (substr(sample_1,5,5) == substr(sample_2,5,5)) &
            (substr(sample_1,5,5) == "N")) {
            group_id <- "LRN-N"
        } else if ((substr(sample_1,1,2) == substr(sample_2,1,2)) & 
            (substr(sample_1,1,2) == "LR") & 
            (substr(sample_1,5,5) != substr(sample_2,5,5))) {
            group_id <- "LRT-N"
        }
        patient_dist_frame <- data.frame(distance = distance[1],group_id = group_id)
        all_dist_frame <- rbind(all_dist_frame,patient_dist_frame)
    }


    if (T) {
        gp <- ggplot(all_dist_frame, aes(x = group_id, y = distance, fill = group_id)) +
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
                comparisons = list(c("PHT-T","LRT-T"),c("PHN-N","LRN-N"),c("PHT-N","LRT-N")),
                step_increase = 0.2, 
                map_signif_level = F,
                #test.args = "greater",
                test = t.test)

        ggsave(filename = "FigureTNDistance.png",plot = gp,width = 5,height = 6,units = c("in"))
        ggsave(filename = "FigureTNDistance.pdf",plot = gp,width = 5,height = 6,units = c("in"))
    }





    for (p in unique(cell_stat_frame$patient_id)) {
        distance <- dist(cell_stat_frame %>% filter(patient_id %in% c(p)) %>% select(immune_celltype_array),method='euclidean')
        patient_dist_frame <- data.frame(patient_id = p,distance = distance[1])
        all_dist_frame <- rbind(all_dist_frame,patient_dist_frame)
    }

    all_dist_frame <- all_dist_frame %>% mutate(group_id = substr(patient_id,1,2))

    if (T) {
        consistent_sample_dist_frame <- dist_frame %>%
            filter(substr(sample1, 1, 4) == substr(sample2, 1, 4)) %>%
            mutate(group = substr(sample1,1,2))
        # r$> head(consistent_sample_dist_frame)
        #   sample1 sample2      dist group
        # 1   LR01N   LR01T 0.2858201    LR
        # 2   LR02N   LR02T 0.4568531    LR
        # 3   LR03N   LR03T 0.3662490    LR
        # 4   LR04N   LR04T 2.4794840    LR
        # 5   LR05N   LR05T 0.9822114    LR
        # 6   LR06N   LR06T 1.5889930    LR
        gp <- ggplot(all_dist_frame, aes(x = group_id, y = distance, fill = group_id)) +
            geom_boxplot() +
            scale_fill_manual(values = main_group_cols) +
            theme(
                axis.text.x = element_text(angle=30, hjust=1, vjust=1),
                panel.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_rect(color = "black", fill = NA)
            ) +
            xlab("") +
            ylab("Euclidean distances") +
            geom_signif(
                comparisons = list(c("PH","LR")),
                step_increase = 0.3, 
                map_signif_level = F,
                #test.args = "greater",
                test = t.test)

        ggsave(filename = "Figure12.png",plot = gp,width = 5,height = 6,units = c("in"))
        ggsave(filename = "Figure12.pdf",plot = gp,width = 5,height = 6,units = c("in"))
    }




    sample_dist_frame <- dist_frame %>%
        mutate(type = case_when(
            ((substr(sample1, 1, 4) == substr(sample2, 1, 4)) & (substr(sample1, 1, 2) == "PH")) ~ "PH_tumor_matched",
            ((substr(sample1, 1, 4) == substr(sample2, 1, 4)) & (substr(sample1, 1, 2) == "LR")) ~ "LR_tumor_matched",
            ((substr(sample1, 1, 4) != substr(sample2, 1, 4)) &
                (substr(sample1, 5, 5) == "T") & (substr(sample2,5,5) == "T")) ~ "tumor_unmatched",
            ((substr(sample1, 1, 4) != substr(sample2, 1, 4)) &
                (substr(sample1, 5, 5) == "N") & (substr(sample2,5,5) == "N")) ~ "normal_unmatched",
            ((substr(sample1, 1, 4) != substr(sample2, 1, 4)) &
                ((substr(sample1, 5, 5) == "T") & (substr(sample2,5,5) == "N") |
                (substr(sample1, 5, 5) == "N") & (substr(sample2,5,5) == "T"))) ~ "tumor_normal_unmatched",
        ))

    compaired <- list(
        c("PH_tumor_matched", "LR_tumor_matched"),
        c("tumor_normal_unmatched","PH_tumor_matched"),
        c("tumor_normal_unmatched","LR_tumor_matched")
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

    ggsave(filename = "Figure12.png",plot = gp,width = 8,height = 6,units = c("in"))
    ggsave(filename = "Figure12.pdf",plot = gp,width = 8,height = 6,units = c("in"))

}





# Figure 14: CellPhoneDB Heatmap of Ligands and Receptor

if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure14")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)












}






# Figure 15: Plot Global Marker Gene Super Modified Heatmap By Cluster
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure15")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    Idents(global_sce) <- "seurat_clusters"
    cell_marker_frame = read.table(file = cell_marker,
                sep = "\t", 
                header = T,
                stringsAsFactors = F)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$geneSymbol),]

    AverageExp <- AverageExpression(global_sce,
        features = cell_marker_frame$geneSymbol, 
        group.by = "clusters"
    )

    expr <- AverageExp$RNA

    heatmapColorRamp = colorRampPalette(c("#0658D6","#FDF8EE","#ED7A03"))(64)

    expr <- ScaleData(expr)

    ha = rowAnnotation(
        labels = cell_marker_frame$cluster,
        col = list(bar = celltype_colors),
        annotation_legend_param = list(title = "", title_gp = gpar(fontsize = 20), 
                                               labels_gp = gpar(fontsize = 20)))
    expr <- expr[, unique(cell_marker_frame$cluster)]
    expr <- expr[,c("Endothelial","EPCAM+","HSC","Mastcell","Myeloid","NK","Plasma/B","Tcell")]
    ht <- Heatmap(as.matrix(expr),
            col = heatmapColorRamp,
            right_annotation = ha,
            rect_gp = gpar(col = "white", lwd = 1),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 22),
            row_split = cell_marker_frame$cluster,
            row_title_side = "right",
            row_title_rot = 0, 
            row_title_gp = gpar(fontsize = 25),
            heatmap_legend_param = list(title = "", 
                title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 30)),
            cluster_columns = F,
            cluster_rows = F)

    png("Figure3.png",width = 1080,height = 1080 * 1.5)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

    pdf("Figure3.pdf",width = 12,height = 16)
    draw(ht,padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

}






# Figure 17: Plot Cell type stat BarPlot
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure17")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    celltype_group <- "Myeloid	immune
        Tcell	immune
        NK	immune
        Plasma/B	immune
        Mastcell	immune
        Endothelial	nonimmune
        Malignant	nonimmune
        HSC	nonimmune
        Epithelial	nonimmune"
    # 读取celltype group
    celltype_group_frame <- read.table(
        text = gsub("        ","",celltype_group),
        sep = "\t",
        header = F,
    )
    # 获取immune的细胞类型
    immune_celltype_array <- (celltype_group_frame %>% filter(V2 %in% c("immune")))$V1

    # 生成细胞数量列联表(样本vs细胞类型)
    immune_cell_stat_frame <- as.data.frame.array(table(
        as.character(global_sce$group),
        as.character(global_sce$clusters)
    ))

    immune_cell_stat_frame <- immune_cell_stat_frame[immune_celltype_array]
    write.table(
        x = immune_cell_stat_frame,
        file = "global_immune_stat.xls",
        sep = "\t",
        col.name = TRUE,
        row.name = TRUE,
        quote=FALSE
    )
    immune_cell_stat_frame <- immune_cell_stat_frame %>% tibble::rownames_to_column("sample_group")

    for (group in names(group_cols)) {
        sub_immune_cell_stat_frame <- immune_cell_stat_frame %>%
            gather(celltype, count, -sample_group) %>%
            filter(sample_group %in% c(group))
        sub_immune_cell_stat_frame <- sub_immune_cell_stat_frame %>%
            mutate(proportion = count / sum(sub_immune_cell_stat_frame$count))
        gp <- ggplot(sub_immune_cell_stat_frame, aes(x = "", y = count, fill = celltype)) +
            geom_bar(stat = "identity", width = 1, color = "white") +
            coord_polar("y", start = 0) +
            scale_fill_manual(values = immune_celltype_colors) +
            geom_text(aes(label = paste0(round(proportion,2) * 100, "%")), 
                position = position_stack(vjust=0.5)) +
            theme(
                panel.background = element_blank(),
                panel.grid = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank()
            ) + ylab("") + xlab("")
        ggsave(paste0(group,"_pie.pdf"),plot = gp, width = 8, height = 8)
        ggsave(paste0(group,"_pie.png"),plot = gp, width = 8, height = 8)
    }


}











# Figure 22 CellphoneDB analysis in PHT and LRT
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure22")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    global_sce <- readRDS("/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds")
    tnk_sce <- readRDS("/P02_TNKCell/P999_TNK_Analysis/Figure0/TNK.rds")
    myeloid_sce <- readRDS("/P03_MyeloidCell/P03_Myeloid_UMAPByCellType/myeloid.rds")
    tnk_celltype_array <- tnk_sce$clusters
    # r$> head(tnk_celltype_array)
    # ZP01N_AAACCTGAGGACGAAA-1 ZP01N_AAACCTGAGTCTTGCA-1 ZP01N_AAACCTGAGTGACTCT-1 ZP01N_AAACCTGCATCTACGA-1 ZP01N_AAACCTGGTCGAATCT-1 ZP01N_AAACCTGGTCGGCACT-1 
    #                 CD8-C1                     MAIT                    NK-C2                     MAIT                   CD8-C6                   CD8-C6 
    # Levels: CD4-C1 CD8-C1 NK-C1 MAIT CD8-C6 NK-C2 CD8-C2 CD4-C2 CD8-C3 CD8-C4 T-Cycling CD8-C5 CD4-C3 NK-C3
    myeloid_celltype_array <- myeloid_sce$clusters
    subcluster_celltype_frame <- as.data.frame(c(tnk_celltype_array, myeloid_celltype_array))
    colnames(subcluster_celltype_frame) <- "subclusters"
    # r$> head(subcluster_celltype_frame)
    #                         subclusters
    # ZP01N_AAACCTGAGGACGAAA-1      CD8-C1
    # ZP01N_AAACCTGAGTCTTGCA-1        MAIT
    # ZP01N_AAACCTGAGTGACTCT-1       NK-C2
    # ZP01N_AAACCTGCATCTACGA-1        MAIT
    # ZP01N_AAACCTGGTCGAATCT-1      CD8-C6
    # ZP01N_AAACCTGGTCGGCACT-1      CD8-C6
    global_sce <- AddMetaData(global_sce, metadata = subcluster_celltype_frame)

    library(ktplots)
    pvals <- read.delim("out/pvalues.txt", check.names = FALSE)
    means <- read.delim("out/means.txt", check.names = FALSE)
    decons <- read.delim("out/deconvoluted.txt", check.names = FALSE)

    interaction_annotation_text <- "
        T cell Inhibition	PDCD1_PDCD1LG2
        T cell Inhibition	CTLA4_CD80
        T cell Inhibition	CTLA4_CD86
        T cell Inhibition	PVR_TIGIT
        T cell Inhibition	TIGIT_NECTIN2
        T cell Inhibition	LGALS9_HAVCR2
        T cell Inhibition	BTLA_TNFRSF14
        T cell Inhibition	SPP1_CD44
        Chemokines	CXCL12_CXCR3
        Chemokines	CXCR6_CXCL16
        Chemokines	CXCL12_CXCR4
        M2-polarization	CSF1R_CSF1
        M2-polarization	CD74_MIF
        M2-polarization	SIRPA_CD74
        Cell Adhesion	CRTAM_CADM1"

    interaction_annotation <- read.table(
        text = gsub("        ","",interaction_annotation_text),
        sep = "\t", 
        col.names = c("role","interaction")
    )
	# control the LR show in the dot plot
	means <- means %>% filter(interacting_pair %in% interaction_annotation$interaction)
	pvals <- pvals %>% filter(interacting_pair %in% interaction_annotation$interaction)

    gp <- plot_cpdb(cell_type1 = 'CD8-C4', cell_type2 = 'Macrophage', scdata = global_sce,
        idents = 'subclusters', # column name where the cell ids are located in the metadata
        means = means, pvals = pvals,
        col_option = "maroon",
        highlight = "#480C5C",
		keep_significant_only = T,
        genes = gene_frame$gene) +
        small_axis(fontsize = 16) + 
		small_grid() + 
		small_guide(keysize=.8) + 
		small_legend(keysize=.8,fontsize = 20) # some helper functions included in ktplots to help with the plotting

    ggsave(file = "Macrophage_cpd_dotplot.png",plot = gp, w = 20 ,h = 13)
    ggsave(file = "Macrophage_cpd_dotplot.pdf",plot = gp, w = 20 ,h = 13)


    gp <- plot_cpdb(cell_type1 = 'CD8-C4', cell_type2 = 'Monocyte', scdata = global_sce,
        idents = 'subclusters', # column name where the cell ids are located in the metadata
        means = means, pvals = pvals,
        col_option = "maroon",
        highlight = "#480C5C",
		keep_significant_only = T,
        genes = gene_frame$gene) +
        small_axis(fontsize = 14) + 
		small_grid() + 
		small_guide(keysize=.8) + 
		small_legend(keysize=.8,fontsize = 20) # some helper functions included in ktplots to help with the plotting

    ggsave(file = "Monocyte_cpd_dotplot.png",plot = gp, w = 20 ,h = 13)
    ggsave(file = "Monocyte_cpd_dotplot.pdf",plot = gp, w = 20 ,h = 13)

    # https://github.com/zktuong/ktplots

    celltype_colors <- c(
        "CD8-C4" = "#84278B",
        "Macrophage-C1" = "#EF7A2A",
        "Macrophage-C2" = "#B968A5",
        "Macrophage-C3" = "#BDD240",
        "Macrophage-C4" = "#E63863",
        "Macrophage-C5" = "#E4C755"
    )

    p <- plot_cpdb2(cell_type1 = 'CD8-C4', cell_type2 = 'Macrophage',
        scdata = as.SingleCellExperiment(global_sce),
        idents = 'subclusters', # column name where the cell ids are located in the metadata
        means = means,
        pvals = pvals,
        deconvoluted = decons, # new options from here on specific to plot_cpdb2
        # desiredInteractions = list(
        #     c('CD8-C4', 'Macrophage-C1'),
        #     c('Macrophage-C1', 'CD8-C4')),
        # 这里不是用colname来访问，是用index来访问
        interaction_grouping = interaction_annotation[,c(2,1)],
        edge_group_colors = c(
            "T cell Inhibition" = "#e15759",
            "Chemokines" = "#59a14f",
            "M2-polarization" = "#4e79a7",
            "Cell Adhesion" = "#9c755f"
            ),
        node_group_colors = celltype_colors,
        keep_significant_only = TRUE,
        standard_scale = TRUE,
        remove_self = TRUE
    ) + small_legend(keysize=.8,fontsize = 20)

    ggsave(file = "cpdb_circos.png", plot = p, w = 18, h = 18)
    ggsave(file = "cpdb_circos.pdf", plot = p, w = 18, h = 18)


    p <- plot_cpdb4(
            interaction = 'CSF1-CSF1R',
            cell_type1 = 'CD8-C4', cell_type2 = 'Macrophage',
            scdata = as.SingleCellExperiment(global_sce),
            idents = 'subclusters',
            means = means,
            pvals = pvals,
            deconvoluted = decons,
            keep_significant_only = TRUE,
            standard_scale = TRUE,
        )

    ggsave(file = "cpdb_circos2.png",plot = p,w = 8, h = 8)
    ggsave(file = "cpdb_circos2.pdf",plot = p,w = 8, h = 8)



}





# Figure 23 CellphoneDB analysis in both PHN and LRN
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure23")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    tnk_celltype = c(
        "CD8-C1",
        "CD8-C2",
        "CD8-C3",
        "CD8-C4",
        "CD8-C5",
        "CD8-C6",
        "CD4-C1",
        "CD4-C2",
        "CD4-C3",
        "MAIT",
        "T-Cycling",
        "NK-C1",
        "NK-C2",
        "NK-C3"
    )
    myeloid_celltype <- c(
        "Macrophage-C1",
        "Macrophage-C2",
        "Macrophage-C3",
        "Macrophage-C4",
        "Macrophage-C5",
        "Macrophage-C6",
        "Monocyte-C1",
        "Monocyte-C2",
        "Monocyte-C3",
        "Monocyte-C4",
        "Monocyte-C5",
        "DC-C1",
        "DC-C2",
        "DC-C3"
    )

    # loading global, tnk, and myeloid seurat rds data
    global_sce <- readRDS("/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds")
    tnk_sce <- readRDS("/P02_TNKCell/P999_TNK_Analysis/Figure0/TNK.rds")
    myeloid_sce <- readRDS("/P03_MyeloidCell/P03_Myeloid_UMAPByCellType/myeloid.rds")
    
    if (T) {
        sample_list <- c(
            "ZP01N",
            "ZP02N",
            "ZP03N",
            "ZP04N",
            "ZP05N",
            "ZP06N",
            "ZP07N",
            "ZP08N",
            "ZP09N",
            "ZP10N",
            "ZP11N",
            "ZR01N",
            "ZR02N",
            "ZR03N",
            "ZR04N",
            "ZR05N",
            "ZR06N",
            "ZR07N"
        )

        sub_global_sce <- subset(global_sce, orig.ident %in% sample_list)
        sub_tnk_sce <- subset(tnk_sce, orig.ident %in% sample_list)
        sub_myeloid_sce <- subset(myeloid_sce, orig.ident %in% sample_list)
    }

    # get subset of tnk
    sub_tnk_sce <- subset(sub_tnk_sce, clusters %in% tnk_celltype)
    tnk_cell_id <- Cells(sub_tnk_sce)
    tnk_celltype_frame <- data.frame(cell_id = tnk_cell_id,celltype = sub_tnk_sce$clusters)
    
    sub_myeloid_sce <- subset(sub_myeloid_sce, clusters %in% myeloid_celltype)
    myeloid_cell_id <- Cells(sub_myeloid_sce)
    myeloid_celltype_frame <- data.frame(cell_id = myeloid_cell_id,celltype = sub_myeloid_sce$clusters)
    # never normalize data!!!
    # sce <- NormalizeData(sce, normalization.method =  "RC", scale.factor = 100000)
    # sce <- ScaleData(sce)
    selected_cell_id <- c(tnk_cell_id,myeloid_cell_id)
    sub_global_sce <- subset(sub_global_sce, cells = selected_cell_id)

    celltype_frame <- rbind(tnk_celltype_frame,myeloid_celltype_frame)

    #cell.names <- sapply(strsplit(cell.names$V1,"-",fixed=T), function(x) paste0(x[1], "-1"))
    data <- data.frame(GetAssayData(sub_global_sce, 'data'), check.names = FALSE)
    data$Gene <- row.names(data)
    
    #ensembl_frame <- read.table(file = ensembl_gene_file,sep = "\t",header = T,stringsAsFactors = F)

    #data <- merge(data,ensembl_frame,by = c("Gene"))

    #data$Gene <- NULL
    #colnames(data)[ colnames(data) == "EnsGeneName" ] = 'Gene'
    data <- data %>% dplyr::select("Gene",everything())

    #data <- data.frame('cellname'=cell.name)
    #data$sample <- sce@meta.data$orig.ident
    #data$cellCluster <- sce@meta.data$clusters
    cat("Generating count data.\n")
    library(data.table)
    fwrite(
        x = data,
        file = "count.txt",
        sep = "\t",
        col.name = TRUE,
        row.name = FALSE,
        quote=FALSE
    )

    fwrite(
        x = celltype_frame,
        file = "celltype.txt",
        sep = "\t",
        col.name = TRUE,
        row.name = FALSE,
        quote = FALSE
    )

    system("/mnt/share01/tools/miniconda/envs/cpdb/bin/cellphonedb method statistical_analysis --counts-data gene_name --threads 32 celltype.txt count.txt")
    system("/mnt/share01/tools/miniconda/envs/cpdb/bin/cellphonedb plot heatmap_plot celltype.txt")
    system("/mnt/share01/tools/miniconda/envs/cpdb/bin/cellphonedb plot dot_plot")
    system("/mnt/share01/tools/miniconda/envs/cpdb/bin/cellphonedb plot dot_plot --rows in/rows.txt --columns in/columns.txt")

}




#VIP Figure 24 DotPlot of Global Cells
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure24")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    global_markers <- "
        Endothelial	VWF
        Endothelial	TM4SF1
        Endothelial	INSR
        EPCAM+	EPCAM
        EPCAM+	APOA2
        EPCAM+	ALB
        EPCAM+	APOA1
        EPCAM+	AMBP
        HSC	RGS5
        HSC	COL1A1
        HSC	ACTA2
        HSC	PDGFRB
        Mastcell	TPSAB1
        Mastcell	MS4A2
        Mastcell	CPA3
        Myeloid	CD68
        Myeloid	CD163
        Myeloid	CD14
        NK	GNLY
        NK	GZMB
        NK	KLRD1
        NK	KLRF1
        Plasma/B	MZB1
        Plasma/B	SDC1
        Plasma/B	CD79A
        Tcell	CD3D
        Tcell	CD3E
        Tcell	CD8A
        Tcell	IL7R"

    global_markers_frame <- read.table(
		text = gsub("        ","",global_markers),
		sep = "\t", 
		col.names = c("CellType", "Gene")
	)

    celltype_order <- c("HSC","Mastcell","Plasma/B","Endothelial","NK","EPCAM+","Tcell","Myeloid")

    g0_plot = DotPlot(global_sce, features=split(global_markers_frame$Gene, global_markers_frame$CellType)[rev(celltype_order)],
                cols= c("lightyellow", "red3"))+
        geom_point(aes(size = avg.exp), shape = 21, colour="black", stroke=0.5) +
        RotatedAxis() +
        scale_color_gradientn(colours = viridis::viridis(20), 
            limits = c(0,3), 
            oob = scales::squish) +
        theme(
            # 各个画板
            panel.border = element_rect(color="black"),#要边框
            panel.spacing = unit(1, "mm"), #画板间距
            axis.title = element_blank(), #去掉 坐标轴 label
        )
    ggsave(file = "DotPlotByCellType.pdf", plot = g0_plot, width = 16, height = 4)


    # 使用 ggplot2 画左边的文字和彩色圆圈。
    df1 = data.frame(x = 0, y= celltype_order, stringsAsFactors = F)
    df1$y = factor(df1$y, levels = rev(celltype_order))
    # df1$y
    g1_left = ggplot(df1, aes(x,y, color=y))+
        geom_point(size=6, show.legend = F)+
        scale_color_manual(values = celltype_colors[celltype_order])+
        theme_classic()+
        scale_x_continuous(expand=c(0,0))+
        theme(
            plot.margin =margin(r=0), #no margin on the right
            axis.title = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.text.y=element_text(size=12)
        )
    ggsave(file = "left.pdf", plot = g1_left, width = 4, height = 4)



    # 拼合图形
    library(cowplot)
    # https://wilkelab.org/cowplot/articles/aligning_plots.html
    # we can align both the bottom and the top axis (axis = "bt").
    g <- plot_grid(g1_left, g0_plot, align ="h", axis="bt", rel_widths = c(1, 9))
    ggsave(file = "DotPlotByCellTypeCombine.pdf", plot = g, width = 16, height = 4)


}




# Figure 27 survival analysis in TCGA dataset
if (F) {

    library(data.table)
    library(magrittr)
    library(biomaRt)
    library(limma)
    library(survival)
    library(survminer)
    library(dplyr)
    library(Seurat)
    library(GSVA)

    fig_outdir <- file.path(outdir, "Figure51")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    # load TCGA data: LIHC_expr: LIHC TCGA expression matrix;
    # LIHC_survival: LIHC TCGA survival information, sample/OS/OS.time;

    LIHC_expr <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.expr.rds")
    # r$> head(LIHC_expr[1:5,1:5])
    #     TCGA-DD-A4NG-01A TCGA-G3-AAV4-01A TCGA-2Y-A9H1-01A TCGA-BC-A10Y-01A TCGA-K7-AAU7-01A
    # A4GNT         0.000000         0.000000        0.0000000        0.0000000        0.2813901
    # AAAS          2.738210         3.981126        2.3862749        3.7305095        3.6792481
    # AACS          0.589942         0.854231        0.1818565        0.6532292        1.7143228
    # AADAC         7.376066        10.259319        4.5717290        6.4303159        5.7398575
    # AADAT         1.255301         1.532797        0.5527291        3.9032253        1.6114995

    LIHC_survival <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.survival.rds")
    #     r$> head(LIHC_survival)
    #             sample OS     OS.time
    # 1 TCGA-FV-A495-11A  0 0.002739726
    # 2 TCGA-FV-A495-01A  0 0.002739726
    # 3 TCGA-ED-A7PZ-01A  0 0.016438356
    # 4 TCGA-ED-A97K-01A  0 0.016438356
    # 5 TCGA-ED-A7PX-01A  0 0.016438356
    # 6 TCGA-BC-A3KF-01A  0 0.021917808

    cpdb_signature <- "
        TAM	VEGFA
        TAM	CCL2
        TAM	CD81
        TAM	HLA-DRA
        M2	CD163
        M2	MSR1
        M2	MRC1
        Exhausted	CD8A
        Exhausted	CD8B
        Exhausted	PDCD1
        Exhausted	CTLA4
        Exhausted	TIGIT
        Exhausted	BTLA
        LR	CD86
        LR	HAVCR2
        LR	LGALS9
        LR	PVR
        LR	SPP1
        LR	SIRPA
        LR	CD47
        LR	CD44
        LR	CRTAM
        LR	TNFRSF14
        LR	ICAM1
        LR	VCAM1
        LR	PDCD1LG2
        LR	CXCL12
        LR	CXCR6
        LR	CXCR3
        LR	CXCL16
        LR	CXCR4
        LR	NECTIN2
        LR	MIF
        LR	CD74
        LR	CSF1
        LR	CSF1R"
    cpdb_signature_frame <- read.table(
        text = gsub("        ","",cpdb_signature),
        sep = "\t", 
        col.names = c("cluster", "marker")
    )
    
    # gsva
    gs_exp <- gsva(as.matrix(LIHC_expr), 
        list("cpdb" = cpdb_signature_frame$marker), 
        kcdf = "Poisson", 
        min.sz = 10)


    data <- as.data.frame(t(gs_exp))
    data$sample <- row.names(data)
    row.names(data) <- NULL

    data <- as.data.frame(t(gs_exp))

    # expression
    gene <- cpdb_signature_frame$marker
    gene <- intersect(gene,row.names(LIHC_expr))
    data <- LIHC_expr[gene,]
    data <- t(data)
    cpdb <- apply(data[,gene],1,mean)
    attr(cpdb,"names") <- NULL
    data <- cbind(data,cpdb)
    data <- data.frame(data)
    data$sample <- row.names(data)
    row.names(data) <- NULL
    


    surdata <- merge(LIHC_survival,data,by = c("sample"))
    surdata$level <- ifelse(surdata[,"cpdb"] > median(surdata[,"cpdb"]),'High','Low')
    #surdata$level <- ifelse(surdata[,"mean_expr"] > as.numeric(quantile(surdata[,"mean_expr"],0.75)),'High','Low')
    fit <- survfit(Surv(OS.time,OS)~level,data = surdata)
    print(surv_pvalue(fit)$pval)

    p <- ggsurvplot(fit,
        title = "Exhausted CD8/TAM interaction signature",
        pval = T, 
        pval.method = T, 
        risk.table = T,
        size = 2, # line size
        linetype = "strata", # line type by groups
        palette = c("#E7B800", "#2E9FDF"),
        conf.int = T # add confidence interval
    )
    
    png(file = "cpdb.survival.png",width = 1024,height = 768)
    print(p)
    dev.off()

    pdf(file = "cpdb.survival.pdf",width = 12,height = 8)
    print(p)
    dev.off()


}











# Figure 28 survival analysis
if (F) {

    library(data.table)
    library(magrittr)
    library(biomaRt)
    library(limma)
    library(survival)
    library(survminer)
    library(dplyr)
    library(Seurat)
    library(GSVA)
    library(GEOquery)

    fig_outdir <- file.path(outdir, "Figure28")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    gset <- getGEO("GSE14520", GSEMatrix =TRUE, getGPL=F)

    ## 表达矩阵提取
    exprSet <- exprs(gset[[1]])
    #     GSM362958 GSM362959 GSM362960
    # 1007_s_at     6.876     7.648     7.915
    # 1053_at       4.651     4.283     4.250
    # 117_at        6.775     3.796     3.380
    ## 分组信息等提取
    pData <- pData(gset[[1]])

    platform <- getGEO("GPL571")

    gene_table <- data.frame(gene = platform@dataTable@table$"Gene Symbol",
        id = platform@dataTable@table$ID)
    gene_array <- platform@dataTable@table$"Gene Symbol"
    names(gene_array) <- platform@dataTable@table$ID

    gene_symbol_array <- as.character(gene_array[row.names(exprSet)])
    row.names(exprSet) <- gene_symbol_array

    meta_info <- read.delim(file = "/mnt/share02/project/P06_public_data/Cancer/Liver/Roessler_GSE14520/raw/GSE14520_Extra_Supplement.txt")
    meta_info <- meta_info %>% dplyr::select(c("Affy_GSM","Survival.status","Survival.months","Recurr.status","Recurr.months"))
    meta_info$Survival.months <- meta_info$Survival.months / 12
    meta_info$Recurr.months <- meta_info$Recurr.months / 12
    # r$> head(meta_info)
    # Affy_GSM Survival.status Survival.months Recurr.status Recurr.months
    # 1 GSM363205               0            58.0             0          58.0
    # 2 GSM363115               0            66.6             0          66.6
    # 3 GSM362970               0            67.3             0          67.3
    # 4 GSM363354               1            10.4             1          10.4
    # 5 GSM363039               0            52.8             0          52.8
    # 6 GSM363209               0            60.8             0          60.8
    meta_info = meta_info %>% rename("sample" = "Affy_GSM")

    cpdb_signature <- "
        TAM	VEGFA
        TAM	CCL2
        TAM	CD81
        TAM	HLA-DRA
        M2	CD163
        M2	MSR1
        M2	MRC1
        Exhausted	CD8A
        Exhausted	CD8B
        Exhausted	PDCD1
        Exhausted	CTLA4
        Exhausted	TIGIT
        Exhausted	BTLA
        LR	CD86
        LR	HAVCR2
        LR	LGALS9
        LR	PVR
        LR	SPP1
        LR	SIRPA
        LR	CD47
        LR	CD44
        LR	CRTAM
        LR	TNFRSF14
        LR	ICAM1
        LR	VCAM1
        LR	PDCD1LG2
        LR	CXCL12
        LR	CXCR6
        LR	CXCR3
        LR	CXCL16
        LR	CXCR4
        LR	NECTIN2
        LR	MIF
        LR	CD74
        LR	CSF1
        LR	CSF1R"
    cpdb_signature_frame <- read.table(
        text = gsub("        ","",cpdb_signature),
        sep = "\t", 
        col.names = c("cluster", "marker")
    )
    

    # expression
    gene <- cpdb_signature_frame$marker
    gene <- intersect(gene,row.names(exprSet))
    data <- exprSet[gene,]
    data <- t(data)
    cpdb <- apply(data[,gene],1,mean)
    attr(cpdb,"names") <- NULL
    data <- cbind(data,cpdb)
    data <- data.frame(data)
    data$sample <- row.names(data)
    row.names(data) <- NULL
    


    surdata <- merge(meta_info,data,by = c("sample"))
    surdata$level <- ifelse(surdata[,"cpdb"] > median(surdata[,"cpdb"]),'High','Low')
    #surdata$level <- ifelse(surdata[,"mean_expr"] > as.numeric(quantile(surdata[,"mean_expr"],0.75)),'High','Low')
    fit <- survfit(Surv(Survival.months,Survival.status)~level,data = surdata)
    print(surv_pvalue(fit)$pval)

    p <- ggsurvplot(fit,
        pval = T, 
        pval.method = T, 
        risk.table = T,
        size = 2, # line size
        linetype = "strata", # line type by groups
        palette = c("#E7B800", "#2E9FDF"),
        conf.int = T # add confidence interval
    )
    
    png(file = "cpdb.survival.png",width = 1024,height = 768)
    print(p)
    dev.off()

    pdf(file = "cpdb.survival.pdf",width = 12,height = 8)
    print(p)
    dev.off()
}




# Figure29 Scissor analysis

if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure29")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)

    library(Seurat)
    # library(dior)
    library(Scissor)
    library(dplyr)

    tnk_sce <- readRDS("/P02_TNKCell/P999_TNK_Analysis/Figure0/TNK.rds")
    myeloid_sce <- readRDS("/P03_MyeloidCell/P03_Myeloid_UMAPByCellType/myeloid.rds")

    tnk_subcluster_frame <- tnk_sce@meta.data["clusters"]
    colnames(tnk_subcluster_frame) <- c("subclusters")
    myeloid_subcluster_frame <- myeloid_sce@meta.data["clusters"]
    colnames(myeloid_subcluster_frame) <- c("subclusters")
    tnk_myeloid_subclusters_frame <- rbind(tnk_subcluster_frame,myeloid_subcluster_frame)

    LIHC_expr <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.expr.rds")
    LIHC_surv <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.survival.rds")
    inter_sample_id_array <- intersect(colnames(LIHC_expr),LIHC_surv$sample)
    LIHC_expr <- LIHC_expr[,inter_sample_id_array]
    LIHC_surv <- LIHC_surv %>% filter(sample %in% inter_sample_id_array)
    # load(url(paste0(location, 'TCGA_LUAD_exp1.RData')))
    # load(url(paste0(location, 'TCGA_LUAD_survival.RData')))

    # saveRDS(bulk_dataset,file = "scissor_TCGA_LUAD_expr_demo_data.rds",compress = F)
    # saveRDS(bulk_survival,file = "scissor_TCGA_LUAD_survival_demo_data.rds",compress = F)

    LIHC_phenotype <- LIHC_surv[,c(3,2)]
    colnames(LIHC_phenotype) <- c("time", "status")
    LIHC_phenotype$time <- LIHC_phenotype$time * 365
    head(LIHC_phenotype)
    source("/mnt/share01/projects/scLiver/script/Scissor/as_matrix.R")
    source("/mnt/share01/projects/scLiver/script/Scissor/scissor_change.R")
    sample_array <- unique(as.character(global_sce$sample_id))
    
    alpha_array <- 2^(-(24:2)/2)
    for (sample in sample_array[grepl("T",sample_array)]) {
        message(sample)
        sample_sce <- subset(global_sce,(sample_id %in% c(sample)))
        message(sprintf("%s cells",length(Cells(sample_sce))))
        message("subset finished")
        LIHC_dataset <- Seurat_preprocessing(sample_sce@assays$RNA@counts, verbose = F)
        message("preprocessing finished")
        # It contains the required preprocessed matrix and constructed cell-cell similarity network, 
        # as well as other helpful dimensionality reduction results, such as the PCA, t-SNE, and UMAP.
        # print(names(LIHC_dataset))
        # infos1 <- Scissor_change(LIHC_expr, LIHC_dataset, LIHC_phenotype, alpha = alpha_array, cutoff = 0.3,
        #             family = "cox", Save_file = paste0(sample,'_Scissor_LIHC_survival.RData'))
        # infos1 <- Scissor_change(LIHC_expr, LIHC_dataset, LIHC_phenotype, alpha = alpha_array, cutoff = 0.3,
        #             family = "cox", Save_file = paste0(sample,'_Scissor_LIHC_survival.RData'))
        # infos1 <- Scissor_change(LIHC_expr, LIHC_dataset, LIHC_phenotype, alpha = alpha_array, cutoff = 0.3,
        #             family = "cox", Load_file = paste0(sample,'_Scissor_LIHC_survival.RData'))
        infos1 <- Scissor_change(LIHC_expr, 
            LIHC_dataset, 
            LIHC_phenotype, 
            alpha = alpha_array, 
            cutoff = 0.3,
            family = "cox")
        saveRDS(infos1, file = paste0(sample,'_Scissor_LIHC_survival.rds'))
        message("scissor finished")
        # saveRDS(infos1,file = paste0(sample_id,".scissor.rds"),compress = F)
    }

    # load("/P01_GlobalCell/P999_GlobalCell_Analysis/Figure29/LR01T_Scissor_LIHC_survival.RData")

    scissor_frame <- data.frame()
    for (s in sample_array[grepl("T",sample_array)]) {
        message(s)
        sci <- readRDS(paste0(s,"_Scissor_LIHC_survival.rds"))
        if (length(sci$Scissor_pos) > 0) {
            scissor_pos_metaframe <- data.frame(cell_id = sci$Scissor_pos)
            scissor_pos_metaframe$Scissor_status <- 1
            scissor_frame <- rbind(scissor_frame,scissor_pos_metaframe)
        }
        if (length(sci$Scissor_neg) > 0) {
            scissor_neg_metaframe <- data.frame(cell_id = sci$Scissor_neg)
            scissor_neg_metaframe$Scissor_status <- 2
            scissor_frame <- rbind(scissor_frame,scissor_neg_metaframe)
        }
    }
    
    row.names(scissor_frame) <- scissor_frame$cell_id
    scissor_frame$cell_id <- NULL
    global_sce <- AddMetaData(global_sce, metadata = scissor_frame)
    global_sce@meta.data$Scissor_status[is.na(global_sce@meta.data$Scissor_status)] <- 0
    gp <- DimPlot(global_sce, 
        reduction = 'umap', 
        group.by = 'Scissor_status', 
        cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
    ggsave(gp, file = "umap_scissor.png",w = 10,h = 10)

    

    celltype_array <- unique(as.character(global_sce$clusters))
    tumor_sce <- subset(global_sce,sample_id %in% sample_array[grepl("T",sample_array)])
    all_celltype_stat_frame <- data.frame()
    for (celltype in celltype_array) {
        message(celltype)
        celltype_sce <- subset(tumor_sce,clusters %in% c(celltype))

        prop_frame <- table(celltype_sce@meta.data$sample_id,celltype_sce@meta.data$Scissor_status) %>% prop.table(margin = 1)
        prop_frame <- prop_frame + 0.01
        log2fd <- log2(mean(prop_frame[,2]) / mean(prop_frame[,3]))
        #long_prop_frame <- prop_frame %>% as.data.frame() %>% filter(Var2 %in% c(1,2))
        wt <- wilcox.test(prop_frame[,2],prop_frame[,3])
        celltype_frame <- data.frame(celltype = c(celltype),log2fd = c(log2fd),pvalue = c(wt$p.value))
        all_celltype_stat_frame <- rbind(all_celltype_stat_frame,celltype_frame)
    }
    
    gp <- ggplot(all_celltype_stat_frame, aes(x=celltype, y=log2fd)) +
         geom_bar(stat="identity") +
         theme_bw()
    ggsave("global_scissor_barplot.png",w = 10,h = 6)

    tnk_myeloid_sce <- subset(global_sce,cells = row.names(tnk_myeloid_subclusters_frame))
    tnk_myeloid_sce <- AddMetaData(tnk_myeloid_sce, metadata = tnk_myeloid_subclusters_frame)


    celltype_array <- unique(as.character(tnk_myeloid_sce$subclusters))
    tnk_myeloid_sce <- subset(tnk_myeloid_sce,sample_id %in% sample_array[grepl("T",sample_array)])
    all_celltype_stat_frame <- data.frame()
    for (celltype in celltype_array) {
        message(celltype)
        celltype_sce <- subset(tnk_myeloid_sce,subclusters %in% c(celltype))

        prop_frame <- table(celltype_sce@meta.data$sample_id,celltype_sce@meta.data$Scissor_status) %>% prop.table(margin = 1)
        prop_frame <- prop_frame + 0.01
        if (ncol(prop_frame) < 3) {
            prop_frame <- cbind(prop_frame,0.01)
            colnames(prop_frame) <- c("0","1","2")
        }
        log2fd <- log2(mean(prop_frame[,2]) / mean(prop_frame[,3]))
        #long_prop_frame <- prop_frame %>% as.data.frame() %>% filter(Var2 %in% c(1,2))
        wt <- wilcox.test(prop_frame[,2],prop_frame[,3])
        celltype_frame <- data.frame(celltype = c(celltype),log2fd = c(log2fd),pvalue = c(wt$p.value))
        all_celltype_stat_frame <- rbind(all_celltype_stat_frame,celltype_frame)
    }
   
    gp <- ggplot(all_celltype_stat_frame, aes(x=celltype, y=log2fd)) +
         geom_bar(stat="identity") +
         theme_bw()
    ggsave("tnk_myeloid_scissor_barplot.png",w = 10,h = 6)
 

}



# Figure30 get Seurat matrix to csv

if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure30")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)

    library(data.table)

    fwrite(x = as.data.frame(global_sce[["RNA"]]@counts), 
        row.names=T,
        file = "LRT.csv")


}












