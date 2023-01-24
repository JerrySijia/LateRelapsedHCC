
# load library
library(Seurat)
library(ggplot2)
library(ggridges)
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)
library(ggsignif)
library(circlize)
library(future)
library(ggrepel)
#library("monocle3")
library(SingleCellExperiment)
library(AUCell)
library(msigdbr)
library(fgsea)


tnk_seurat_rds <- "/P02_TNKCell/P999_TNK_Analysis/Figure0/TNK.rds"
ref_seurat_rds <- "/P00_RefData/global_seurat.RDS" # nolint
sample_list <- "/P02_TNKCell/P999_TNK_Analysis/script/sample_list.txt"
celltype_group <- "/P01_GlobalCell/P999_GlobalCell_Analysis/celltype_group.txt"
hypoxia_siganture_file = "/signature/hypoxia_signature.txt"
tcell_signature_file = "/P02_TNKCell/P999_TNK_Analysis/script/TCell_signature_gene.txt"
outdir <- "/P02_TNKCell/P999_TNK_Analysis"

# Read TNK Seurat RDS Data
cat("Reading TNK Seurat RDS data...\n")
tnk_sce <- readRDS(tnk_seurat_rds)
cat("Complete.\n")

if (F) {
    cat("Reading Ref Seurat RDS data...\n")
    ref_sce <- readRDS(ref_seurat_rds)
    cat("Complete.\n")
}

plan("multicore", workers = 4)
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
    ZR04N = "#DCC1DD",
    ZR05T = "#CCE0F5",
    ZR05N = "#CCC9E6",
    ZR06T = "#625D9E",
    ZR06N = "#68A180",
    ZR07T = "#3A6963",
    ZR07N = "#968175"
)

group_cols <- c(
    ZPN = "#ef9020",
    ZPT = "#00af3e", 
    ZRN = "#0081b4", 
    ZRT = "#CD1E24"
)

group_cols <- c(
    PHT = "#00af3e", 
	PHN = "#ef9020",
    LRT = "#CD1E24",
    LRN = "#0081b4" 
)


celltype_colors <- c(
    "CD8-C1" = "#262C68",
    "CD8-C2" = "#CD1E24",
    "CD8-C3" = "#1E843F",
    "CD8-C4" = "#84278B",
    "CD8-C5" = "#0B6E78",
    "CD8-C6" = "#E4C755",
    "CD4-C1" = "#EF7A2A",
    "CD4-C2" = "#B968A5",
    "CD4-C3" = "#BDD240",
    "MAIT" = "#E63863",
    "T-Cycling" = "#D6E7A3",
    "NK-C1" = "#DCC1DD",
    "NK-C2" = "#58A4C3",
    "NK-C3" = "#AB3282"
)


#VIP Figure 0: Add Celltype information to Seurat Object
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure0")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)

    celltype <- "
        0	CD4-C1
        1	CD8-C1
        2	NK-C1
        3	MAIT
        4	CD8-C6
        5	NK-C2
        6	CD8-C2
        7	CD4-C2
        8	CD8-C3
        9	CD8-C4
        11	T-Cycling
        12	CD8-C5
        13	CD4-C3
        15	NK-C3
        16	NK-C1"

    # sce <- eval(parse(text = load(seurat_rdata)))
    tnk_sce <- readRDS(tnk_seurat_rds)

    sample_list <- data.frame(orig.ident = names(sample_colors))
    sample_list$sample_id <- gsub("ZP","PH",sample_list$orig.ident)
    sample_list$sample_id <- gsub("ZR","LR",sample_list$sample_id)
    sample_list$group_id <- paste0(substr(sample_list$sample_id,1,2),substr(sample_list$sample_id,5,5))

    meta_frame <- tnk_sce@meta.data %>% rownames_to_column("cell_id") %>%  merge(sample_list,by = "orig.ident") %>% select(c("sample_id","group_id","cell_id"))
    row.names(meta_frame) <- meta_frame$cell_id
    meta_frame$cell_id <- NULL
    tnk_sce[["sample_id"]] <- NULL
    tnk_sce[["group_id"]] <- NULL
    # 这里只要有cell_id作为row.names就可以用AddMetaData进行合并
    tnk_sce <- AddMetaData(tnk_sce,meta_frame)
    saveRDS(tnk_sce,file = "TNK.rds",compress = F)

    # 根据cell type 文件（第一列cluster_name,第二列细胞类型）标注细胞类型到seurat object
    cell_type_frame <- read.table(
        text = celltype,
        sep = "\t",
        col.names = c("index", "cell_type")
    )

    # get the cluster index and give each index corresponding cell types
    tnk_sce <- tnk_sce[,(tnk_sce@meta.data$seurat_clusters %in% cell_type_frame$index)]
    # 取第二列转成character类型的向量
    new.cluster.ids <- as.character(cell_type_frame$cell_type)
    # 
    names(new.cluster.ids) <- levels(tnk_sce)
    tnk_sce <- RenameIdents(tnk_sce, new.cluster.ids)
    tnk_sce@meta.data$clusters = tnk_sce@active.ident

    data = tnk_sce@reductions$umap@cell.embeddings %>% as.data.frame() %>% cbind(clusters = tnk_sce@meta.data$clusters)

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
        geom_text(aes(label = clusters), data = class_avg, size = 5)+
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

    ggsave("Figure0.pdf",plot = umap,width = 12,height = 8)
    ggsave("Figure0.png",plot = umap,width = 12,height = 8)

    # save
    saveRDS(tnk_sce, file = paste0('TNK.rds'), compress = FALSE)

}


# Figure 1: Plot TNK UMAP by group
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure1")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    setwd(fig_outdir)

    ## add group information

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

    tnk_sce@meta.data$group <- as.character(group_array[as.character(tnk_sce$orig.ident)])

    ### do not use DimPlot Function, change to ggplot2
    ### step1 get umap coordinates
    data = tnk_sce@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        cbind(group = tnk_sce@meta.data$group)

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
            plot.margin=unit(rep(1,4),'cm'),
            legend.text = element_text(size = 18),
            legend.key.size = unit(0.4, "inches")
        ) +
        guides(colour = guide_legend(override.aes = list(size = 6)))

    
    ggsave("Figure1.pdf",plot = umap,width = 12,height = 12)
    ggsave("Figure1.png",plot = umap,width = 12,height = 12)

}


# Figure 2: Plot TNK UMAP by sample
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure2")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    ## Plot Figures
    data = tnk_sce@reductions$umap@cell.embeddings %>%
        as.data.frame() %>%
        cbind(clusters = tnk_sce@meta.data$orig.ident)

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


#VIP Figure 3: Plot TNK Marker Gene Heatmap By Celltype
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure3")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    tnk_cellmarker = "
		Naive	TCF7
		Naive	LEF1
		Naive	CCR7
		Naive	SELL
		Naive	MAL
		T cell markers	CD3D
		T cell markers	CD3E
		T cell markers	CD4
		T cell markers	CD8A
		T cell markers	CD8B
		Proliferation	MKI67
		Proliferation	CDK1
		Proliferation	STMN1
		CD4-FOXP3	FOXP3
		Proliferation	STMN1
		Proliferation	STMN1
		NK cell markers	NCAM1
		NK cell markers	FCGR3A
		MAIT markers	SLC4A10
		MAIT markers	RORC"

	cell_marker_frame <- read.table(
        text = gsub("\t\t","",tnk_cellmarker),
        sep = "\t", 
        col.names = c("celltype", "gene")
    )


    Idents(tnk_sce) <- "clusters"

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$gene),]

    AverageExp <- AverageExpression(tnk_sce,
        features = cell_marker_frame$gene, 
        group.by = "clusters"
    )

    expr <- AverageExp$RNA

    heatmapColorRamp = colorRamp2(c(-3, 0, 3), c("#6E6BA9","#FDF8EE","#9D1F28"))

    expr <- ScaleData(expr)
    color_array <- c("#9F1F5C",
                    "#EF9020",
                    "#00AF3E",
                    "#29245C",
                    "#E5352B",
                    "#FFD616",
                    "#E990AB",
                    "#0081B4")

    expr <- expr[, sort(colnames(expr))]

    celltype_group <- c()
    for (celltype in colnames(expr)) {
        celltype_group <- c(celltype_group, unlist(strsplit(celltype, "-"))[1])
    }
    

    celltype_color_structure <- structure(color_array[1:(length(unique(cell_marker_frame$celltype)))],
                                        names = unique(cell_marker_frame$celltype)) 

    ha = rowAnnotation(labels = cell_marker_frame$celltype,
					show_legend = F,
                    col = list(bar = celltype_color_structure))

    ht <- Heatmap(as.matrix(expr),
            col = heatmapColorRamp,
            right_annotation = ha,
            rect_gp = gpar(col = "white", lwd = 1),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 16),
            column_split = celltype_group,
            # remove the title showed on the top of the heatmap (due to column_split will show the title)
            column_title = NULL,
            row_split = cell_marker_frame$celltype,
            row_title_side = "right",
            row_title_rot = 0, 
            row_title_gp = gpar(fontsize = 25),
			heatmap_legend_param = list(title = "",
				grid_width = unit(1, "cm"),
				legend_height = unit(8, "cm"),
                title_gp = gpar(fontsize = 20), 
                labels_gp = gpar(fontsize = 30)),
            cluster_columns = F,
            cluster_rows = F)

    png("TNKMarkerGeneHeatmapByCelltype.png",width = 1080,height = 1080 )
    draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

    pdf("TNKMarkerGeneHeatmapByCelltype.pdf",width = 12,height = 12)
    draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

}


# Figure 4: Plot TNK UMAP by marker
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure4")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    marker_frame <- read.table(file = cell_marker, sep = "\t", header = T, stringsAsFactors = F)

    marker_list <- split(marker_frame,marker_frame$celltype)
    for (marker in marker_list){
        cell <- marker[, 1][1]
        cat("Parsing ",cell,"\n")
        markers <- marker[, 2]
        
        plot_array <- c()
        
        p1 <- FeaturePlot(tnk_sce,
            features = markers, 
            cols = c("lightgrey", "#e32119"), 
            pt.size = 0.5,
            raster = F
        )

        ggsave(plot = p1, filename = paste0(cell, "_Figure4.pdf"), w = 12, h = 8, units = c("in"))
        ggsave(plot = p1, filename = paste0(cell, "_Figure4.png"), w = 12, h = 8, units = c("in"))
    }

}

#Figure 5: Plot Cell type stat BarPlot
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure5")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cell_stat_frame <- as.data.frame.array(table(tnk_sce$orig.ident, tnk_sce$clusters))

    write.table(
        x = cell_stat_frame,
        file = "global_celltype_stat.xls",
        sep = "\t",
        col.name = TRUE,
        row.name = TRUE,
        quote=FALSE
    )

    cell_stat_frame$group <- paste0(
        substr(cell_stat_frame$Sample, 0, 2),
        substr(cell_stat_frame$Sample, 5, 6)
    )

    cell_stat_frame %>%
        group_by(group) %>%
        summarise(
            Frequency = sum(Total_cell),
            Proportion = sum(Total_cell) / sum(cell_stat_frame$Total_cell)
        )
}




# Figure 6: CD8 frequency BoxPlot by ZPT, ZRT, ZPN and ZRN
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure6")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)



    cell_stat_frame <- as.data.frame.array(table(
        as.character(tnk_sce$sample_id),
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
                #test.args = list(var.equal = TRUE,alternative = "greater"), # one side or two side?
                test = wilcox.test) +
    facet_grid(~variable)

    ggsave("CD8FrequencyBoxPlot.png",p1,width = 10,height = 4.5,units = "in")
    ggsave("CD8FrequencyBoxPlot.pdf",p1,width = 10,height = 4.5,units = "in")

}



# Figure 7: CD8 hypoxia signature VlnPlot and DotPlot for each celltype
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure7")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    hypoxia_signature_frame <- read.table(file = hypoxia_siganture_file, sep = "\t",header = F,stringsAsFactors = F)
    tnk_sce <- AddModuleScore(object = tnk_sce, list(hypoxia_signature_frame$V2),name = "hypoxia_score")
    
    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)
    
    plot <- VlnPlot(cd8_sce,
        features = "hypoxia_score1",
        group.by = "clusters",
        col = celltype_colors
    )

    ggsave("Figure7-VlnPlot.png",plot = plot,width = 8,height = 4.5,units = "in")
    ggsave("Figure7-VlnPlot.pdf",plot = plot,width = 8,height = 4.5,units = "in")


    p <- DotPlot(
        object = cd8_sce,
        features = "hypoxia_score1",
        dot.scale = 12,
        cols = c("lightgrey", "red"),
        group.by = "clusters",
    ) + theme(axis.text.x = element_text(angle = 90)) + theme_bw()

    ggsave("Figure7-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure7-DotPlot.pdf",p,width = 8,height = 8,units = "in")

}


# Figure 8 Find Anchors for each cell type
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure8")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)
    ref_tnk_sce <- subset(ref_sce,celltype_global %in% c("Lymphoid-NK","Lymphoid-T"))
    anchors <- FindTransferAnchors(reference = ref_tnk_sce, query = tnk_sce, dims = 1:30)
    predictions <- TransferData(anchorset = anchors, refdata = ref_tnk_sce$celltype_sub, dims = 1:30)

    tnk_sce <- AddMetaData(tnk_sce, metadata = predictions)
    my36colors <- c(
        "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87",
        "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658",
        "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398",
        "#AA9A59", "#E63863", "#E39A35", "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B",
        "#712820", "#DCC1DD", "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963", "#968175"
    )
    colors.num <- length(levels(tnk_sce))
    colors <- DiscretePalette(colors.num, palette = "glasbey")
    plot1 <- DimPlot(tnk_sce,
        reduction = "umap",
        label = TRUE, 
        label.size = 5, 
        pt.size = 0.4,
        group.by = "predicted.id", 
        cols = my36colors
    )

    pdf(paste0("Figure8_UMAP.pdf"),w=16,h=12)
    print(plot1)
    dev.off()
    png(paste0("Figure8_UMAP.png"),width=1600,height=1200)
    print(plot1)
    dev.off()

    celltype_matrix <- table(tnk_sce$clusters,tnk_sce$predicted.id)
    for (row in 1:nrow(celltype_matrix)) {
        sum <- sum(celltype_matrix[row, ])
        celltype_matrix[row, ] <- celltype_matrix[row, ] / sum
    }

    png("Figure8_Heatmap.png",width = 1024,height = 1024)
    col_fun = colorRamp2(c(0, 0.2, 1), c("blue", "white", "red"))
    Heatmap(celltype_matrix,
        rect_gp = gpar(col = "white", lwd = 2),
        col = col_fun
    )
    dev.off()

}


# Figure 9 CD8 T cell signature VlnPlot and DotPlot for each cell type
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure9")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

	cd8_signature <- "
		Cytotoxic	PRF1
		Cytotoxic	IFNG
		Cytotoxic	GNLY
		Cytotoxic	NKG7
		Cytotoxic	GZMB
		Cytotoxic	GZMA
		Cytotoxic	HZMH
		Cytotoxic	KLRK1
		Cytotoxic	KLRKB1
		Cytotoxic	KLRD1
		Cytotoxic	CTSW
		Cytotoxic	CST7
		Cytotoxic	CX3CR1
		Cytotoxic	FGFBP2
		Cytotoxic	S1PR5
		Cytotoxic	FCGR3A
		Cytotoxic	PLAC8
		Cytotoxic	C1orf21
		Cytotoxic	TGFBR3
		Cytotoxic	PLEK
		Cytotoxic	FGR
		Cytotoxic	KLRF1
		Cytotoxic	SPON2
		Cytotoxic	CD300A
		Cytotoxic	S1PR1
		Cytotoxic	FAM65B
		Cytotoxic	EFHD2
		Cytotoxic	STK38
		Cytotoxic	C1orf162
		Cytotoxic	SORL1
		Cytotoxic	EMP3
		Cytotoxic	ARL4C
		Cytotoxic	BIN2
		Cytotoxic	CCND3
		Cytotoxic	FCRL6
		Cytotoxic	SAMD3
		Cytotoxic	TRDC
		Cytotoxic	TYROBP
		Cytotoxic	GNLY
		Cytotoxic	KLRG1
		Exhaustion	FCRL3
		Exhaustion	CD27
		Exhaustion	PRKCH
		Exhaustion	B2M
		Exhaustion	ITM2A
		Exhaustion	TIGIT
		Exhaustion	ID3
		Exhaustion	GBP2
		Exhaustion	PDCD1
		Exhaustion	KLRK1
		Exhaustion	HSPA1A
		Exhaustion	SRGN
		Exhaustion	TNFRSF9
		Exhaustion	TMBIM6
		Exhaustion	TNFRSF1B
		Exhaustion	CADM1
		Exhaustion	ACTB
		Exhaustion	CD8A
		Exhaustion	RGS2
		Exhaustion	FAIM3
		Exhaustion	EID1
		Exhaustion	HSPB1
		Exhaustion	RNF19A
		Exhaustion	IFI16
		Exhaustion	LYST
		Exhaustion	PRF1
		Exhaustion	STAT1
		Exhaustion	UBC
		Exhaustion	CD74
		Exhaustion	IL2RG
		Exhaustion	FYN
		Exhaustion	PTPN6
		Exhaustion	HLA-DRB1
		Exhaustion	HNRNPC
		Exhaustion	UBB
		Exhaustion	CD8B
		Exhaustion	HAVCR2
		Exhaustion	IRF8
		Exhaustion	LAG3
		Exhaustion	ATP5B
		Exhaustion	STAT3
		Exhaustion	IGFLR1
		Exhaustion	MGEA5
		Exhaustion	HSPA1B
		Exhaustion	COTL1
		Exhaustion	VCAM1
		Exhaustion	HLA-DMA
		Exhaustion	PDE7B
		Exhaustion	TBC1D4
		Exhaustion	SNAP47
		Exhaustion	RGS4
		Exhaustion	CBLB
		Exhaustion	TOX
		Exhaustion	CALM2
		Exhaustion	ATHL1
		Exhaustion	SPDYE5
		Exhaustion	DDX5
		Exhaustion	SLA
		Exhaustion	PTPRCAP
		Exhaustion	IRF9
		Exhaustion	MATR3
		Exhaustion	LITAF
		Exhaustion	TPI1
		Exhaustion	ETV1
		Exhaustion	PAM
		Exhaustion	ARID4B
		Exhaustion	NAB1
		Exhaustion	RAPGEF6
		Exhaustion	LDHA
		Exhaustion	WARS
		Exhaustion	RASSF5
		Exhaustion	OSBPL3
		Exhaustion	FAM3C
		Exhaustion	TAP1
		Exhaustion	HLA-DRB6
		Exhaustion	FABP5
		Exhaustion	CD200
		Exhaustion	CTLA4
		Exhaustion	SNX9
		Exhaustion	ETNK1
		Exhaustion	MALAT1
		Exhaustion	ZDHHC6
		Exhaustion	ARL6IP5
		Exhaustion	DUSP2
		Exhaustion	HLA-DQB1
		Exhaustion	HNRNPK
		Exhaustion	DGKH
		Exhaustion	LRMP
		Exhaustion	H3F3B
		Exhaustion	IDH2
		Exhaustion	TRAF5
		Exhaustion	TBL1XR1
		Exhaustion	ANKRD10
		Exhaustion	ALDOA
		Exhaustion	LSP1
		Exhaustion	PTPN7
		Exhaustion	NSUN2
		Exhaustion	RNF149
		Exhaustion	CD2
		Exhaustion	SRSF1
		Exhaustion	GOLPH3
		Exhaustion	HLA-A
		Exhaustion	LIMS1
		Exhaustion	SDF4
		Exhaustion	ROCK1
		Exhaustion	EDEM1
		Exhaustion	APLP2
		Exhaustion	ITK
		Exhaustion	TRIM22
		Exhaustion	SPRY2
		Exhaustion	ACTG1
		Exhaustion	HLA-DPA1
		Exhaustion	EWSR1
		Exhaustion	SRSF4
		Exhaustion	ESYT1
		Exhaustion	LUC7L3
		Exhaustion	ARNT
		Exhaustion	GNAS
		Exhaustion	ARF6
		Exhaustion	ARPC5L
		Exhaustion	NCOA3
		Exhaustion	PAPOLA
		Exhaustion	GFOD1
		Exhaustion	GPR174
		Exhaustion	DDX3X
		Exhaustion	CAPRIN1
		Exhaustion	ARPC2
		Exhaustion	PDIA6
		Exhaustion	SEMA4A
		Exhaustion	CSDE1
		Exhaustion	PSMB9
		Exhaustion	NFATC1
		Exhaustion	PTPN11
		Exhaustion	AGFG1
		Exhaustion	PCED1B
		Exhaustion	CCL4L1
		Exhaustion	CCND2
		Exhaustion	CCL4L2
		Exhaustion	CXCR6
		Exhaustion	AKAP5
		Exhaustion	IFNG
		Exhaustion	MIR155HG
		Exhaustion	ENTPD1
		Exhaustion	TOX2
		Exhaustion	CD7
		Exhaustion	RAB27A
		Exhaustion	ITGAE
		Exhaustion	PHLDA1
		Exhaustion	PAG1
		Exhaustion	CSF1
		Exhaustion	NBL1
		Exhaustion	CCL3
		Exhaustion	ILRB
		Exhaustion	FASLG
		Exhaustion	ZC3H12C
		Exhaustion	MYO7A
		Exhaustion	SIRPG
		Exhaustion	GALNT1
		Exhaustion	UBE2F
		Exhaustion	DUSP4
		Exhaustion	SYT11
		Exhaustion	TRAC
		Exhaustion	TNS3
		Exhaustion	RDH10
		Exhaustion	PTMS
		Exhaustion	CXCL13
		Exhaustion	KIR2DL4
		Naive	CCR7
		Naive	TCF7
		Naive	LEF1
		Naive	SELL"

		
	tcell_signature_frame <- read.table(
		text = gsub("\t\t","",cd8_signature),
		sep = "\t", 
		col.names = c("V1", "V2")
	)

	cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(tcell_signature_frame, tcell_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        cd8_sce <- AddModuleScore(
            object = cd8_sce,
            features = list(gene_array),
            name = signature
        )
        colnames(cd8_sce@meta.data) <- gsub(colnames(cd8_sce@meta.data),
            pattern = paste0(signature, 1),
			replacement = signature
        )
    }

    p <- DotPlot(
        object = cd8_sce,
        features = names(split_signature_list),
        dot.scale = 12,
        cols = c("#0E01FA", "#FB010F"),
        group.by = "clusters",
    ) + theme_bw() + 
	theme(axis.text.x = element_text(size = 15),
			axis.text.y = element_text(size = 15)) +
        xlab("") +
        ylab("")

    ggsave("Figure9-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure9-DotPlot.pdf",p,width = 8,height = 8,units = "in")


    plot <- VlnPlot(cd8_sce,
        features = names(split_signature_list),
        group.by = "clusters",
        col = celltype_colors,
        ncol = 2
    )
    
    ggsave("Figure9-VlnPlot.png",plot = plot,width = 8,height = 4.5,units = "in")
    ggsave("Figure9-VlnPlot.pdf",plot = plot,width = 8,height = 4.5,units = "in")


}


#VIP Figure 10 CD4 T cell signature VlnPlot and DotPlot for each cell type
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure10")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd4_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD4", tnk_sce@meta.data$clusters)]))

	cd4_signature <- "
		Cytotoxic	PRF1
		Cytotoxic	IFNG
		Cytotoxic	GNLY
		Cytotoxic	NKG7
		Cytotoxic	GZMB
		Cytotoxic	GZMA
		Cytotoxic	HZMH
		Cytotoxic	KLRK1
		Cytotoxic	KLRKB1
		Cytotoxic	KLRD1
		Cytotoxic	CTSW
		Cytotoxic	CST7
		Cytotoxic	CX3CR1
		Cytotoxic	FGFBP2
		Cytotoxic	S1PR5
		Cytotoxic	FCGR3A
		Cytotoxic	PLAC8
		Cytotoxic	C1orf21
		Cytotoxic	TGFBR3
		Cytotoxic	PLEK
		Cytotoxic	FGR
		Cytotoxic	KLRF1
		Cytotoxic	SPON2
		Cytotoxic	CD300A
		Cytotoxic	S1PR1
		Cytotoxic	FAM65B
		Cytotoxic	EFHD2
		Cytotoxic	STK38
		Cytotoxic	C1orf162
		Cytotoxic	SORL1
		Cytotoxic	EMP3
		Cytotoxic	ARL4C
		Cytotoxic	BIN2
		Cytotoxic	CCND3
		Cytotoxic	FCRL6
		Cytotoxic	SAMD3
		Cytotoxic	TRDC
		Cytotoxic	TYROBP
		Cytotoxic	GNLY
		Cytotoxic	KLRG1
		Exhaustion	FCRL3
		Exhaustion	CD27
		Exhaustion	PRKCH
		Exhaustion	B2M
		Exhaustion	ITM2A
		Exhaustion	TIGIT
		Exhaustion	ID3
		Exhaustion	GBP2
		Exhaustion	PDCD1
		Exhaustion	KLRK1
		Exhaustion	HSPA1A
		Exhaustion	SRGN
		Exhaustion	TNFRSF9
		Exhaustion	TMBIM6
		Exhaustion	TNFRSF1B
		Exhaustion	CADM1
		Exhaustion	ACTB
		Exhaustion	CD8A
		Exhaustion	RGS2
		Exhaustion	FAIM3
		Exhaustion	EID1
		Exhaustion	HSPB1
		Exhaustion	RNF19A
		Exhaustion	IFI16
		Exhaustion	LYST
		Exhaustion	PRF1
		Exhaustion	STAT1
		Exhaustion	UBC
		Exhaustion	CD74
		Exhaustion	IL2RG
		Exhaustion	FYN
		Exhaustion	PTPN6
		Exhaustion	HLA-DRB1
		Exhaustion	HNRNPC
		Exhaustion	UBB
		Exhaustion	CD8B
		Exhaustion	HAVCR2
		Exhaustion	IRF8
		Exhaustion	LAG3
		Exhaustion	ATP5B
		Exhaustion	STAT3
		Exhaustion	IGFLR1
		Exhaustion	MGEA5
		Exhaustion	HSPA1B
		Exhaustion	COTL1
		Exhaustion	VCAM1
		Exhaustion	HLA-DMA
		Exhaustion	PDE7B
		Exhaustion	TBC1D4
		Exhaustion	SNAP47
		Exhaustion	RGS4
		Exhaustion	CBLB
		Exhaustion	TOX
		Exhaustion	CALM2
		Exhaustion	ATHL1
		Exhaustion	SPDYE5
		Exhaustion	DDX5
		Exhaustion	SLA
		Exhaustion	PTPRCAP
		Exhaustion	IRF9
		Exhaustion	MATR3
		Exhaustion	LITAF
		Exhaustion	TPI1
		Exhaustion	ETV1
		Exhaustion	PAM
		Exhaustion	ARID4B
		Exhaustion	NAB1
		Exhaustion	RAPGEF6
		Exhaustion	LDHA
		Exhaustion	WARS
		Exhaustion	RASSF5
		Exhaustion	OSBPL3
		Exhaustion	FAM3C
		Exhaustion	TAP1
		Exhaustion	HLA-DRB6
		Exhaustion	FABP5
		Exhaustion	CD200
		Exhaustion	CTLA4
		Exhaustion	SNX9
		Exhaustion	ETNK1
		Exhaustion	MALAT1
		Exhaustion	ZDHHC6
		Exhaustion	ARL6IP5
		Exhaustion	DUSP2
		Exhaustion	HLA-DQB1
		Exhaustion	HNRNPK
		Exhaustion	DGKH
		Exhaustion	LRMP
		Exhaustion	H3F3B
		Exhaustion	IDH2
		Exhaustion	TRAF5
		Exhaustion	TBL1XR1
		Exhaustion	ANKRD10
		Exhaustion	ALDOA
		Exhaustion	LSP1
		Exhaustion	PTPN7
		Exhaustion	NSUN2
		Exhaustion	RNF149
		Exhaustion	CD2
		Exhaustion	SRSF1
		Exhaustion	GOLPH3
		Exhaustion	HLA-A
		Exhaustion	LIMS1
		Exhaustion	SDF4
		Exhaustion	ROCK1
		Exhaustion	EDEM1
		Exhaustion	APLP2
		Exhaustion	ITK
		Exhaustion	TRIM22
		Exhaustion	SPRY2
		Exhaustion	ACTG1
		Exhaustion	HLA-DPA1
		Exhaustion	EWSR1
		Exhaustion	SRSF4
		Exhaustion	ESYT1
		Exhaustion	LUC7L3
		Exhaustion	ARNT
		Exhaustion	GNAS
		Exhaustion	ARF6
		Exhaustion	ARPC5L
		Exhaustion	NCOA3
		Exhaustion	PAPOLA
		Exhaustion	GFOD1
		Exhaustion	GPR174
		Exhaustion	DDX3X
		Exhaustion	CAPRIN1
		Exhaustion	ARPC2
		Exhaustion	PDIA6
		Exhaustion	SEMA4A
		Exhaustion	CSDE1
		Exhaustion	PSMB9
		Exhaustion	NFATC1
		Exhaustion	PTPN11
		Exhaustion	AGFG1
		Exhaustion	PCED1B
		Exhaustion	CCL4L1
		Exhaustion	CCND2
		Exhaustion	CCL4L2
		Exhaustion	CXCR6
		Exhaustion	AKAP5
		Exhaustion	IFNG
		Exhaustion	MIR155HG
		Exhaustion	ENTPD1
		Exhaustion	TOX2
		Exhaustion	CD7
		Exhaustion	RAB27A
		Exhaustion	ITGAE
		Exhaustion	PHLDA1
		Exhaustion	PAG1
		Exhaustion	CSF1
		Exhaustion	NBL1
		Exhaustion	CCL3
		Exhaustion	ILRB
		Exhaustion	FASLG
		Exhaustion	ZC3H12C
		Exhaustion	MYO7A
		Exhaustion	SIRPG
		Exhaustion	GALNT1
		Exhaustion	UBE2F
		Exhaustion	DUSP4
		Exhaustion	SYT11
		Exhaustion	TRAC
		Exhaustion	TNS3
		Exhaustion	RDH10
		Exhaustion	PTMS
		Exhaustion	CXCL13
		Exhaustion	KIR2DL4
		Naive	CCR7
		Naive	TCF7
		Naive	LEF1
		Naive	SELL
		Treg	NT5E
		Treg	CD3D
		Treg	CD3G
		Treg	CD3E
		Treg	CD4
		Treg	CD5
		Treg	ENTPD1
		Treg	CTLA4
		Treg	IZUMO1R
		Treg	TNFRSF18
		Treg	IL2RA
		Treg	ITGAE
		Treg	LAG3
		Treg	TGFB1
		Treg	LRRC32
		Treg	TNFRSF4
		Treg	SELL
		Treg	FOXP3
		Treg	STAT5A
		Treg	STAT5B
		Treg	LGALS1
		Treg	IL10
		Treg	IL12A
		Treg	EBI3
		Treg	TGFB1"

    tcell_signature_frame <- read.table(
		text = gsub("\t\t","",cd4_signature),
		sep = "\t", 
		col.names = c("V1", "V2")
	)

    cd4_sce <- subset(tnk_sce, clusters %in% cd4_celltype_array)

    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(tcell_signature_frame, tcell_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        cd4_sce <- AddModuleScore(
            object = cd4_sce,
            features = list(gene_array),
            name = signature
        )
        colnames(cd4_sce@meta.data) <- gsub(colnames(cd4_sce@meta.data),
			pattern = paste0(signature, 1),
			replacement = signature
		)
    }
    
    p <- DotPlot(
        object = cd4_sce,
        features = names(split_signature_list),
        dot.scale = 12,
        cols = c("#0E01FA", "#FB010F"),
        group.by = "clusters",
    ) + theme_bw() + theme(axis.text.x = element_text(size = 15),
			axis.text.y = element_text(size = 15)) +
            xlab("") +
            ylab("")

    ggsave("Figure10-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure10-DotPlot.pdf",p,width = 8,height = 8,units = "in")


    plot <- VlnPlot(cd4_sce,
        features = paste0(names(split_signature_list),"1"),
        group.by = "clusters",
        col = celltype_colors,
        ncol = 2
    )
    
    ggsave("Figure10-VlnPlot.png",plot = plot,width = 8,height = 4.5,units = "in")
    ggsave("Figure10-VlnPlot.pdf",plot = plot,width = 8,height = 4.5,units = "in")


}


#VIP Figure 11 CD4 frequency BoxPlot by ZPT, ZRT, ZPN and ZRN
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure11")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

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
    cd4_celltype_array = colnames(cell_stat_frame %>% select(starts_with("CD4")))
    cell_stat_frame <- cell_stat_frame[, c(cd4_celltype_array, "Sample", "sample_group")]

    cell_stat_frame[c("Total_cell")] <- NULL

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

    ggsave("Figure11.png",p1,width = 8,height = 4.5,units = "in")
    ggsave("Figure11.pdf",p1,width = 8,height = 4.5,units = "in")

}


#VIP Figure 12 NK frequency BoxPlot by ZPT, ZRT, ZPN and ZRN
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure12")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

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
    cd4_celltype_array = colnames(cell_stat_frame %>% dplyr::select(starts_with("NK")))
    cell_stat_frame <- cell_stat_frame[, c(cd4_celltype_array, "Sample", "sample_group")]

    cell_stat_frame[c("Total_cell")] <- NULL

    melt_cell_stat_frame = melt(cell_stat_frame,id = c("Sample","sample_group"))
    melt_cell_stat_frame$sample_group <- factor(melt_cell_stat_frame$sample_group,levels = c("PHT","PHN","LRT","LRN"))
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
                map_signif_level = F,
                #test.args = list(var.equal = TRUE,alternative = "greater"), # one side or two side?
                test = wilcox.test) +
    facet_grid(~variable)

    ggsave("Figure12.png",p1,width = 8,height = 4.5,units = "in")
    ggsave("Figure12.pdf",p1,width = 8,height = 4.5,units = "in")

}


# Figure 13 CD8 Vocalno Plot by ZPT and ZRT
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure13")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)
    cd8_tumor_sce <- subset(cd8_sce, group %in% c("ZPT","ZRT"))
    diff_marker_frame <- FindMarkers(cd8_tumor_sce,
        ident.1 = "ZPT", 
        group.by = "group", 
        assay = "RNA", 
        slot = "data", # must use data not count
        logfc.threshold = 0, 
        min.pct = 0,
        test.use="DESeq2"
    )

    write.table(
        x = diff_marker_frame,
        file = "cd8_tumor_diffexp.gene.txt",
        sep = "\t",
        col.name = TRUE,
        row.name = F,
        quote=FALSE
    )

    log2FC_threshold = 0.8 
    p_value_threshold = 0.01
    diff_marker_frame[which(diff_marker_frame$p_val_adj < 0.01 & diff_marker_frame$avg_log2FC <= -log2FC_threshold),'sig'] <- 'CD8+ ZRT'
    diff_marker_frame[which(diff_marker_frame$p_val_adj < 0.01 & diff_marker_frame$avg_log2FC >= log2FC_threshold),'sig'] <- 'CD8+ ZPT'
    diff_marker_frame[which(diff_marker_frame$p_val_adj >= 0.01 | abs(diff_marker_frame$avg_log2FC) < log2FC_threshold),'sig'] <- 'None'
    # 横轴 log2FC，纵轴 -log10(adj.P.Val)，颜色表示差异

    diff_marker_frame$Gene <- row.names(diff_marker_frame)

    p <- ggplot(diff_marker_frame, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
        geom_point(alpha = 0.8, size = 0.6) +
        scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('CD8+ ZPT', 'CD8+ ZRT', 'None')) +
        theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
        theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
        geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), color = 'gray', size = 0.3) +
        geom_hline(yintercept = -log(0.01, 10), color = 'gray', size = 0.3) +
        xlim(-2, 2) + 
        ylim(0, 150) +
        labs(x = '\nLog2 Fold Change', y = '-log10(pvalue)\n', color = '', title = 'ZPT vs ZRT\n')

    up <- subset(diff_marker_frame, sig == 'CD8+ ZPT')
    up <- up[order(up$p_val_adj), ][1:40, ]
    down <- subset(diff_marker_frame, sig == 'CD8+ ZRT')
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

    ggsave('Figure13.png', p1, width = 5, height = 6)

}


# Figure 14 CD8 Vocalno Plot by ZPN and ZRN
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure14")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)
    cd8_normal_sce <- subset(cd8_sce, group %in% c("ZPN", "ZRN"))
    
    allCells = colnames(cd8_normal_sce)
    allType = names(table(cd8_normal_sce@meta.data$orig.ident))

    choose_Cells = unlist(lapply(allType, function(x){
        cgCells = allCells[cd8_normal_sce@meta.data$orig.ident == x]
        num = ceiling(dim(as.data.frame(cgCells))[1]*0.3)
        cg = sample(cgCells,num)
        cg  
    }))

    cd8_normal_sce = cd8_normal_sce[, allCells %in% choose_Cells]

    diff_marker_frame <- FindMarkers(cd8_normal_sce,
        ident.1 = "ZPN", 
        group.by = "group", 
        assay = "RNA", 
        slot = "data", 
        logfc.threshold = 0, 
        min.pct = 0
    )

    write.table(
        x = diff_marker_frame,
        file = "cd8_normal_diffexp.gene.txt",
        sep = "\t",
        col.name = TRUE,
        row.name = F,
        quote=FALSE
    )

    log2FC_threshold = 0.8 
    p_value_threshold = 0.01
    diff_marker_frame[which(diff_marker_frame$p_val_adj < 0.01 & diff_marker_frame$avg_log2FC <= -log2FC_threshold),'sig'] <- 'CD8+ ZRN'
    diff_marker_frame[which(diff_marker_frame$p_val_adj < 0.01 & diff_marker_frame$avg_log2FC >= log2FC_threshold),'sig'] <- 'CD8+ ZPN'
    diff_marker_frame[which(diff_marker_frame$p_val_adj >= 0.01 | abs(diff_marker_frame$avg_log2FC) < log2FC_threshold),'sig'] <- 'None'
    # 横轴 log2FC，纵轴 -log10(adj.P.Val)，颜色表示差异

    diff_marker_frame$Gene <- row.names(diff_marker_frame)

    p <- ggplot(diff_marker_frame, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
        geom_point(alpha = 0.8, size = 0.6) +
        scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('CD8+ ZPN', 'CD8+ ZRN', 'None')) +
        theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
        theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
        geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), color = 'gray', size = 0.3) +
        geom_hline(yintercept = -log(0.01, 10), color = 'gray', size = 0.3) +
        xlim(-2, 2) + 
        ylim(0, 10) +
        labs(x = '\nLog2 Fold Change', y = '-log10(pvalue)\n', color = '', title = 'ZPN vs ZRN\n')

    up <- subset(diff_marker_frame, sig == 'CD8+ ZPN')
    up <- up[order(up$p_val_adj), ][1:40, ]
    down <- subset(diff_marker_frame, sig == 'CD8+ ZRN')
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

    ggsave('Figure14.png', p1, width = 5, height = 6)

}


#VIP Figure 15 Marker VlnPlot by CD8 Cell type
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure15")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    # gene_array <- c("TCF7","ENTPD1","PDCD1","CSF1","CCL3","CCL5","CD80","CD86","HLA-A")
    if (F) {
        markers <- "
    ResidentMemory	CD69
    ResidentMemory	ZNF683
    ResidentMemory	ITGAE
    ResidentMemory	KLRB1
    ResidentMemory	CD160
    ResidentMemory	ITGA1"
    
        markers <- "
			NKT	NCAM1
			NKT	CD3D
			NKT	CD3E
			NKT	KLRD1
			NKT	CD8A
			NKT	CD8B"

        markers <- "
			LR	PDCD1
			LR	CTLA4
			LR	TIGIT
			LR	HAVCR2
			LR	BTLA
			LR	SPP1
			LR	CXCR3
			LR	CXCR6
			LR	CXCR4
			LR	CSF1
			LR	MIF
			LR	CD47"
        markers <- "FABP5"
	}

    marker_frame <- read.table(
        text = gsub("\t\t\t","",markers),
        sep = "\t", 
        col.names = c("Celltype", "Marker")
    )


    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

    cd8_sce$clusters <- factor(cd8_sce$clusters,levels = sort(cd8_celltype_array))

    p <- VlnPlot(cd8_sce, markers,slot = "data")
    ggsave(file = "test.png",plot = p, w = 9,h = 6)
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"), ...) {
        compaired <- list(
            c("CD8-C1", "CD8-C4"),
            c("CD8-C2", "CD8-C4"),
            c("CD8-C3", "CD8-C4"),
            c("CD8-C5", "CD8-C4"),
            c("CD8-C6", "CD8-C4")
        )
        p <- VlnPlot(obj, features = feature, pt.size = pt.size, group.by = "clusters", log = TRUE, ...) +
            xlab("") + ylab(feature) + ggtitle("") +
            theme(
                legend.position = "none",
                # axis.text.x = element_text(size = rel(1), angle = 30),
                axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(),
                axis.title.y = element_text(size = rel(2), angle = 0, vjust = 0.5),
                plot.margin = plot.margin
            ) + geom_signif(
                comparisons = compaired,
                step_increase = 0.3, 
                map_signif_level = T,
                #test.args = "greater", 
                test = wilcox.test
            ) 
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, ...))
        plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
            theme(
                axis.text.x = element_text(size = rel(1.5), angle = 30),
                axis.ticks.x = element_line()
            )
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
        return(p)
    }

    p <- StackedVlnPlot(cd8_sce, markers,
        pt.size = 0, 
        slot = "data",
        cols = celltype_colors
    )

    ggsave(filename = "markers.png", plot = p, width = 8, height = 8, units = c("in"))
    # ggsave(filename = "Figure15-LR.pdf", plot = p, width = 8, height = 18, units = c("in"))


    ggsave(filename = "Figure15-LR.png", plot = p, width = 8, height = 18, units = c("in"))
    ggsave(filename = "Figure15-LR.pdf", plot = p, width = 8, height = 18, units = c("in"))


}


# Figure 16 Marker of CD8 cell type VlnPlot by group
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure16")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)
    markers <- read.table(
        file = "/P02_TNKCell/P999_TNK_Analysis/Figure16/markers.txt",
        header = F, sep = "\t"
    )[, "V1"]

    #orig.ident
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, -1, -0.75, -1), "cm"),...) {
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

    p <- StackedVlnPlot(cd8_sce,
        markers,
        pt.size = 0, 
        cols = group_cols
    )

    # 这里高度必须足够，否则无法展示出VlnPlot
    ggsave(filename = "Figure16.png", plot = p, width = 8, height = 18, units = c("in"))
    ggsave(filename = "Figure16.pdf", plot = p, width = 8, height = 18, units = c("in"))


}


#VIP Figure 17 Marker VlnPlot by CD4-C1 Cell type
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure17")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    gene_array <- c("CD8A","CD8B")
    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD4-C1", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

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
                axis.title.y = element_text(size = rel(2), angle = 0, vjust = 0.5),
                plot.margin = plot.margin
            )
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, ...))
        plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
            theme(
                axis.text.x = element_text(size = rel(1.5), angle = 30),
                axis.ticks.x = element_line()
            )
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
        return(p)
    }

    p <- StackedVlnPlot(cd8_sce, gene_array,
        slot = "scale.data", 
        pt.size = 0.3, 
        cols = celltype_colors
    )

    ggsave(filename = "CD4MarkerVlnPlot.png", plot = p, width = 6, height = 5, units = c("in"))
    ggsave(filename = "CD4MarkerVlnPlot.pdf", plot = p, width = 6, height = 5, units = c("in"))

    cd8_sce = AddMetaData(cd8_sce, cd8_sce[["RNA"]]@scale.data["CD8A", ],col.name = "CD8A")

    cd8_sce = AddMetaData(cd8_sce, cd8_sce[["RNA"]]@scale.data["CCR7", ],col.name = "CCR7")

    meta_frame <- cd8_sce@meta.data %>% filter(CD8A > 1 & CCR7 > 1)

    gp <- ggplot(meta_frame, aes(x = CD8A, y = CCR7)) +
        geom_point(size = 2,shape = 21,fill = "#0B6E78") +
        theme_bw()

    ggsave("CD8CCR7ScatterPlot.png",plot = gp,width = 8,height = 8)
    ggsave("CD8CCR7ScatterPlot.pdf",plot = gp,width = 8,height = 8)

}


#VIP Figure3E CD8 T cell signature VlnPlot by group
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure18")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    tcell_signature_frame <- read.table(file = "/P02_TNKCell/P999_TNK_Analysis/Figure18/CD8_signature.txt", sep = "\t",header = F,stringsAsFactors = F)
    tcell_signature_frame <- tcell_signature_frame %>% filter(!(V1 %in% c("hypoxia")))
    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)
    cd8_sce@meta.data$group_id<- factor(cd8_sce@meta.data$group_id,levels = names(group_cols))
    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(tcell_signature_frame, tcell_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        cd8_sce <- AddModuleScore(
                        object = cd8_sce,
                        features = list(gene_array),
                        name = signature)
		colnames(cd8_sce@meta.data) <- gsub(colnames(cd8_sce@meta.data),
      		pattern = paste0(signature, 1),
			replacement = signature
		)
	}
    
    modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
        compaired <- list(
            c("LRT", "PHT"),
            c("LRN", "PHN"),
            c("PHT", "PHN"),
            c("LRT", "LRN")
        )

        p <- VlnPlot(obj, features = feature, pt.size = pt.size,group.by = 'group_id', log=TRUE,... ) +
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
                # test.args = "less",
                test = wilcox.test
            ) +
                # add boxplot to the violin plot
                geom_boxplot(width = 0.15, fill = "white")
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.25, 0, -0.25, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
                plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
                theme(axis.text.x=element_text(), axis.ticks.x = element_line())
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 4)
        return(p)
    }

    p <- StackedVlnPlot(cd8_sce,
        names(split_signature_list),
        pt.size = 0, 
        cols = group_cols
    )


    ggsave(filename = "CD8SignatureVlnPlotByGroup.png", plot = p, width = 18, height = 6, units = c("in"))
    ggsave(filename = "CD8SignatureVlnPlotByGroup.pdf", plot = p, width = 18, height = 6, units = c("in"))

}


# Figure 19 CD8 T cell DiffusionMap trajectories
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure19")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_sce <- subset(tnk_sce, clusters %in% c(paste0("CD8-C",1:5),"CD4-C1"))

    # Running diffusion map
    cat("Running DiffusionMap...\n")
    cd8_sce_exp <- as.SingleCellExperiment(cd8_sce)
    dm <- DiffusionMap(as.SingleCellExperiment(cd8_sce), verbose = TRUE)
    cat("Complete\n")
    saveRDS(dm, file = "DM.rds", compress = F)


}



# Figure 20 CD8 T Diffusion map
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure20")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- c("CD8-C1","CD8-C2","CD4-C1","CD8-C4","CD8-C5","CD8-C6")
    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

    if (F) {
        sce_array <- c()
        for (cluster in unique(cd8_sce$clusters)) {
            sub_sce <- subset(cd8_sce, clusters %in% c(cluster))
            print(length(Cells(sub_sce)))
            if (length(Cells(sub_sce)) < 1000) {
                cell_id_array <- Cells(sub_sce)
            } else {
                cell_id_array <- sample(Cells(sub_sce), 1000, replace = FALSE)
            }
            sub_sce <- subset(sub_sce, cells = cell_id_array)
            sce_array <- c(sce_array, sub_sce)
        }

        sub_cd8_sce <- sce_array[[1]]
        for (i in 2:length(sce_array)) {
            sub_cd8_sce <- merge(sub_cd8_sce, y = sce_array[[i]])
        }
    }

    # Running diffusion map
    cat("Running DiffusionMap...\n")
    dm <- DiffusionMap(as.SingleCellExperiment(sub_cd8_sce), verbose = TRUE)
    cat("Complete\n")
    saveRDS(dm, file = "DM.rds", compress = F)

    dm <- readRDS("/P02_TNKCell/P999_TNK_Analysis/Figure20/DM.rds")
    # Ploting the diffusion map
    # cd8_sce <- subset(cd8_sce, clusters %in% cd8_celltype_array)
    dc_frame <- data.frame(
        DC1 = eigenvectors(dm)[, 1],
        DC2 = eigenvectors(dm)[, 2],
        cluster = as.character((sub_cd8_sce$clusters),
            group_info = paste0(
                substr(sub_cd8_sce$orig.ident, 1, 2),
                substr(sub_cd8_sce$orig.ident, 5, 6)
            )
        )
    )
    gp <- ggplot(dc_frame, aes(x = DC1, y = DC2, colour = cluster)) +
        geom_point(size = 1,alpha = 0.6,shape = 21) +
        scale_color_manual(values = celltype_colors[cd8_celltype_array]) +
        xlab("Diffusion component 1") +
        ylab("Diffusion component 2") +
        theme_classic() +
        guides(color = guide_legend(override.aes = list(size = 5)))
    ggsave(filename = "DC1_DC2.pdf", plot = gp, width = 8, height = 8, units = c("in"))
    ggsave(filename = "DC1_DC2.png", plot = gp, width = 8, height = 8, units = c("in"))



	# scatter plot by group
	gp <- ggplot(dc_frame, aes(x = DC1, y = DC2, colour = group_info)) +
		geom_point(size = 0.1) +
		xlab("Diffusion component 1") +
		ylab("Diffusion component 2") +
		theme_classic() +
		scale_color_manual(values = c("ZPT" = rgb(45, 67,121, 0, maxColorValue = 255),
									"ZPN" = rgb(45, 67,121, 0, maxColorValue = 255),
									"ZRT" = rgb(45, 67,121, 255, maxColorValue = 255),
									"ZRN" = rgb(45, 67,121, 0, maxColorValue = 255)))

    ggsave(filename = "DC1_DC2_group.pdf", plot = gp, width = 8, height = 8, units = c("in"))
    ggsave(filename = "DC1_DC2_group.png", plot = gp, width = 8, height = 8, units = c("in"))






    # Plotting cell progression along the diffusion map components
    cd8_sce$pseud_dm1 <- rank(eigenvectors(dm)[,1])      # ramyeloid cells by their dpt dm1
    cd8_sce$pseud_dm2 <- rank(eigenvectors(dm)[,2])      # ramyeloid cells by their dpt dm2
    cd8_sce$pseud_dm1R <- rank(-eigenvectors(dm)[,1])    # ramyeloid cells by their dpt dm1 reverse order
    cd8_sce$pseud_dm2R <- rank(-eigenvectors(dm)[,2])    # ramyeloid cells by their dpt dm2 reverse order
    
    group_id_array <- paste0(substr(cd8_sce@meta.data$orig.ident, 1, 2), substr(cd8_sce@meta.data$orig.ident, 5, 5))
    cd8_sce <- AddMetaData(cd8_sce,
        group_id_array,
        col.name = "group_id"
    )

    cd8_sce_exp <- as.SingleCellExperiment(cd8_sce)

    SortedDM1 <- data.frame(
        DM1Sort = as.data.frame(colData(cd8_sce_exp))$pseud_dm1,
        Cluster = as.data.frame(colData(cd8_sce_exp))$clusters,
        Groups = as.data.frame(colData(cd8_sce_exp))$group_id
    )
    SortedDM2 <- data.frame(
        DM2Sort = as.data.frame(colData(cd8_sce_exp))$pseud_dm2,
        Cluster = as.data.frame(colData(cd8_sce_exp))$clusters,
        Groups = as.data.frame(colData(cd8_sce_exp))$group_id
    )
    SortedDM1R <- data.frame(
        DM1SortR = as.data.frame(colData(cd8_sce_exp))$pseud_dm1R,
        Cluster = as.data.frame(colData(cd8_sce_exp))$clusters,
        Groups = as.data.frame(colData(cd8_sce_exp))$group_id
    )
    SortedDM2R <- data.frame(
        DM2SortR = as.data.frame(colData(cd8_sce_exp))$pseud_dm2R,
        Cluster = as.data.frame(colData(cd8_sce_exp))$clusters,
        Groups = as.data.frame(colData(cd8_sce_exp))$group_id
    )
    # By clusters
    gp1 <- ggplot(SortedDM1, aes(x = SortedDM1[, 1], y = Cluster, fill = Cluster)) +
        geom_density_ridges(alpha = 0.5, height = 0.5) +
        geom_boxplot(width = 0.1,alpha = 0.5,outlier.alpha = 0) +
        scale_fill_manual(values = celltype_colors[cd8_celltype_array]) +
        xlab("Diffusion component 1 (DC1)") +
        ylab("Cluster") +
        ggtitle("Cells ordered by DC1") +
        theme_bw()
    ggsave(filename = "SortedDM1.pdf", plot = gp1, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedDM1.png", plot = gp1, width = 8, height = 8, units = c("in"))

    gp2 <- ggplot(SortedDM2, aes(x = SortedDM2[, 1], y = Cluster, fill = Cluster)) +
        geom_density_ridges(alpha = 0.5, height = 0.5) +
        geom_boxplot(width = 0.1,alpha = 0.5,outlier.alpha = 0) +
        scale_fill_manual(values = celltype_colors[cd8_celltype_array]) +
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

    gp1 <- ggplot(SortedDM1, aes(x = SortedDM1[, 1], y = Groups, fill = Groups)) +
        geom_density_ridges(alpha = 0.5, height = 0.5) +
        geom_boxplot(width = 0.1,alpha = 0.5,outlier.alpha = 0) +
        scale_fill_manual(values = group_cols) +
        xlab("Diffusion component 1 (DC1)") +
        ylab("Groups") +
        ggtitle("Cells ordered by DC1") +
        theme_bw()
    ggsave(filename = "SortedGroupsDM1.pdf", plot = gp1, width = 8, height = 4, units = c("in"))
    ggsave(filename = "SortedGroupsDM1.png", plot = gp1, width = 8, height = 4, units = c("in"))

    gp2 <- ggplot(SortedDM2, aes(x = SortedDM2[, 1], y = Groups, fill = Groups)) +
        geom_density_ridges(alpha = 0.5, height = 0.5) +
        geom_boxplot(width = 0.1,alpha = 0.5,outlier.alpha = 0) +
        scale_fill_manual(values = group_cols) +
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
 

}


# Figure21: Find All Markers in TNK types
if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure21")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    marker_frame <- FindAllMarkers(tnk_sce,group.by = "clusters", only.pos = TRUE)
    AllMakers <- 'tnk_all_markers.csv'
    marker_frame <- marker_frame %>% group_by(cluster)
    write.csv(marker_frame, file=AllMakers, quote=F)

    marker_frame <- marker_frame %>%
        group_by(cluster) %>%
        top_n(3, avg_log2FC)
    marker_frame <- marker_frame[!duplicated(marker_frame$gene),]
    #DoHeatmap(tnk_sce, features = marker_frame$gene, group.by = "clusters", label = TRUE)
    dp <- DotPlot(tnk_sce,
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
        filename = "Figure21.png",
        plot = dp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "Figure21.pdf", 
       plot = dp, 
       width = 8, 
       height = 12, 
       units = c("in"))
}


# Figure 22: Marker of NK cell type VlnPlot
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure22")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    gene_array <- c("IL7R","IFNG","FCGR3A","CD69","CX3CR1","CD160","HSPA1A","KLRC1","STMN1","ITGA1")
    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("NK", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

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
                axis.title.y = element_text(size = rel(2), angle = 0, vjust = 0.5),
                plot.margin = plot.margin
            )
        return(p)
    }

    StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
        plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj, feature = x, ...))
        plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
            theme(
                axis.text.x = element_text(size = rel(1.5), angle = 30),
                axis.ticks.x = element_line()
            )
        p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
        return(p)
    }

    p <- StackedVlnPlot(cd8_sce, gene_array,
        slot = "scale.data", 
        pt.size = 0.3, 
        cols = celltype_colors
    )

    ggsave(filename = "Figure22.png", plot = p, width = 6, height = 10, units = c("in"))
    ggsave(filename = "Figure22.pdf", plot = p, width = 6, height = 10, units = c("in"))

}


# Figure23: Find All Markers in NK types and grouping volcano plot
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure23")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("NK", tnk_sce@meta.data$clusters)]))
    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

    marker_frame <- FindAllMarkers(cd8_sce,group.by = "clusters", only.pos = TRUE)
    AllMakers <- 'nk_all_markers.csv'
    marker_frame <- marker_frame %>% group_by(cluster)
    write.csv(marker_frame, file=AllMakers, quote=F)

    marker_frame <- marker_frame %>%
        group_by(cluster) %>%
        top_n(3, avg_log2FC)
    marker_frame <- marker_frame[!duplicated(marker_frame$gene),]
    #DoHeatmap(tnk_sce, features = marker_frame$gene, group.by = "clusters", label = TRUE)
    dp <- DotPlot(cd8_sce,
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
        filename = "Figure23.png",
        plot = dp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "Figure23.pdf", 
       plot = dp, 
       width = 8, 
       height = 12, 
       units = c("in"))


    # Grouping volcano plot
    marker_frame <- read.csv(
        file = "/P02_TNKCell/P999_TNK_Analysis/Figure23/nk_all_markers.csv",
        header = T,
        row.names = 1
    )

    padj <- 0.01
    # 将frame分成两类，一类Up regulate,一类Not significant
    marker_frame$change <- ifelse(marker_frame$avg_log2FC >= 2 & marker_frame$p_val_adj < padj, 
                        "Up regulate","Not significant")
    # 筛选显著marker frame作为标签显示
    label_data <- marker_frame[marker_frame$change != "Not significant", ]
    # 可以根据FC(Fold Change)进行排序，取前50个基因进行展示
    # label_data <- label_data[order(label_data$avg_log2FC, decreasing = T)[c(1:50)], ]
    label_data <- label_data[order(label_data$avg_log2FC, decreasing = T),]

    # 确定显示的基因赋值给label字段
    label_data$label <- label_data$gene

    # 将不显著的基因筛选出，赋值给另一个geom_point
    nolabel_data <- marker_frame[marker_frame$change == "Not significant", ]


    pos <- position_jitter(width = 0.4,seed = 2)

    gp <- ggplot(marker_frame) +
        # 显示标签，设置两个geom_point，一个显示不显著的point，另一个显示显著的point，并且这两个都根据需求设置position，这个是用来和repel协同调整标签位置
        geom_point(data = nolabel_data, aes(cluster, avg_log2FC, color = change),
                    size = 0.85, 
                    #width = 0.4, 
                    alpha = 0.8, 
                    position = pos) +
        geom_point(data = label_data, aes(cluster, avg_log2FC, color = change),
                    size = 0.85, 
                    #width = 0.4, 
                    alpha = 0.8, 
                    position = pos) +
        # 分组方块：
        geom_tile(aes(cluster, 0, fill = cluster),
                    height=0.8,
                    color = "black",
                    alpha = 0.5,
                    show.legend = F,
                    width=0.85) +
        # 文字：
        geom_text(data = marker_frame[!duplicated(marker_frame$cluster), ],
                    aes(cluster, 0, label = cluster),
                    size =2,
                    color ="black") + 
        # 基因标签: 注意这里的postion参数，可以用来调整标签和点的一致性位置
        geom_label_repel(
            data = label_data,
            aes(cluster, avg_log2FC, label = gene),
            size = 2, max.overlaps = 100, box.padding = unit(0.35, "lines"),
            point.padding = unit(0.5, "lines"), 
            segment.colour = "grey50",position = pos
        ) +
        xlab("Cell Subtype")+
        ylab("log2FoldChange")+
        # 颜色模式
        scale_fill_manual(values = brewer.pal(8, "Set3"))+
        scale_color_manual(name = "Regulate", values = c("Up regulate" = "red","Not significant" = "grey"))+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
                legend.position = "top")

    ggsave(
        filename = "GroupingVolcano.png",
        plot = gp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "GroupingVolcano.pdf", 
       plot = gp, 
       width = 8, 
       height = 12, 
       units = c("in"))


}


# Figure 24: NK Diffusion map
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure24")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    nk_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("NK", tnk_sce@meta.data$clusters)]))
    nk_sce <- subset(tnk_sce, clusters %in% nk_celltype_array)

    # Running diffusion map
    cat("Running DiffusionMap...\n")
    nk_sce_exp <- as.SingleCellExperiment(nk_sce)
    dm <- DiffusionMap(as.SingleCellExperiment(nk_sce), verbose = TRUE)
    cat("Complete\n")
    saveRDS(dm, file = "DM.rds", compress = F)

    dm <- readRDS("/P02_TNKCell/P999_TNK_Analysis/Figure20/DM.rds")
    # Ploting the diffusion map
    dc_frame <- data.frame(DC1 = eigenvectors(dm)[, 1],
                    DC2 = eigenvectors(dm)[, 2],
                    DC3 = eigenvectors(dm)[, 3],
                    DC4 = eigenvectors(dm)[, 4],
                    cluster = nk_sce$clusters)
    gp <- ggplot(dc_frame, aes(x = DC1, y = DC2, colour = cluster)) +
        geom_point(size = 0.5)  + 
        scale_color_manual(values = celltype_colors[nk_celltype_array]) +
        xlab("Diffusion component 1") + 
        ylab("Diffusion component 2") +
        theme_classic()
    ggsave(filename = "DC1_DC2.pdf", plot = gp, width = 8, height = 8, units = c("in"))
    ggsave(filename = "DC1_DC2.png", plot = gp, width = 8, height = 8, units = c("in"))

    # Plotting cell progression along the diffusion map components
    nk_sce$pseud_dm1 <- rank(eigenvectors(dm)[,1])      # rank cells by their dpt dm1
    nk_sce$pseud_dm2 <- rank(eigenvectors(dm)[,2])      # rank cells by their dpt dm2
    nk_sce$pseud_dm1R <- rank(-eigenvectors(dm)[,1])    # rank cells by their dpt dm1 reverse order
    nk_sce$pseud_dm2R <- rank(-eigenvectors(dm)[,2])    # rank cells by their dpt dm2 reverse order

    SortedDM1 <- data.frame(DM1Sort = as.data.frame(colData(nk_sce_exp))$pseud_dm1,
                            Samples = as.data.frame(colData(nk_sce_exp))$ident)
    SortedDM2 <- data.frame(DM2Sort = as.data.frame(colData(nk_sce_exp))$pseud_dm2,
                            Samples = as.data.frame(colData(nk_sce_exp))$ident)
    SortedDM1R <- data.frame(DM1SortR = as.data.frame(colData(nk_sce_exp))$pseud_dm1R,
                            Samples = as.data.frame(colData(nk_sce_exp))$ident)
    SortedDM2R <- data.frame(DM2SortR = as.data.frame(colData(nk_sce_exp))$pseud_dm2R,
                            Samples = as.data.frame(colData(nk_sce_exp))$ident)

    gp1 <- ggplot(SortedDM1, aes(x = SortedDM1[, 1], y = Samples, color = Samples)) +
        geom_boxplot() +
        scale_color_manual(values = celltype_colors[nk_celltype_array]) +
        xlab("Diffusion component 1 (DC1)") +
        ylab("Samples") +
        ggtitle("Cells ordered by DC1") +
        theme_bw()
    ggsave(filename = "SortedDM1.pdf", plot = gp1, width = 8, height = 8, units = c("in"))
    ggsave(filename = "SortedDM1.png", plot = gp1, width = 8, height = 8, units = c("in"))

    gp2 <- ggplot(SortedDM2, aes(x = SortedDM2[, 1], y = Samples, color = Samples)) +
        geom_boxplot() +
        scale_color_manual(values = celltype_colors[nk_celltype_array]) +
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


# Figure 26 CD8 Vocalno Plot by Tex vs Tother
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure26")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)
    
    allCells = colnames(cd8_sce)
    allType = names(table(cd8_sce@meta.data$orig.ident))

    if (F) {
        choose_Cells = unlist(lapply(allType, function(x){
            cgCells = allCells[cd8_sce@meta.data$orig.ident == x]
            num = ceiling(dim(as.data.frame(cgCells))[1]*0.3)
            cg = sample(cgCells,num)
            cg  
        }))
    }
    cd8_sce = cd8_sce[, allCells %in% choose_Cells]

    markers <- FindMarkers(cd8_sce,
        ident.1 = "CD8-C4", 
        ident.2 = paste0("CD8-C",1:6)[-4],
        group.by = "clusters", 
        assay = "RNA", 
        slot = "data", 
        logfc.threshold = 0, 
        min.pct = 0
    )

    # GSEA analysis
    mdb_h <- msigdbr(species = "Homo sapiens", category = "C7")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)
    markers$genes = rownames(markers)
    cluster0.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
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



# Figure 27 CD8 T cell signature correlated with marker expression
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure27")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8-C4", tnk_sce@meta.data$clusters)]))

    exhausted_markers <- "
		Exhaustion	FCRL3
		Exhaustion	CD27
		Exhaustion	PRKCH
		Exhaustion	B2M
		Exhaustion	ITM2A
		Exhaustion	TIGIT
		Exhaustion	ID3
		Exhaustion	GBP2
		Exhaustion	PDCD1
		Exhaustion	KLRK1
		Exhaustion	HSPA1A
		Exhaustion	SRGN
		Exhaustion	TNFRSF9
		Exhaustion	TMBIM6
		Exhaustion	TNFRSF1B
		Exhaustion	CADM1
		Exhaustion	ACTB
		Exhaustion	CD8A
		Exhaustion	RGS2
		Exhaustion	FAIM3
		Exhaustion	EID1
		Exhaustion	HSPB1
		Exhaustion	RNF19A
		Exhaustion	IFI16
		Exhaustion	LYST
		Exhaustion	PRF1
		Exhaustion	STAT1
		Exhaustion	UBC
		Exhaustion	CD74
		Exhaustion	IL2RG
		Exhaustion	FYN
		Exhaustion	PTPN6
		Exhaustion	HLA-DRB1
		Exhaustion	HNRNPC
		Exhaustion	UBB
		Exhaustion	CD8B
		Exhaustion	HAVCR2
		Exhaustion	IRF8
		Exhaustion	LAG3
		Exhaustion	ATP5B
		Exhaustion	STAT3
		Exhaustion	IGFLR1
		Exhaustion	MGEA5
		Exhaustion	HSPA1B
		Exhaustion	COTL1
		Exhaustion	VCAM1
		Exhaustion	HLA-DMA
		Exhaustion	PDE7B
		Exhaustion	TBC1D4
		Exhaustion	SNAP47
		Exhaustion	RGS4
		Exhaustion	CBLB
		Exhaustion	TOX
		Exhaustion	CALM2
		Exhaustion	ATHL1
		Exhaustion	SPDYE5
		Exhaustion	DDX5
		Exhaustion	SLA
		Exhaustion	PTPRCAP
		Exhaustion	IRF9
		Exhaustion	MATR3
		Exhaustion	LITAF
		Exhaustion	TPI1
		Exhaustion	ETV1
		Exhaustion	PAM
		Exhaustion	ARID4B
		Exhaustion	NAB1
		Exhaustion	RAPGEF6
		Exhaustion	LDHA
		Exhaustion	WARS
		Exhaustion	RASSF5
		Exhaustion	OSBPL3
		Exhaustion	FAM3C
		Exhaustion	TAP1
		Exhaustion	HLA-DRB6
		Exhaustion	FABP5
		Exhaustion	CD200
		Exhaustion	CTLA4
		Exhaustion	SNX9
		Exhaustion	ETNK1
		Exhaustion	MALAT1
		Exhaustion	ZDHHC6
		Exhaustion	ARL6IP5
		Exhaustion	DUSP2
		Exhaustion	HLA-DQB1
		Exhaustion	HNRNPK
		Exhaustion	DGKH
		Exhaustion	LRMP
		Exhaustion	H3F3B
		Exhaustion	IDH2
		Exhaustion	TRAF5
		Exhaustion	TBL1XR1
		Exhaustion	ANKRD10
		Exhaustion	ALDOA
		Exhaustion	LSP1
		Exhaustion	PTPN7
		Exhaustion	NSUN2
		Exhaustion	RNF149
		Exhaustion	CD2
		Exhaustion	SRSF1
		Exhaustion	GOLPH3
		Exhaustion	HLA-A
		Exhaustion	LIMS1
		Exhaustion	SDF4
		Exhaustion	ROCK1
		Exhaustion	EDEM1
		Exhaustion	APLP2
		Exhaustion	ITK
		Exhaustion	TRIM22
		Exhaustion	SPRY2
		Exhaustion	ACTG1
		Exhaustion	HLA-DPA1
		Exhaustion	EWSR1
		Exhaustion	SRSF4
		Exhaustion	ESYT1
		Exhaustion	LUC7L3
		Exhaustion	ARNT
		Exhaustion	GNAS
		Exhaustion	ARF6
		Exhaustion	ARPC5L
		Exhaustion	NCOA3
		Exhaustion	PAPOLA
		Exhaustion	GFOD1
		Exhaustion	GPR174
		Exhaustion	DDX3X
		Exhaustion	CAPRIN1
		Exhaustion	ARPC2
		Exhaustion	PDIA6
		Exhaustion	SEMA4A
		Exhaustion	CSDE1
		Exhaustion	PSMB9
		Exhaustion	NFATC1
		Exhaustion	PTPN11
		Exhaustion	AGFG1
		Exhaustion	PCED1B
		Exhaustion	CCL4L1
		Exhaustion	CCND2
		Exhaustion	CCL4L2
		Exhaustion	CXCR6
		Exhaustion	AKAP5
		Exhaustion	IFNG
		Exhaustion	MIR155HG
		Exhaustion	ENTPD1
		Exhaustion	TOX2
		Exhaustion	CD7
		Exhaustion	RAB27A
		Exhaustion	ITGAE
		Exhaustion	PHLDA1
		Exhaustion	PAG1
		Exhaustion	CSF1
		Exhaustion	NBL1
		Exhaustion	CCL3
		Exhaustion	ILRB
		Exhaustion	FASLG
		Exhaustion	ZC3H12C
		Exhaustion	MYO7A
		Exhaustion	SIRPG
		Exhaustion	GALNT1
		Exhaustion	UBE2F
		Exhaustion	DUSP4
		Exhaustion	SYT11
		Exhaustion	TRAC
		Exhaustion	TNS3
		Exhaustion	RDH10
		Exhaustion	PTMS
		Exhaustion	CXCL13
		Exhaustion	KIR2DL4"
    exhausted_frame <- read.table(
		text = gsub("\t\t","",exhausted_markers),
		sep = "\t", 
		col.names = c("V1", "V2")
	)
    # cd8_sce <- subset(tnk_sce, clusters %in% c("CD8-C1","CD8-C4"))
    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }

    split_signature_list = split(exhausted_frame, exhausted_frame$V1)
    
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        cd8_sce <- AddModuleScore(
                        object = cd8_sce,
                        features = list(gene_array),
                        name = signature)
    }

    cd8_sce = AddMetaData(cd8_sce, cd8_sce[["RNA"]]@scale.data["FABP5", ],col.name = "FABP5")
    cd8_sce = AddMetaData(cd8_sce, cd8_sce[["RNA"]]@scale.data["PDCD1", ],col.name = "PDCD1")
 	#meta_frame <- cd8_sce@meta.data %>% filter(CSF1 > 0)
 	meta_frame <- cd8_sce@meta.data
    meta_frame <- meta_frame %>% filter((FABP5 > 0) & (PDCD1 > 0))
    stat <- cor.test(meta_frame$PDCD1,meta_frame$FABP5)
    gp <- ggplot(meta_frame, aes(x = PDCD1,y = FABP5)) + 
        geom_point(col = "blue",size = 2,alpha=.05,shape = 19) + 
        geom_smooth(method = lm, colour='#D25565', fill='#fe5f55', size = 1) + 
        theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA)) +
        ggtitle(paste0("R square=",signif(stat$estimate,2)," P value=",signif(stat$p.value,2)))

    ggsave(filename = "FigureFABP5.png", plot = gp, width = 8, height = 10, units = c("in"))
    ggsave(filename = "FigureFABP5.pdf", plot = gp, width = 8, height = 10, units = c("in"))

}



# Figure 28 CD8 T AUCell Signature
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure28")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

    cells_rankings <- AUCell_buildRankings(cd8_sce@assays$RNA@data,nCores = 16)
    c5 <- read.gmt("c5.cc.v7.1.symbols.gmt")

    mdb_h <- msigdbr(species = "Homo sapiens", category = "C7")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    CSF1_geneset <- fgsea_sets[names(fgsea_sets)[grepl("CSF1", names(fgsea_sets))]]
    cells_AUC <- AUCell_calcAUC(CSF1_geneset, 
        cells_rankings,
        aucMaxRank = nrow(cells_rankings) * 0.1
    )

    geneSet <- "GSE11864_UNTREATED_VS_CSF1_IFNG_IN_MAC_DN"
    
    aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
    cd8_sce$AUC <- aucs


    p <- DotPlot(
        object = cd8_sce,
        features = c("AUC"),
        dot.scale = 12,
        cols = c("#0E01FA", "#FB010F"),
        group.by = "clusters",
    ) + theme(axis.text.x = element_text(angle = 90)) +
        theme_bw()

    ggsave("Figure28-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure28-DotPlot.pdf",p,width = 8,height = 8,units = "in")


    plot <- VlnPlot(cd8_sce,
        features = c("AUC"),
        group.by = "clusters",
        pt.size = 0,
        col = celltype_colors,
        ncol = 2
    )
    
    ggsave("Figure28-VlnPlot.png",plot = plot,width = 8,height = 4.5,units = "in")
    ggsave("Figure28-VlnPlot.pdf",plot = plot,width = 8,height = 4.5,units = "in")



    df <- data.frame(cd8_sce@meta.data, cd8_sce@reductions$umap@cell.embeddings)

    class_avg <- df %>%
        group_by(clusters) %>%
        summarise(
            UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2)
        )


    gp <- ggplot(df, aes(UMAP_1, UMAP_2))  +
        geom_point(aes(colour  = AUC)) + viridis::scale_color_viridis(option="A") +
        ggrepel::geom_label_repel(aes(label = clusters),
                                    data = class_avg,
                                    size = 6,
                                    label.size = 0,
                                    segment.color = NA
        )+   theme(legend.position = "none") + theme_bw()

    ggsave(file = "geneSet_UMAP.png", plot = gp, width = 8, height = 8, units = c("in"))
    ggsave(file = "geneSet_UMAP.pdf", plot = gp, width = 8, height = 8, units = c("in"))

}




# Figure 29 CD8 Vocalno Plot by CD8-C2 vs CD8-C3
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure29")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))

    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)
    
    allCells = colnames(cd8_sce)
    allType = names(table(cd8_sce@meta.data$orig.ident))

    if (F) {
        choose_Cells = unlist(lapply(allType, function(x){
            cgCells = allCells[cd8_sce@meta.data$orig.ident == x]
            num = ceiling(dim(as.data.frame(cgCells))[1]*0.3)
            cg = sample(cgCells,num)
            cg  
        }))
        cd8_sce = cd8_sce[, allCells %in% choose_Cells]
    }

    deg <- FindMarkers(cd8_sce,
        ident.1 = "CD8-C2", 
        ident.2 = "CD8-C3",
        group.by = "clusters", 
        assay = "RNA", 
        slot = "data", 
        logfc.threshold = 0, 
        min.pct = 0
    )

    if (F) {
        log2FC_threshold = 0.8
        p_value_threshold = 0.01
        deg[which(deg$p_val_adj < 0.05 & deg$avg_log2FC <= -log2FC_threshold),'sig'] <- 'CD8-C2'
        deg[which(deg$p_val_adj < 0.05 & deg$avg_log2FC >= log2FC_threshold),'sig'] <- 'CD8-C3'
        deg[which(deg$p_val_adj >= 0.05 | abs(deg$avg_log2FC) < log2FC_threshold),'sig'] <- 'None'

        # 横轴 log2FC，纵轴 -log10(adj.P.Val)，颜色表示差异
        deg$Gene <- row.names(deg)
        deg$percent_difference <- deg$pct.1 - deg$pct.2
        p <- ggplot(deg, aes(x = percent_difference, y = avg_log2FC, color = sig)) +
            geom_point(alpha = 0.8, size = 0.6) +
            scale_colour_manual(values = c("#00af3e", "#CD1E24", "gray"), 
                limits = c("CD8-C2", "CD8-C3", "None")) +
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

        up <- subset(deg, sig == 'CD8-C2')
        up <- up[order(up$p_val_adj), ][1:40, ]
        down <- subset(deg, sig == 'CD8-C3')
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

        ggsave('C2_C3_VolcanoPlot.png', p1, width = 6 * 2, height = 4 * 2)
        ggsave('C2_C3_VolcanoPlot.pdf', p1, width = 6 * 2, height = 4 * 2)

    }


    mdb_h <- msigdbr(species = "Homo sapiens", category = "C7")
    fgsea_sets <- mdb_h %>% split(x = .$gene_symbol, f = .$gs_name)

    # GSEA BarPlot
    deg$genes = rownames(deg)
    cluster0.genes <- deg %>%
        arrange(desc(avg_log2FC)) %>%
        dplyr::select(genes, avg_log2FC)
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
            labs(x = "C7", y = "Normalized Enrichment Score", 
            title = "C7 gene sets NES from GSEA") +
            theme_bw()

    ggsave(filename = "C2_C3_gsea_barplot.png", width = 12, height = 6, units = c("in"))
    ggsave(filename = "C2_C3_gsea_barplot.pdf", width = 12, height = 6, units = c("in"))

}





# Figure 30 Summary tissue resident memory T cells

if (F) {
	fig_outdir <- paste0(outdir, "/", "Figure30")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    cd8_cell_stat_frame <- as.data.frame.array(table(
        as.character(cd8_sce$orig.ident),
        as.character(cd8_sce$clusters)
    ))

    write.table(
        x = tnk_cell_stat_frame,
        file = "tnk_immune_stat.xls",
        sep = "\t",
        col.name = TRUE,
        row.name = TRUE,
        quote=FALSE
    )

    cd8_cell_stat_frame$sample <- row.names(cd8_cell_stat_frame)

    cd8_cell_stat_frame$sample_group <- paste0(
        substr(row.names(cd8_cell_stat_frame), 1, 2),
        substr(row.names(cd8_cell_stat_frame), 5, 6)
    )


    bp <- cd8_cell_stat_frame %>%
        gather(celltype, count, -sample) %>%
        mutate(count = as.numeric(count))  %>% 
        ggplot(aes(
                x = sample,
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
                axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 15),
                legend.text = element_text(size = 25),
                axis.text.y = element_text(size = 20),
                axis.title.y = element_text(size = 30)
            )


	ggsave("Figure30.pdf", plot = bp, width = 16, height = 8)
	ggsave("Figure30.png", plot = bp, width = 16, height = 8)


}

# Figure 31 DotPlot for NK
if (F) {
	fig_outdir <- paste0(outdir, "/", "Figure31")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    Idents(tnk_sce) <- "seurat_clusters"


    nk_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("NK", tnk_sce@meta.data$clusters)]))

    nk_sce <- subset(tnk_sce, clusters %in% nk_celltype_array)

	nk_marker <- "
		NK	CD160
		NK	NR4A1
		NK	GZMK
		NK	CXCR6
		NK	GNLY
		NK	FCGR3A"

    cell_marker_frame = read.table(text = gsub("\t\t","",nk_marker),
                sep = "\t", 
                header = F,
    			col.names = c("V1","V2"),
                stringsAsFactors = F)

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$V2),]

    dp <- DotPlot(nk_sce,
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



# Figure 32 NK Signature
if (F) {
	fig_outdir <- paste0(outdir, "/", "Figure32")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

	nk_signature <- "
		Cytotoxic	PRF1
		Cytotoxic	IFNG
		Cytotoxic	GNLY
		Cytotoxic	NKG7
		Cytotoxic	GZMB
		Cytotoxic	GZMA
		Cytotoxic	HZMH
		Cytotoxic	KLRK1
		Cytotoxic	KLRKB1
		Cytotoxic	KLRD1
		Cytotoxic	CTSW
		Cytotoxic	CST7
		Cytotoxic	CX3CR1
		Cytotoxic	FGFBP2
		Cytotoxic	S1PR5
		Cytotoxic	FCGR3A
		Cytotoxic	PLAC8
		Cytotoxic	C1orf21
		Cytotoxic	TGFBR3
		Cytotoxic	PLEK
		Cytotoxic	FGR
		Cytotoxic	KLRF1
		Cytotoxic	SPON2
		Cytotoxic	CD300A
		Cytotoxic	S1PR1
		Cytotoxic	FAM65B
		Cytotoxic	EFHD2
		Cytotoxic	STK38
		Cytotoxic	C1orf162
		Cytotoxic	SORL1
		Cytotoxic	EMP3
		Cytotoxic	ARL4C
		Cytotoxic	BIN2
		Cytotoxic	CCND3
		Cytotoxic	FCRL6
		Cytotoxic	SAMD3
		Cytotoxic	TRDC
		Cytotoxic	TYROBP
		Cytotoxic	GNLY
		Cytotoxic	KLRG1
		Exhaustion	FCRL3
		Exhaustion	CD27
		Exhaustion	PRKCH
		Exhaustion	B2M
		Exhaustion	ITM2A
		Exhaustion	TIGIT
		Exhaustion	ID3
		Exhaustion	GBP2
		Exhaustion	PDCD1
		Exhaustion	KLRK1
		Exhaustion	HSPA1A
		Exhaustion	SRGN
		Exhaustion	TNFRSF9
		Exhaustion	TMBIM6
		Exhaustion	TNFRSF1B
		Exhaustion	CADM1
		Exhaustion	ACTB
		Exhaustion	CD8A
		Exhaustion	RGS2
		Exhaustion	FAIM3
		Exhaustion	EID1
		Exhaustion	HSPB1
		Exhaustion	RNF19A
		Exhaustion	IFI16
		Exhaustion	LYST
		Exhaustion	PRF1
		Exhaustion	STAT1
		Exhaustion	UBC
		Exhaustion	CD74
		Exhaustion	IL2RG
		Exhaustion	FYN
		Exhaustion	PTPN6
		Exhaustion	HLA-DRB1
		Exhaustion	HNRNPC
		Exhaustion	UBB
		Exhaustion	CD8B
		Exhaustion	HAVCR2
		Exhaustion	IRF8
		Exhaustion	LAG3
		Exhaustion	ATP5B
		Exhaustion	STAT3
		Exhaustion	IGFLR1
		Exhaustion	MGEA5
		Exhaustion	HSPA1B
		Exhaustion	COTL1
		Exhaustion	VCAM1
		Exhaustion	HLA-DMA
		Exhaustion	PDE7B
		Exhaustion	TBC1D4
		Exhaustion	SNAP47
		Exhaustion	RGS4
		Exhaustion	CBLB
		Exhaustion	TOX
		Exhaustion	CALM2
		Exhaustion	ATHL1
		Exhaustion	SPDYE5
		Exhaustion	DDX5
		Exhaustion	SLA
		Exhaustion	PTPRCAP
		Exhaustion	IRF9
		Exhaustion	MATR3
		Exhaustion	LITAF
		Exhaustion	TPI1
		Exhaustion	ETV1
		Exhaustion	PAM
		Exhaustion	ARID4B
		Exhaustion	NAB1
		Exhaustion	RAPGEF6
		Exhaustion	LDHA
		Exhaustion	WARS
		Exhaustion	RASSF5
		Exhaustion	OSBPL3
		Exhaustion	FAM3C
		Exhaustion	TAP1
		Exhaustion	HLA-DRB6
		Exhaustion	FABP5
		Exhaustion	CD200
		Exhaustion	CTLA4
		Exhaustion	SNX9
		Exhaustion	ETNK1
		Exhaustion	MALAT1
		Exhaustion	ZDHHC6
		Exhaustion	ARL6IP5
		Exhaustion	DUSP2
		Exhaustion	HLA-DQB1
		Exhaustion	HNRNPK
		Exhaustion	DGKH
		Exhaustion	LRMP
		Exhaustion	H3F3B
		Exhaustion	IDH2
		Exhaustion	TRAF5
		Exhaustion	TBL1XR1
		Exhaustion	ANKRD10
		Exhaustion	ALDOA
		Exhaustion	LSP1
		Exhaustion	PTPN7
		Exhaustion	NSUN2
		Exhaustion	RNF149
		Exhaustion	CD2
		Exhaustion	SRSF1
		Exhaustion	GOLPH3
		Exhaustion	HLA-A
		Exhaustion	LIMS1
		Exhaustion	SDF4
		Exhaustion	ROCK1
		Exhaustion	EDEM1
		Exhaustion	APLP2
		Exhaustion	ITK
		Exhaustion	TRIM22
		Exhaustion	SPRY2
		Exhaustion	ACTG1
		Exhaustion	HLA-DPA1
		Exhaustion	EWSR1
		Exhaustion	SRSF4
		Exhaustion	ESYT1
		Exhaustion	LUC7L3
		Exhaustion	ARNT
		Exhaustion	GNAS
		Exhaustion	ARF6
		Exhaustion	ARPC5L
		Exhaustion	NCOA3
		Exhaustion	PAPOLA
		Exhaustion	GFOD1
		Exhaustion	GPR174
		Exhaustion	DDX3X
		Exhaustion	CAPRIN1
		Exhaustion	ARPC2
		Exhaustion	PDIA6
		Exhaustion	SEMA4A
		Exhaustion	CSDE1
		Exhaustion	PSMB9
		Exhaustion	NFATC1
		Exhaustion	PTPN11
		Exhaustion	AGFG1
		Exhaustion	PCED1B
		Exhaustion	CCL4L1
		Exhaustion	CCND2
		Exhaustion	CCL4L2
		Exhaustion	CXCR6
		Exhaustion	AKAP5
		Exhaustion	IFNG
		Exhaustion	MIR155HG
		Exhaustion	ENTPD1
		Exhaustion	TOX2
		Exhaustion	CD7
		Exhaustion	RAB27A
		Exhaustion	ITGAE
		Exhaustion	PHLDA1
		Exhaustion	PAG1
		Exhaustion	CSF1
		Exhaustion	NBL1
		Exhaustion	CCL3
		Exhaustion	ILRB
		Exhaustion	FASLG
		Exhaustion	ZC3H12C
		Exhaustion	MYO7A
		Exhaustion	SIRPG
		Exhaustion	GALNT1
		Exhaustion	UBE2F
		Exhaustion	DUSP4
		Exhaustion	SYT11
		Exhaustion	TRAC
		Exhaustion	TNS3
		Exhaustion	RDH10
		Exhaustion	PTMS
		Exhaustion	CXCL13
		Exhaustion	KIR2DL4"

		
	nk_signature_frame <- read.table(
		text = gsub("\t\t","",nk_signature),
		sep = "\t", 
		col.names = c("V1", "V2")
	)

	nk_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("NK", tnk_sce@meta.data$clusters)]))

    nk_sce <- subset(tnk_sce, clusters %in% nk_celltype_array)

    addscore <- function(obj, feature, name) {
        sce <- AddModuleScore(
            object = obj,
            features = feature,
            name = name
        )
        return(sce)
    }
    split_signature_list = split(nk_signature_frame, nk_signature_frame$V1)
    for(signature in names(split_signature_list)) {
        gene_array <- split_signature_list[[signature]]$V2
        nk_sce <- AddModuleScore(
            object = nk_sce,
            features = list(gene_array),
            name = signature
        )
        colnames(nk_sce@meta.data) <- gsub(colnames(nk_sce@meta.data),
            pattern = paste0(signature, 1),
			replacement = signature
        )
    }

    p <- DotPlot(
        object = nk_sce,
        features = names(split_signature_list),
        dot.scale = 12,
        cols = c("#0E01FA", "#FB010F"),
        group.by = "clusters",
    ) + theme_bw() + 
	theme(axis.text.x = element_text(size = 15),
			axis.text.y = element_text(size = 15)) +
        xlab("") +
        ylab("")

    ggsave("Figure9-DotPlot.png",p,width = 8,height = 8,units = "in")
    ggsave("Figure9-DotPlot.pdf",p,width = 8,height = 8,units = "in")


    plot <- VlnPlot(nk_sce,
        features = names(split_signature_list),
        group.by = "clusters",
        col = celltype_colors,
        ncol = 2
    )
    
    ggsave("Figure9-VlnPlot.png",plot = plot,width = 8,height = 4.5,units = "in")
    ggsave("Figure9-VlnPlot.pdf",plot = plot,width = 8,height = 4.5,units = "in")

}







# Figure 33 T or NK subpopulation survival analysis
if (F) {

    library(data.table)
    library(magrittr)
    library(biomaRt)
    library(limma)
    library(survival)
    library(survminer)
    library(dplyr)
    library(Seurat)

    fig_outdir <- paste0(outdir, "/", "Figure33")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    # input parameter
    global_seurat_rds = "/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds"
    subcluster_seurat_rds = "/P02_TNKCell/P999_TNK_Analysis/Figure0/TNK.rds"
    celltype = "CD8-C1"

    # load TCGA data: LIHC_expr: LIHC TCGA expression matrix;
    # LIHC_survival: LIHC TCGA survival information, sample/OS/OS.time;

    LIHC_expr <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.expr.rds")
    LIHC_survival <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.survival.rds")

    # load global Seurat object
    global_sce <- readRDS(global_seurat_rds)

    # load subtype Seurat object
    subcluster_sce <- readRDS(subcluster_seurat_rds)

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





# Figure 34 Specify gene set to get survival results
if (F) {


    library(data.table)
    library(magrittr)
    library(biomaRt)
    library(limma)
    library(survival)
    library(survminer)
    library(dplyr)
    library(Seurat)

    fig_outdir <- paste0(outdir, "/", "Figure34")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    gene_array <- c("PMCH",
		"HLA-DRA",
		"DUSP4",
		"CD74",
		"VCAM1",
		"HLA-DRB1",
		"RGS1",
		"HLA-DQB1",
		"HLA-DRB5",
		"RGS2",
		"HLA-DPA1",
		"HLA-DPB1",
		"TNFRSF9",
		"HLA-DQA1",
		"MTRNR2L8",
		"LAG3",
		"GZMK",
		"CRTAM",
		"ENC1",
		"HMOX1",
		"MS4A6A",
		"NKG7",
		"HLA-DMA",
		"SNAP47",
		"FAM3C",
		"FXYD2",
		"CST7",
		"BLVRA",
		"TSC22D1",
		"BHLHE40-AS1",
		"CD8A",
		"PDCD1",
		"MAP1LC3A",
		"SEMA4A",
		"CMC1",
		"HAVCR2",
		"FAM174B",
		"PMAIP1",
		"CD27",
		"CD8B",
		"LYST",
		"ITM2A",
		"PON2",
		"RASSF4",
		"AC069363.1",
		"PTMS",
		"HLA-DQA2",
		"CAV1",
		"TMEM155",
		"UBE2F",
		"CREM",
		"TOX",
		"CCL5",
		"SLA",
		"CD82",
		"PHLDA1",
		"PON3",
		"APOBEC3C",
		"TTN",
		"SYNGR2",
		"CLECL1",
		"RASD1",
		"CCL3",
		"CXCL13",
		"PVR",
		"CADM1",
		"TIGIT",
		"SPP1",
		"LAMP1",
		"ADGRG5",
		"ICAM1",
		"TNFRSF1B",
		"CD44",
		"TFRC",
		"TNFRSF1A",
		"TNFRSF10B",
		"TNF",
		"COPA",
		"ICOSLG",
		"LTA",
		"TNFSF13",
		"CXCR3",
		"CD46",
		"CXCL12",
		"LGALS9",
		"SIRPA",
		"TNFRSF10C",
		"PLXNB2",
		"CD72",
		"MIF",
		"BTLA",
		"TNFSF10",
		"CXCR6",
		"PLD2",
		"CSF1",
		"ADORA2B",
		"TNFSF9",
		"CD2",
		"CSF1R",
		"CD55",
		"CTLA4",
		"CD99",
		"SPN",
		"CD6",
		"LTBR",
		"LILRA4",
		"CD47",
		"FN1",
		"PLAUR",
		"IL10RA",
		"IL10RB",
		"IFNG",
		"CD96",
		"NECTIN2",
		"PTGER4",
		"PDCD1LG2",
		"AREG",
		"GRN",
		"HBEGF",
		"TNFSF13B",
		"ITGAL",
		"FASLG",
		"VSIR",
		"SORT1",
		"CLEC2D",
		"ICOS",
		"FAS",
		"RIPK1",
		"CCL20",
		"JAG1",
		"SLC1A5",
		"SEMA4D",
		"TNFRSF14",
		"TNFRSF10D",
		"CXCL16",
		"CXCR4",
		"ARF1",
		"PLXND1",
		"ENTPD1",
		"CD58",
		"ADGRE5",
		"CD80",
		"CD86",
		"PILRA",
		"SIGLEC1",
		"ALCAM",
		"LTB",
		"BST2",
		"SIRPG",
		"IL10",
		"IL15",
		"ITGA4",
		"ITGB1",
		"ITGB7",
		"ITGAM",
		"ITGB2",
		"ITGAX",
		"IFNGR1",
		"IFNGR2",
		"IL15RA",
		"IL2RB",
		"IL2RG")
    
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
    # surdata$level <- ifelse(surdata[,"mean_expr"] > as.numeric(quantile(surdata[,"mean_expr"],0.75)),'High','Low')
    fit <- survfit(Surv(OS.time,OS)~level,data = surdata)
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
    
    png(file = paste0(celltype,".survival.png"),width = 1024,height = 768)
    print(p)
    dev.off()

    pdf(file = paste0(celltype,".survival.pdf"),width = 12,height = 8)
    print(p)
    dev.off()


}


# Figure 35 CSF gene expression correlated with CD8 Exhausted T cell infiltration
if (F) {
    
    fig_outdir <- paste0(outdir, "/", "Figure35")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    LIHC_expr <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.expr.rds")
    LIHC_transpose_expr <- as.data.frame(t(LIHC_expr))
    LIHC_transpose_expr <- LIHC_transpose_expr %>% select("CSF1","CCL3","CCL4","CCL5")


    cell_marker <- "/mnt/share01/projects/single_cell/analysis/RNA/data/cell_type_markers.csv"


    d = log2(LIHC_expr + 1)

    d = t(d)
    e = d
    ## cell type scores:
    cellgenes = read.csv(cell_marker,stringsAsFactors = F)

    lost.genes = setdiff(cellgenes[,1],colnames(d))
    kept.genes = intersect(cellgenes[,1],colnames(d))
    cellgenes$kept = NA; 
    cellgenes$kept[is.element(cellgenes$Gene,kept.genes)]=1
    cellgenes$kept[is.element(cellgenes$Gene,lost.genes)]=0
    #write.csv(cellgenes,file="cell type markers missing from Chen2016.csv")
    cellgenes = cellgenes[is.element(cellgenes[,1],colnames(d)),]

    cells = unique(cellgenes[,2])

    cellscores = matrix(NA,nrow(e),length(cells)); dimnames(cellscores) = list(rownames(e),cells)

    for(cell in cells)
    {
    cellscores[,cell] = rowMeans(e[,cellgenes[cellgenes[,2]==cell,1],drop=F])
    }


    # calc allTotTils as elsewhere:
    high.cor.w.cd45 = cor(cellscores)[,"CD45"]>0.6
    alltottils = rowMeans(cellscores[,high.cor.w.cd45])
    cellscores = cbind(alltottils,cellscores)

    write.table(cellscores,file = "escc.unscaled_cellscore_matrix.txt",sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)



    cellscores <- scale(cellscores,center = TRUE,scale = TRUE)
    gene_cellscores_frame <- merge(LIHC_transpose_expr,cellscores,by = "row.names")


    stat <- cor.test(gene_cellscores_frame$CCL5, gene_cellscores_frame$"Exhausted CD8")


    gp <- ggplot(gene_cellscores_frame,aes(x = CCL5,y = gene_cellscores_frame$"Exhausted CD8")) + 
        geom_point(col = "blue",size = 2,alpha=.05,shape = 19) + 
        geom_smooth(method = lm, colour='#D25565', fill='#fe5f55', size = 1) + 
        theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA),
        axis.text = element_blank()) +
        ggtitle(paste0("R square=", signif(stat$estimate, 2), " P value=", signif(stat$p.value, 2))) +
            xlab("CSF1 expression") +
            ylab("Exhausted CD8 score")

    ggsave("CSF1_ExhaustedCD8_cor.png",plot = gp,width = 8,height = 8)
    ggsave("CSF1_ExhaustedCD8_cor.pdf",plot = gp,width = 8,height = 8)

}




# Figure 36: Euclidean distance between ZPT&ZPN vs ZRTvsZRN

if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure36")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    cell_stat_frame <- as.data.frame.array(table(
        as.character(tnk_sce$orig.ident),
        as.character(tnk_sce$clusters)
    ))

    immune_celltype_array <- as.character(unique(tnk_sce$clusters))

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
                (substr(sample1, 5, 5) == "T") & (substr(sample2,5,5) == "T")) ~ "tumor_uunspliced_matched",
            ((substr(sample1, 1, 4) != substr(sample2, 1, 4)) &
                (substr(sample1, 5, 5) == "N") & (substr(sample2,5,5) == "N")) ~ "normal_uunspliced_matched",
            ((substr(sample1, 1, 4) != substr(sample2, 1, 4)) &
                ((substr(sample1, 5, 5) == "T") & (substr(sample2,5,5) == "N") |
                (substr(sample1, 5, 5) == "N") & (substr(sample2,5,5) == "T"))) ~ "tumor_normal_uunspliced_matched",
        ))

    compaired <- list(
        c("ZP_tumor_matched", "ZR_tumor_matched"),
        c("tumor_normal_uunspliced_matched","ZP_tumor_matched"),
        c("tumor_normal_uunspliced_matched","ZR_tumor_matched")
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

    ggsave(filename = "PHTvsLRT_TNK_DistanceBoxPlot.png",plot = gp,width = 8,height = 6,units = c("in"))
    ggsave(filename = "PHTvsLRT_TNK_DistanceBoxPlot.pdf",plot = gp,width = 8,height = 6,units = c("in"))

}



#VIP Figure37: Find All Markers in CD8 types: DotPlot and Grouping volcano plot
if (F) {

    fig_outdir <- paste0(outdir, "/", "Figure37")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8", tnk_sce@meta.data$clusters)]))
    cd8_sce <- subset(tnk_sce, clusters %in% cd8_celltype_array)

    marker_frame <- FindAllMarkers(cd8_sce,group.by = "clusters", only.pos = TRUE)
    AllMakers <- 'cd8_all_markers.csv'
    marker_frame <- marker_frame %>% group_by(cluster)
    write.csv(marker_frame, file=AllMakers, quote=F)
    marker_frame <- read.csv(file = "/P02_TNKCell/P999_TNK_Analysis/Figure37/cd8_all_markers.csv", header = T,row.names = 1)
    marker_frame <- marker_frame %>%
        group_by(cluster) %>%
        top_n(15, avg_log2FC)
    # marker_frame <- marker_frame[!duplicated(marker_frame$gene),]
    #DoHeatmap(tnk_sce, features = marker_frame$gene, group.by = "clusters", label = TRUE)
    dp <- DotPlot(cd8_sce,
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
        filename = "CD8FindAllMarkersDotPlot.png",
        plot = dp,
        width = 8,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "CD8FindAllMarkersDotPlot.pdf", 
       plot = dp, 
       width = 8, 
       height = 12, 
       units = c("in"))

    # Grouping volcano plot
    # 将frame分成两类，一类Up regulate,一类Not significant

    label_frame <- marker_frame %>%
        group_by(cluster) %>%
        top_n(15, avg_log2FC)

    # 使用anti_join获取剩下的数据
    nolabel_frame <- anti_join(marker_frame,label_frame,by = c("cluster","gene"))

    label_frame <- label_frame[order(label_frame$avg_log2FC, decreasing = T),]
    label_frame$change <- "Up Regulated"
    # 确定显示的基因赋值给label字段

    nolabel_frame$change <- "Not Significant"

    pos <- position_jitter(width = 0.4,seed = 2)

    gp <- ggplot(marker_frame) +
        # 显示标签，设置两个geom_point，一个显示不显著的point，另一个显示显著的point，并且这两个都根据需求设置position，这个是用来和repel协同调整标签位置
        geom_point(data = nolabel_frame, aes(cluster, avg_log2FC, color = change),
                    size = 0.85, 
                    #width = 0.4, 
                    alpha = 0.8, 
                    position = pos) +
        geom_point(data = label_frame, aes(cluster, avg_log2FC, color = change),
                    size = 0.85, 
                    #width = 0.4, 
                    alpha = 0.8, 
                    position = pos) +
        # 分组方块：
        geom_tile(aes(cluster, 0, fill = cluster),
                    height=0.8,
                    color = "black",
                    alpha = 0.5,
                    show.legend = F,
                    width=0.85) +
        # 文字：
        geom_text(data = marker_frame[!duplicated(marker_frame$cluster), ],
                    aes(cluster, 0, label = cluster),
                    size = 8,
                    color ="white") + 
        # 基因标签: 注意这里的postion参数，可以用来调整标签和点的一致性位置
        geom_label_repel(
            data = label_frame,
            aes(cluster, avg_log2FC, label = gene),
            size = 4, max.overlaps = 100, box.padding = unit(0.35, "lines"),
            point.padding = unit(0.5, "lines"), 
            segment.colour = "grey50",position = pos
        ) +
        xlab("")+
        ylab("log2FoldChange")+
        # 颜色模式
        scale_fill_manual(values = celltype_colors)+
        scale_color_manual(name = "Regulate", values = c("Up Regulated" = "red",
            "Not Significant" = "grey"))+
        theme_bw()+
        theme(axis.text.x = element_blank(),
                axis.text.y = element_text(size = 15),
                axis.title.y = element_text(size = 20),
                legend.position = "None")

    ggsave(
        filename = "CD8GroupingVolcano.png",
        plot = gp,
        width = 16,
        height = 12,
        units = c("in")
    )

    ggsave(filename = "CD8GroupingVolcano.pdf", 
       plot = gp, 
       width = 16, 
       height = 12, 
       units = c("in"))


}




# Figure38: demo velocyto analysis
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure38")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    # radian --r-binary=/mnt/share01/tools/miniconda/envs/velocyto/bin/R
    library(velocyto.R)
    library(SeuratWrappers)
    library(stringr)
    library(ggplot2)
    library(Seurat)

    ldat <- ReadVelocity(file = "SCG71.loom")
    bm <- as.Seurat(x = ldat)

    bm <- SCTransform(object = bm, assay = "spliced") %>%
        RunPCA(verbose = FALSE) %>%
        FindNeighbors(dims = 1:20) %>%
        FindClusters() %>%
        RunUMAP(dims = 1:20)

    dp <- DimPlot(bm)
    ggsave("demo.png",plot = dp,width = 8, height = 8)


    bm <- RunVelocity(
        object =bm,
        deltaT = 1,
        kcells = 25,
        fit.quantile = 0.02) 

    ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
    names(x = ident.colors) <- levels(x = bm)
    cell.colors <- ident.colors[Idents(object = bm)]
    names(x = cell.colors) <- colnames(x = bm)

    png("demo.png",width = 1080,height = 1080)
    show.velocity.on.embedding.cor(
        emb = Embeddings(object = bm, reduction = "umap"),
        vel = Tool(object = bm, slot = "RunVelocity"),
        n = 200,
        scale = "sqrt",
        cell.colors = ac(x = cell.colors, alpha = 0.5),
        cex = 0.8,
        arrow.scale = 3,
        show.grid.flow = TRUE,
        min.grid.flow = TRUE,
        min.grid.cell.mass = 0.5,
        grid.n = 40,
        arrow.lwd = 1,
        do.par = FALSE,
        cell.border.alpha = 0.1
    )

    dev.off()

}




# Figure39: CD8-C5, CD8-C1, CD8-C4 velocyto analysis
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure39")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    # get loom python script
    if (F) {
        import sys
        import os
        import logging
        import pandas as pd
        import numpy as np
        import argparse
        from argparse import RawTextHelpFormatter

        if __name__ == '__main__':
            parser = argparse.ArgumentParser(description = "",formatter_class=RawTextHelpFormatter)
            parser.add_argument('--infile',help = "celltype file",required = True)
            parser.add_argument('--outdir',help = "",required = True)
            #parser.add_argument('--prefix',help = "",required = True)
            argv = parser.parse_args()
            samples = argv.infile    
            outdir = argv.outdir

            with open(samples) as f:
                for line in f:
                    sample, path = line.strip().split('\t')
                    if sample.startswith("#"):
                        continue
                    os.system(('/mnt/share01/tools/miniconda/envs/velocyto/bin/velocyto run -U \
                                -b {0}/filtered_feature_bc_matrix/barcodes.tsv.gz -o {1}/{2} \
                                -m /seurat_xiefei/Analysis_34_samples/velocity/files/human_hg19_rmsk.gtf \
                                {0}/possorted_genome_bam.bam /seurat_xiefei/Analysis_34_samples/velocity/files/human_hg19.gtf').format(path,outdir,sample))
    }


    # radian --r-binary=/mnt/share01/tools/miniconda/envs/velocyto/bin/R
    library(velocyto.R)
    library(SeuratWrappers)
    library(stringr)
    library(ggplot2)
    library(Seurat)
    library(dior)

    #tnk_sce <- readRDS(tnk_seurat_rds)
    cd8_h5 <- "/P02_TNKCell/P999_TNK_Analysis/Figure41/cd8_diffusion.h5"
    cd8_sce <- dior::read_h5(cd8_h5)


    data = cd8_sce@reductions$diffmap_@cell.embeddings %>%
        as.data.frame() %>%
        cbind(clusters = tnk_sce@meta.data$clusters)

    # step2 获取要添加标签的位置
    class_avg <- data %>% 
        dplyr::group_by(clusters) %>% 
        dplyr::summarise(
            DIFFMAP_1 = median(DIFFMAP_1),
            DIFFMAP_2 = median(DIFFMAP_2)
        )

    # step2 绘图
    umap <- ggplot(data ,aes(x=DIFFMAP_1,y=DIFFMAP_2))+
        #geom_point(aes(color = clusters),size = 0.2, alpha = 0.2)+ 
        geom_point(aes(color = clusters),size = 0.2)+ 
        scale_color_manual(values = celltype_colors)+
        #geom_text(aes(label = clusters), data = class_avg, size = 5)+
        xlab('DC dimension 1') +
        ylab('DC dimension 2') +
        theme_bw() +
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

    ggsave("CD8_DiffusionMap.pdf",plot = umap,width = 12,height = 8)
    ggsave("CD8_DiffusionMap.png",plot = umap,width = 12,height = 8)






    #sub_tnk_sce <- subset(tnk_sce,group %in% c("PHT","PHN")) # nolint
    #seurat.object <- tnk_sce[,tnk_sce@meta.data$clusters %in% c('CD8-C5','CD8-C1','CD8-C4','CD8-C2','CD8-C3','CD8-C6')] # nolint
    seurat.object <- tnk_sce
    #ldat <- ReadVelocity(file = '/loom/velocity/merged.loom')
    #saveRDS(ldat, file = "merged.loom.rds", compress = F)
    
    ldat <- readRDS("/P02_TNKCell/P999_TNK_Analysis/Figure39/merged.loom.rds")
    # 获取spliced矩阵和unspliced_mat矩阵
    spliced_mat <- ldat$spliced
    unspliced_mat <- ldat$unspliced

    # 统一loom对象和Seurat的细胞名与基因名
    colnames(spliced_mat) <- paste0(colnames(spliced_mat), "-1")
    colnames(spliced_mat) <- gsub("x","",colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_5W6L2:", "ZP01T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_A5XED:", "ZP01N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_P12RO:", "ZP02T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_AJV1V:", "ZP02N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_QQYJ1:", "ZP03T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_PUH87:", "ZP03N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_PFCBY:", "ZP04T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_Q7WF1:", "ZP04N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_JTTWU:", "ZP05T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_UPJNB:", "ZP05N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_JWGIC:", "ZP06T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_X5LH5:", "ZP06N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_D3B9J:", "ZP07T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_7RKKV:", "ZP07N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_GT5XR:", "ZP08T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_RDA3X:", "ZP08N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_5K19I:", "ZP09T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_PQYBV:", "ZP09N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_M8KZA:", "ZP10T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_SCNS4:", "ZP10N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_H3L3G:", "ZP11T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_Z80IH:", "ZP11N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_O9NOQ:", "ZR01T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_O98YO:", "ZR01N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_NGKZG:", "ZR02T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_E3G36:", "ZR02N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_9D780:", "ZR03T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_TGQ3W:", "ZR03N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_R1UDI:", "ZR04T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_TC8IF:", "ZR04N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_65V7P:", "ZR05T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_OH191:", "ZR05N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_M6INY:", "ZR06T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_QH8IO:", "ZR06N_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_BPFCD:", "ZR07T_", colnames(spliced_mat))
    colnames(spliced_mat) <- gsub("possorted_genome_bam_35QSI:", "ZR07N_", colnames(spliced_mat))

    # 统一splice和unspliced对象的列名
    colnames(unspliced_mat) <- colnames(spliced_mat)

    # loom的cell id与Seurat的cell id取交集(有时loom不会包含所有的seurat的cell id)
    # 这里需要注意查看intersect之后获得的cell id的数量
    cell_id_array <- intersect(rownames(seurat.object@meta.data), colnames(spliced_mat))
    spliced_mat <- spliced_mat[, cell_id_array]
    unspliced_mat <- unspliced_mat[, cell_id_array]
    sub_sce <- subset(seurat.object,cells = cell_id_array)

    cat("total cell:", length(Cells(sub_sce)))

    sub_sce[["spliced"]] <- CreateAssayObject(counts = as.matrix(spliced_mat))
    sub_sce[["unspliced"]] <- CreateAssayObject(counts = as.matrix(unspliced_mat))

    # cell type as character
    sub_sce$clusters <- as.character(sub_sce$clusters)

    sub_sce <- subset(sub_sce, cells = Cells(cd8_sce))
   
    #kCells：用于斜率平滑度计算最近邻细胞数量，越大越平滑，越小越能反映局部特征
    #fit.quantile：0.02代表对基因表达量最高2%与最低2%的值执行gamma拟合
    #spliced.average：过滤低表达丰度基因的标准，计算的是基因在cluster内的平均counts值
    #unspliced.average：同上

    velo <- RunVelocity(sub_sce, deltaT = 1, kCells = 25, fit.quantile = 0.02, 
        spliced.average = 0.2, unspliced.average = 0.05, ncores = 32)

    emb = Embeddings(cd8_sce, reduction = "diffmap_")
    vel = Tool(velo, slot = "RunVelocity")

    pdf("velocity.umap.pdf",w=12,h=8)
    show.velocity.on.embedding.cor(
        emb = emb,
        vel = vel,
        n = 200,
        scale = "sqrt",
        cell.colors = ac(colors, alpha = 0.5),
        cex = 0.8,
        arrow.scale = 3,
        show.grid.flow = FALSE,
        min.grid.cell.mass = 0.5,
        grid.n = 40,
        arrow.lwd = 1,
        do.par = FALSE, 
        cell.border.alpha = 0.1,
        n.cores = 32
    )
    dev.off()



    # 提取原有的用seurat做的UMAP的图
    emb <- sub_sce@reductions$umap@cell.embeddings
    # Estimate the cell-cell distances 
    cell.dist <- as.dist(1-armaCor(t(sub_sce@reductions$umap@cell.embeddings)))

    # velocity 计算 
    fit.quantile <- 0.02
    # Main velocity estimation
    rvel.cd <- gene.relative.velocity.estimates(spliced_mat,
        unspliced_mat,
        deltaT = 2,
        kCells = 10,
        cell.dist = cell.dist,
        fit.quantile = fit.quantile,
        n.cores = 12
    )

    saveRDS(rvel.cd,file = "rvel_cd.rds",compress = F)

    # 提取原seurat对象的UMAP图的细胞颜色
    gg <- UMAPPlot(sub_sce)
    ggplot_build(gg)$data
    colors <- as.list(ggplot_build(gg)$data[[1]]$colour)
    names(colors) <- rownames(emb)

    # 将velocity画在原有的seurat的UMAP上
    pdf("velocity.umap.pdf",w=12,h=8)
    show.velocity.on.embedding.cor(emb,
                rvel.cd,
                n = 300, # n = 200箭头太密
                scale = 'sqrt',
                cell.colors=ac(colors,alpha = 0.8),
                cex = 0.3,
                arrow.scale = 3,
                show.grid.flow = T,
                min.grid.flow = TRUE,
                min.grid.cell.mass = 1,
                grid.n = 50,
                arrow.lwd = 1.5,
                do.par=F,
                cell.border.alpha = 0.1,
                n.cores = 12,
                main="T Cell Velocity"
    )
    dev.off()
















    pdf("velocity.umap.pdf",w=12,h=8)
    p1
    dev.off()

    png("velocity.umap.png",width=1200,height=800)
    p1
    dev.off()




    emb <- RunVelocity(
        object = emb,
        deltaT = 1,
        kcells = 25,
        fit.quantile = 0.02) 

    ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
    names(x = ident.colors) <- levels(x = bm)
    cell.colors <- ident.colors[Idents(object = bm)]
    names(x = cell.colors) <- colnames(x = bm)

    png("demo.png",width = 1080,height = 1080)
    show.velocity.on.embedding.cor(
        emb = Embeddings(object = bm, reduction = "umap"),
        vel = Tool(object = bm, slot = "RunVelocity"),
        n = 200,
        scale = "sqrt",
        cell.colors = ac(x = cell.colors, alpha = 0.5),
        cex = 0.8,
        arrow.scale = 3,
        show.grid.flow = TRUE,
        min.grid.flow = TRUE,
        min.grid.cell.mass = 0.5,
        grid.n = 40,
        arrow.lwd = 1,
        do.par = FALSE,
        cell.border.alpha = 0.1
    )

    dev.off()


}




#VIP Figure40: OR of tissue preference
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure40")
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

    cellInfo.tb <- tnk_sce@meta.data
    cellInfo.tb = data.table(cellInfo.tb)
    meta.cluster <- cellInfo.tb$clusters

    cellInfo.tb$loc <- as.factor(cellInfo.tb$group)
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
    out.prefix = "cd8"
    sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)


}


# Figure 41: seurat to scanpy
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure41")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    #读入seurat处理后的rds文件
    library(Seurat)
    library(loomR)
    library(hdf5r)

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8|CD4-C1", tnk_sce@meta.data$clusters)]))
    cd8_sce <- subset(tnk_sce,clusters %in% cd8_celltype_array)

    write_h5(cd8_sce,file = "cd8.h5", 
                object.type = 'seurat', 
                assay.name = 'RNA', 
                save.graphs = TRUE, 
                save.scale=FALSE)









}




# Figure42 CD8-C5, CD8-C1, CD8-C4 velocyto analysis
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure42")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    library(SeuratWrappers)
    library(stringr)
    library(ggplot2)
    library(Seurat)
    library(dior)

    c5_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8-C4|CD8-C5", tnk_sce@meta.data$clusters)]))

    c5_sce <- subset(tnk_sce, clusters %in% c5_celltype_array)
    write_h5(c5_sce, 
        file = "c5.h5", 
        object.type = 'seurat', 
        assay.name = 'RNA', 
        save.graphs = TRUE, 
        save.scale=FALSE)

}



# Figure43 CD4 velocyto analysis
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure43")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    library(SeuratWrappers)
    library(stringr)
    library(ggplot2)
    library(Seurat)
    library(dior)

    cd4_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD4", tnk_sce@meta.data$clusters)]))

    cd4_sce <- subset(tnk_sce, clusters %in% cd4_celltype_array)
    write_h5(cd4_sce, 
        file = "cd4.h5", 
        object.type = 'seurat', 
        assay.name = 'RNA', 
        save.graphs = TRUE, 
        save.scale=FALSE)

}





# Figure 45 T cell gsea enrichment
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure45")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    library(clusterProfiler)

    tcell_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8|CD4", tnk_sce@meta.data$clusters)]))

    tcell_sce <- subset(tnk_sce, clusters %in% tcell_celltype_array)

    all_gsea_frame <- data.frame()

    for (current_celltype in tcell_celltype_array) {
        cat(current_celltype,"\n")
        markers <- FindMarkers(tcell_sce,
            ident.1 = current_celltype, 
            ident.2 = tcell_celltype_array[which(tcell_celltype_array != current_celltype)],
            group.by = "clusters", 
            assay = "RNA", 
            slot = "data", 
            logfc.threshold = 0, 
            min.pct = 0
        )

        mdb <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:REACTOME")
        fgsea_sets <- mdb %>% split(x = .$gene_symbol, f = .$gs_name)
        markers$genes = rownames(markers)
        cluster0.genes<- markers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
        ranks<- deframe(cluster0.genes)

        fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
        fgseaRes <- fgseaRes %>% select(c("pathway","padj","NES"))
        fgseaRes$cluster <- current_celltype
        all_gsea_frame <- rbind(all_gsea_frame,fgseaRes)
    }


    saveRDS(all_gsea_frame,file = "gsea.rds",compress = F)

    TCR_signaling_pathway <- c("REACTOME_PHOSPHORYLATION_OF_CD3_AND_TCR_ZETA_CHAINS",
			 "REACTOME_TRANSLOCATION_OF_ZAP_70_TO_IMMUNOLOGICAL_SYNAPSE",
			 "REACTOME_GENERATION_OF_SECOND_MESSENGER_MOLECULES",
			 "REACTOME_DOWNSTREAM_TCR_SIGNALING",
			 "REACTOME_TCR_SIGNALING")
    

    spread_NES_gsea <- all_gsea_frame %>% select(-c("padj")) %>% spread("cluster","NES")
    spread_padj_gsea <- all_gsea_frame %>% select(-c("NES")) %>% spread("cluster","padj")

    row.names(spread_NES_gsea) <- spread_NES_gsea$pathway
    spread_NES_gsea$pathway <- NULL
    row.names(spread_padj_gsea) <- spread_padj_gsea$pathway
    spread_padj_gsea$pathway <- NULL

    # get the top 80% standard deviation reactome pathway
    row.var <- rowSds(as.matrix(spread_NES_gsea))
    f.goi <- row.var >= quantile(row.var,0.8,na.rm = T)
    spread_NES_gsea <- spread_NES_gsea[f.goi,]
    spread_padj_gsea <- spread_padj_gsea[f.goi,]

    # make the NES matrix between -3 and 3
    zz <- 3.0
    spread_gsea[ spread_gsea > zz] <- zz
    spread_gsea[ spread_gsea < -zz] <- -zz

    Heatmap(spread_gsea)



}


# Figure 46 T cell GEO survival
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure45")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }

    library(GEOquery)
    gset = getGEO(GEO='GSE10143', destdir=".",getGPL = F) 

    saveRDS(gset, file = "GSE10143.rdata", compress = F)

    ## 表达矩阵提取
    exprSet <- exprs(gset[[1]])
    ## 分组信息等提取
    pData <- pData(gset[[1]])





}


# Figure 47 CD8 Heatmap
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure47")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    tnk_cellmarker = "
        Group1	GZMK
        Group1	CXCR5
        Group1	CCR4
        Group1	CD28
        Group1	CXCR3
        Group1	GZMH
        Group1	CD27
        Group1	HLA-DRB1
        Group2	KLRD1
        Group2	TYROBP
        Group2	KIR2DL3
        Group2	KIR2DL1
        Group2	KIR3DL1
        Group2	KIR3DL2
        Group2	CD160
        Group2	EOMES
        Group2	TXK
        Group2	KLRC1
        Group2	KIR2DL4
        Group3	PDCD1
        Group3	CXCL13
        Group3	LAYN
        Group4	ZNF683
        Group4	CD52
        Group4	HOPX
        Group4	ID2
        Group4	CXCR6
        Group4	XCL1
        Group4	XCL2
        Group5	CX3CR1
        Group5	FGFBP2
        Group5	FCGR3A
        Group6	SLC4A10
        Group6	KLRB1
        Group6	TMIGD2
        Group6	RORA
        Group6	RORC
        Group6	ZBTB16
        Group6	IL26
        Group6	IL17A
        Group6	IL23R"

	
	cell_marker_frame <- read.table(
        text = gsub("\t\t","",tnk_cellmarker),
        sep = "\t", 
        col.names = c("celltype", "gene")
    )

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8|MAIT", tnk_sce@meta.data$clusters)]))
    cd8_sce <- subset(tnk_sce,clusters %in% cd8_celltype_array)

    Idents(cd8_sce) <- "seurat_clusters"

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$gene),]

    AverageExp <- AverageExpression(cd8_sce,
        features = cell_marker_frame$gene, 
        group.by = "clusters"
    )

    expr <- AverageExp$RNA

    col_fun = colorRamp2(c(-2, 0, 2), c("#6E6BA9","#FDF8EE","#9D1F28"))

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

    expr <- expr[, sort(colnames(expr))]

    celltype_group <- c()
    for (celltype in colnames(expr)) {
        celltype_group <- c(celltype_group, unlist(strsplit(celltype, "-"))[1])
    }
    

    celltype_color_structure <- structure(color_array[1:(length(unique(cell_marker_frame$celltype)))],
                                        names = unique(cell_marker_frame$celltype)) 

    ha = rowAnnotation(labels = cell_marker_frame$celltype,
					show_legend = F,
                    col = list(bar = celltype_color_structure))


    ht <- Heatmap(as.matrix(expr),
            col = col_fun,
            right_annotation = ha,
            rect_gp = gpar(col = "white", lwd = 3),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 16),
            column_split = celltype_group,
            # remove the title showed on the top of the heatmap (due to column_split will show the title)
            column_title = NULL,
            row_split = cell_marker_frame$celltype,
            row_title_side = "right",
            row_title_rot = 0, 
            row_title_gp = gpar(fontsize = 0),
			heatmap_legend_param = list(title = "",
				grid_width = unit(1, "cm"),
				legend_height = unit(8, "cm"),
                title_gp = gpar(fontsize = 20), 
                labels_gp = gpar(fontsize = 30)),
            cluster_columns = F,
            cluster_rows = F)

    png("TNKMarkerGeneHeatmapByCelltype.png",width = 1080,height = 1080 * 1.5)
    draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

    pdf("TNKMarkerGeneHeatmapByCelltype.pdf",width = 8,height = 18)
    draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()



}



# Figure 48 CD4 Heatmap
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure48")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    tnk_cellmarker = "
        Group0	TCF7
        Group0	LEF1
        Group0	CCR7
        Group0	SELL
        Group0	MAL
        Group1	GZMK
        Group1	CXCR5
        Group1	CCR4
        Group1	CD28
        Group1	CXCR3
        Group1	GZMH
        Group1	CD27
        Group1	HLA-DRB1
        Group2	KLRD1
        Group2	TYROBP
        Group2	KIR2DL3
        Group2	KIR2DL1
        Group2	KIR3DL1
        Group2	KIR3DL2
        Group2	CD160
        Group2	EOMES
        Group2	TXK
        Group2	KLRC1
        Group2	KIR2DL4
        Group3	PDCD1
        Group3	CXCL13
        Group3	LAYN
        Group4	ZNF683
        Group4	CD52
        Group4	HOPX
        Group4	ID2
        Group4	CXCR6
        Group4	XCL1
        Group4	XCL2
        Group5	TBX21
        Group5	ASCL2
        Group5	CX3CR1
        Group5	KLRG1
        Group6	SLC4A10
        Group6	KLRB1
        Group6	TMIGD2"

	
	cell_marker_frame <- read.table(
        text = gsub("\t\t","",tnk_cellmarker),
        sep = "\t", 
        col.names = c("celltype", "gene")
    )

    cd8_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8|MAIT", tnk_sce@meta.data$clusters)]))
    cd8_sce <- subset(tnk_sce,clusters %in% cd8_celltype_array)

    Idents(cd8_sce) <- "seurat_clusters"

    cell_marker_frame <- cell_marker_frame[!duplicated(cell_marker_frame$gene),]

    AverageExp <- AverageExpression(cd8_sce,
        features = cell_marker_frame$gene, 
        group.by = "clusters"
    )

    expr <- AverageExp$RNA

    col_fun = colorRamp2(c(-2, 0, 2), c("#6E6BA9","#FDF8EE","#9D1F28"))

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

    expr <- expr[, sort(colnames(expr))]

    celltype_group <- c()
    for (celltype in colnames(expr)) {
        celltype_group <- c(celltype_group, unlist(strsplit(celltype, "-"))[1])
    }
    

    celltype_color_structure <- structure(color_array[1:(length(unique(cell_marker_frame$celltype)))],
                                        names = unique(cell_marker_frame$celltype)) 

    ha = rowAnnotation(labels = cell_marker_frame$celltype,
					show_legend = F,
                    col = list(bar = celltype_color_structure))


    ht <- Heatmap(as.matrix(expr),
            col = heatmapColorRamp,
            right_annotation = ha,
            rect_gp = gpar(col = "white", lwd = 3),
            show_row_names = T,
            row_names_side = "left",
            row_names_gp = gpar(fontsize = 20),
            column_names_gp = gpar(fontsize = 16),
            column_split = celltype_group,
            # remove the title showed on the top of the heatmap (due to column_split will show the title)
            column_title = NULL,
            row_split = cell_marker_frame$celltype,
            row_title_side = "right",
            row_title_rot = 0, 
            row_title_gp = gpar(fontsize = 0),
			heatmap_legend_param = list(title = "",
				grid_width = unit(1, "cm"),
				legend_height = unit(8, "cm"),
                title_gp = gpar(fontsize = 20), 
                labels_gp = gpar(fontsize = 30)),
            cluster_columns = F,
            cluster_rows = F)

    png("TNKMarkerGeneHeatmapByCelltype.png",width = 1080,height = 1080 * 1.5)
    draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()

    pdf("TNKMarkerGeneHeatmapByCelltype.pdf",width = 8,height = 18)
    draw(ht, padding = unit(c(20, 20, 20, 20), "mm"))
    dev.off()



}




#VIP Figure 49: T cell shanno equitability index bar plot
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure49")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)


    diversity.norm <- function(x)
    {
        x$Sample <- NULL
        x <- x %>% as.numeric()
	    -sum(ifelse(x>0,x*log2(x),0))/log2(length(x))
    }

    cell_stat_frame <- as.data.frame.array(table(
        as.character(tnk_sce$sample_id),
        as.character(tnk_sce$clusters)
    ))

    t_celltype_array = colnames(cell_stat_frame %>% dplyr::select(starts_with("CD8")|starts_with("CD4")|starts_with("MAIT")))
    cell_stat_frame <- cell_stat_frame[t_celltype_array]
    cell_stat_frame$Total_cell <- rowSums(cell_stat_frame)

    cell_stat_frame$Sample <- row.names(cell_stat_frame)
    row.names(cell_stat_frame) <- NULL
    # r$> str(cell_stat_frame)
    # 'data.frame':   36 obs. of  13 variables:
    # $ CD8-C1    : int  362 124 704 800 497 115 196 123 378 30 ...
    # $ CD8-C2    : int  99 63 363 170 131 48 279 126 129 24 ...
    # $ CD8-C3    : int  94 21 636 140 90 9 175 28 140 10 ...
    # $ CD8-C4    : int  47 46 71 191 79 20 20 22 77 4 ...
    # $ CD8-C5    : int  53 42 111 71 34 14 63 33 118 9 ...
    # $ CD8-C6    : int  180 37 254 129 89 70 274 26 261 19 ...
    # $ CD4-C1    : int  296 561 505 648 318 259 133 529 332 91 ...
    # $ CD4-C2    : int  29 914 83 467 60 100 39 255 53 76 ...
    # $ CD4-C3    : int  16 35 22 52 3 157 19 107 9 41 ...
    # $ MAIT      : int  429 284 641 237 105 73 153 591 139 90 ...
    # $ T-Cycling : int  31 63 55 63 30 14 86 50 39 10 ...
    # $ Total_cell: num  1636 2190 3445 2968 1436 ...
    # $ Sample     : chr  "ZP01N" "ZP01T" "ZP02N" "ZP02T" ...

    for(celltype in t_celltype_array){
        cell_stat_frame[celltype] = cell_stat_frame[celltype] / cell_stat_frame["Total_cell"]
    }

    #cell_stat_frame %>% group_by("Sample") %>% summarise(eh = diversity.norm())
    cell_stat_frame$Total_cell <- NULL
    eh_array <- sapply(split(cell_stat_frame,cell_stat_frame$Sample),diversity.norm)
    # r$> eh_array
    #     ZP01N     ZP01T     ZP02N     ZP02T     ZP03N     ZP03T     ZP04N     ZP04T     ZP05N     ZP05T     ZP06N     ZP06T     ZP07N     ZP07T     ZP08N 
    # 0.8594364 0.6930948 0.8800887 0.8551630 0.8299485 0.8659766 0.8830055 0.7532835 0.8881614 0.8395323 0.8564430 0.8461261 0.8661506 0.8158090 0.8354188 
    #     ZP08T     ZP09N     ZP09T     ZP10N     ZP10T     ZP11N     ZP11T     ZR01N     ZR01T     ZR02N     ZR02T     ZR03N     ZR03T     ZR04N     ZR04T 
    # 0.8828011 0.8317987 0.8134837 0.8813491 0.8832040 0.8667049 0.9117895 0.8856375 0.8809125 0.8526674 0.7525652 0.8545989 0.7389573 0.8946581 0.8443773 
    #     ZR05N     ZR05T     ZR06N     ZR06T     ZR07N     ZR07T 
    # 0.8597354 0.7563940 0.8610771 0.8186599 0.8474037 0.6895852
    eh_frame <- data.frame(sample = names(eh_array),eh = as.numeric(eh_array))
    eh_frame$Group <- paste0(substr(eh_frame$sample,1,2),substr(eh_frame$sample,5,5))
    # r$> head(eh_frame)
    #   sample         eh Group
    # 1  ZP01N  -4871.978   ZPN
    # 2  ZP01T  -6778.816   ZPT
    # 3  ZP02N -11291.538   ZPN
    # 4  ZP02T  -9550.084   ZPT
    # 5  ZP03N  -4201.028   ZPN
    # 6  ZP03T  -2397.894   ZPT

    p1 <- ggplot(eh_frame,aes(x = Group,y = eh,fill = Group)) + 
        geom_boxplot() +
        scale_fill_manual(values = group_cols) +
        theme_bw() +
    xlab("") + 
    ylab("") +
    geom_signif(comparisons = list(c("PHT", "LRT"), c("PHN", "LRN"),c("PHT", "PHN"),c("LRN", "LRT")),
                step_increase = 0.1, 
                map_signif_level = T,
                test.args = "greater",
                test = t.test)

    ggsave("diversity_boxplot.png",p1,width = 8,height = 8,units = "in")
    ggsave("diversity_boxplot.pdf",p1,width = 8,height = 8,units = "in")


    p1 <- ggplot(eh_frame,aes(x = celltype,y = eh,fill = celltype)) + 
        geom_bar(stat = "identity") +
        scale_fill_manual(values = group_cols) +
        geom_text(aes(x = celltype, y = eh, label = eh,vjust = 1),size = 1) +
        theme_bw() +
    xlab("") + 
    ylab("")

    ggsave("diversity_barplot.png",p1,width = 3,height = 2,units = "in")
    ggsave("diversity_barplot.pdf",p1,width = 3,height = 2,units = "in")

}




# Figure 50: Monocle3 analysis
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure50")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    t_celltype_array <- unique(as.character(tnk_sce@meta.data$clusters[grepl("CD8|CD4|MAIT", tnk_sce@meta.data$clusters)]))

    t_sce <- subset(tnk_sce, clusters %in% t_celltype_array)


    # 提取Seurat表达矩阵
    mat = GetAssayData(t_sce,assay = "RNA", slot = "counts")
    cell_metadata = t_sce@meta.data
    gene_annotation = data.frame(gene_short_name = row.names(mat))
    rownames(gene_annotation) <- rownames(mat)
    # r$> head(gene_annotation)
    #           gene_short_name
    # AL627309.1         AL627309.1
    # AP006222.2         AP006222.2
    # RP4-669L17.10   RP4-669L17.10
    # RP11-206L10.3   RP11-206L10.3
    # RP11-206L10.2   RP11-206L10.2
    # RP11-206L10.9   RP11-206L10.9

    cds = new_cell_data_set(expression_data = mat, 
                        cell_metadata = cell_metadata, 
                        gene_metadata = gene_annotation)

    # r$> cds
    # class: cell_data_set 
    # dim: 22664 46571 
    # metadata(1): cds_version
    # assays(1): counts
    # rownames(22664): AL627309.1 AP006222.2 ... LGALS7B AP001626.2
    # rowData names(1): gene_short_name
    # colnames(46571): ZP01N_AAACCTGAGGACGAAA-1 ZP01N_AAACCTGAGTCTTGCA-1 ... ZR07T_TTTGTCAGTGTTTGGT-1 ZR07T_TTTGTCATCGCTAGCG-1
    # colData names(15): orig.ident nCount_RNA ... clusters Size_Factor
    # reducedDimNames(0):
    # mainExpName: NULL
    # altExpNames(0):

    colData(cds) = colData(cds)[,!grepl("^RNA_",colnames(colData(cds)))]
    # r$> colnames(colData(cds))
    # [1] "orig.ident"       "nCount_RNA"       "nFeature_RNA"     "percent.mt"       "percent.ribo"     "log10GenesPerUMI" "S.Score"          "G2M.Score"       
    # [9] "Phase"            "old.ident"        "RNA_snn_res.0.8"  "seurat_clusters"  "RNA_snn_res.1"    "clusters"         "Size_Factor"
    
    # 注意这里增加了PCA字段
    cds = preprocess_cds(cds, num_dim = 50)
    #r$> cds
    # class: cell_data_set 
    # dim: 22664 46571 
    # metadata(1): cds_version
    # assays(1): counts
    # rownames(22664): AL627309.1 AP006222.2 ... LGALS7B AP001626.2
    # rowData names(1): gene_short_name
    # colnames(46571): ZP01N_AAACCTGAGGACGAAA-1 ZP01N_AAACCTGAGTCTTGCA-1 ... ZR07T_TTTGTCAGTGTTTGGT-1 ZR07T_TTTGTCATCGCTAGCG-1
    # colData names(13): orig.ident nCount_RNA ... clusters Size_Factor
    # reducedDimNames(1): PCA
    # mainExpName: NULL
    # altExpNames(0):

    png("plot_pc_variance_explained.png",width = 1080, height = 1080)
    plot_pc_variance_explained(cds)
    dev.off()

    # 这个地方是重点代码，将已经得到的UMAP坐标导入到cds对象中，这样展示的轨迹就是基于已有的UMAP坐标
    reducedDim(cds, "UMAP") = Embeddings(object = t_sce[["umap"]])

    gp <- plot_cells(cds,
            reduction_method="UMAP", 
            color_cells_by = "clusters")
    ggsave("plot_cells.png",plot = gp, width = 8, height = 8)

    cds = cluster_cells(cds)
    cds = learn_graph(cds, 
        use_partition = T,
        verbose = T, 
        learn_graph_control=list(
            minimal_branch_len=10
        ))

    p = plot_cells(cds,
           color_cells_by = "clusters",
           trajectory_graph_color="black",
           trajectory_graph_segment_size=1,
           label_groups_by_cluster=F,
           label_roots=F,
           label_leaves=F,
           label_branch_points=F,
           cell_size=0.5, 
           #graph_label_size=2,
           group_label_size=3,
           rasterize=T) + 
    scale_colour_manual(values = celltype_colors)
    ggsave(sprintf("%s.umap_Monocle3line.pdf", "TCell"),plot = p, width=8, height=8)
    ggsave(sprintf("%s.umap_Monocle3line.png", "TCell"),plot = p, width=8, height=8)

    start = c("CD4-C1")
    #
    closest_vertex = cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex = as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes = igraph::V(principal_graph(cds)[["UMAP"]])$name
    #
    flag = closest_vertex[as.character(colData(cds)$clusters) %in% start,]
    flag = as.numeric(names(which.max(table( flag ))))
    root_pr_nodes = root_pr_nodes[flag]
    #
    cds = order_cells(cds, root_pr_nodes=root_pr_nodes )
    ##
    width=5.5
    height=4
    #
    p = plot_cells(cds,
            color_cells_by = "pseudotime",
            label_cell_groups=F,
            label_groups_by_cluster=F,
            label_roots=F,
            label_leaves=F,
            label_branch_points=F,
            cell_size=0.5, 
            group_label_size=3,
            rasterize=F)
    #
    ggsave(p, file = "Tcell_monocle3_pseudotime.png", width=8, height=8)
    ggsave(p, file = "Tcell_monocle3_pseudotime.pdf", width=8, height=8)



}





# Figure 51 TRM and TEX survival analysis from TCGA
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

    # input parameter
    global_seurat_rds = "/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds"
    subcluster_seurat_rds = "/P02_TNKCell/P999_TNK_Analysis/Figure0/TNK.rds"
    celltype_1 = "CD8-C4"
    celltype_2 = "CD8-C5"

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

    # load global Seurat object
    global_sce <- readRDS(global_seurat_rds)

    # load subtype Seurat object
    subcluster_sce <- readRDS(subcluster_seurat_rds)

    # get cell id in the specific subpopulations 1 from global cell seurat
    specific_cells <- Cells(subcluster_sce[,subcluster_sce@meta.data$clusters == celltype_1])
    global_sce@meta.data$yesORno <- ifelse(Cells(global_sce) %in% specific_cells, "yes", "no")

    # find markers (specific subpopulations vs global cells)
    markers <- FindMarkers(global_sce, ident.1 = "yes", 
                        group.by = 'yesORno',
                        only.pos = TRUE)
    markers$Gene.name.uniq <- rownames(markers)

    # r$> head(markers)
    #         p_val avg_log2FC pct.1 pct.2 p_val_adj Gene.name.uniq
    # TNFRSF9     0   3.725037 0.588 0.047         0        TNFRSF9
    # NBL1        0   2.209794 0.225 0.036         0           NBL1
    # LCK         0   1.358431 0.926 0.395         0            LCK
    # CSF1        0   2.226701 0.437 0.100         0           CSF1
    # CD2         0   1.599154 0.932 0.371         0            CD2
    # TTC24       0   1.916794 0.201 0.022         0          TTC24

    celltype_1_global_top50 <- markers %>% top_n(n = 50, wt = avg_log2FC)
    # r$> head(celltype_1_global_top50)
    #         p_val avg_log2FC pct.1 pct.2 p_val_adj Gene.name.uniq
    # TNFRSF9     0   3.725037 0.588 0.047         0        TNFRSF9
    # NBL1        0   2.209794 0.225 0.036         0           NBL1
    # CSF1        0   2.226701 0.437 0.100         0           CSF1
    # TTC24       0   1.916794 0.201 0.022         0          TTC24
    # TNFSF4      0   2.636545 0.197 0.014         0         TNFSF4
    # RGS13       0   1.926395 0.131 0.014         0          RGS13
    saveRDS(celltype_1_global_top50,file = "CD8-C4_top50.rds")



    # get cell id in the specific subpopulations 2 from global cell seurat
    specific_cells <- Cells(subcluster_sce[,subcluster_sce@meta.data$clusters == celltype_2])
    global_sce@meta.data$yesORno <- ifelse(Cells(global_sce) %in% specific_cells, "yes", "no")
    #                          RNA_snn_res.0.8 seurat_clusters RNA_snn_res.1 clusters yesORno
    # ZP01N_AAACCTGAGGACGAAA-1               3               1             1    Tcell      no
    # ZP01N_AAACCTGAGTCTTGCA-1               9               2             2    Tcell      no
    # ZP01N_AAACCTGAGTGACTCT-1               7               8             8       NK      no
    # ZP01N_AAACCTGCATCTACGA-1               9               2             2    Tcell      no
    # ZP01N_AAACCTGGTCCGACGT-1               8               7             7  Myeloid      no
    # ZP01N_AAACCTGGTCGAATCT-1               6               9             9       NK      no

    # find markers (specific subpopulations vs global cells)
    markers <- FindMarkers(global_sce, ident.1 = "yes", 
                        group.by = 'yesORno',
                        only.pos = TRUE)
    markers$Gene.name.uniq <- rownames(markers)

    celltype_2_global_top50 <- markers %>% top_n(n = 50, wt = avg_log2FC)
    # r$> celltype_2_global_top50
    #                   p_val avg_log2FC pct.1 pct.2     p_val_adj Gene.name.uniq
    # ZNF683         0.000000e+00   3.916805 0.624 0.050  0.000000e+00         ZNF683
    # FCRL6          0.000000e+00   1.747362 0.639 0.192  0.000000e+00          FCRL6
    # XCL2           0.000000e+00   1.727555 0.694 0.192  0.000000e+00           XCL2
    # XCL1           0.000000e+00   2.416668 0.735 0.152  0.000000e+00           XCL1
    # AC092580.4     0.000000e+00   2.048075 0.471 0.118  0.000000e+00     AC092580.4
    # CD8A           0.000000e+00   1.626264 0.682 0.240  0.000000e+00           CD8A

    saveRDS(celltype_2_global_top50,file = "CD8-C5_top50.rds")

    #celltype_1_global_top50 <- readRDS("CD8-C4_top50.rds")

    #celltype_2_global_top50 <- readRDS("CD8-C5_top50.rds")

    # 注意第一个参数LIHC_expr必须转成矩阵
    # 注意gsva第二个参数是一个list: key是细胞类型，值是基因数组
    celltype_1_gs_exp <- gsva(as.matrix(LIHC_expr), list("CD8-C4" = celltype_1_global_top50$Gene.name.uniq), kcdf = "Poisson", min.sz = 10)
    celltype_2_gs_exp <- gsva(as.matrix(LIHC_expr), list("CD8-C5" = celltype_2_global_top50$Gene.name.uniq), kcdf = "Poisson", min.sz = 10)

    saveRDS(celltype_1_gs_exp,file = "CD8-C4_gsva.rds")
    saveRDS(celltype_2_gs_exp,file = "CD8-C5_gsva.rds")

    # 对gsva结果进行标准化和归一化
    gsva_celltype1 <- t(scale(t(celltype_1_gs_exp)))

    normalization<-function(x){
        return((x-min(x))/(max(x)-min(x)))
    }
    gsva_celltype1 <- normalization(gsva_celltype1)

    # 计算一下CD8-C4和CD8-C5在TCGA中GSVA的相关性

    gsva_celltype1 = as.data.frame(t(gsva_celltype1))
    # r$> head(gsva_celltype1)
    #                     CD8-C4
    # TCGA-DD-A4NG-01A 0.2203436
    # TCGA-G3-AAV4-01A 0.2242275
    # TCGA-2Y-A9H1-01A 0.4248803
    # TCGA-BC-A10Y-01A 0.1634615
    # TCGA-K7-AAU7-01A 0.5643864
    # TCGA-BC-A10W-01A 0.2021251

    # 对gsva结果进行标准化和归一化
    gsva_celltype2 <- t(scale(t(celltype_2_gs_exp)))

    normalization<-function(x){
        return((x-min(x))/(max(x)-min(x)))
    }
    gsva_celltype2 <- normalization(gsva_celltype2)

    # 计算一下CD8-C4和CD8-C5在TCGA中GSVA的相关性

    gsva_celltype2 = as.data.frame(t(gsva_celltype2))

    # gsva_frame <- cbind(gsva_celltype1,gsva_celltype2)

    # cor.test(gsva_frame[["CD8-C4"]],gsva_frame[["CD8-C5"]])







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









# Figure 52 TRM and TEX survival analysis from GEO
if (F) {

    library(data.table)
    library(magrittr)
    library(biomaRt)
    library(limma)
    library(survival)
    library(survminer)
    library(dplyr)
    library(Seurat)

    fig_outdir <- paste0(outdir, "/", "Figure33")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    # input parameter
    global_seurat_rds = "/P01_GlobalCell/P999_GlobalCell_Analysis/Figure0/HCC.rds"
    subcluster_seurat_rds = "/P02_TNKCell/P999_TNK_Analysis/Figure0/TNK.rds"
    celltype = "CD8-C1"

    # load TCGA data: LIHC_expr: LIHC TCGA expression matrix;
    # LIHC_survival: LIHC TCGA survival information, sample/OS/OS.time;

    LIHC_expr <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.expr.rds")
    LIHC_survival <- readRDS("/P02_TNKCell/P999_TNK_Analysis/data/LIHC.survival.rds")

    # load global Seurat object
    global_sce <- readRDS(global_seurat_rds)

    # load subtype Seurat object
    subcluster_sce <- readRDS(subcluster_seurat_rds)

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






# Figure 53: Euclidean distance between PHT&PHN vs LRTvsLRN

if (T) {
    fig_outdir <- paste0(outdir, "/", "Figure53")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    celltype <- "
        CD8-C2
        CD8-C3
        CD8-C4
        CD8-C5
        CD8-C6
        CD4-C1
        CD4-C2
        CD4-C3
        MAIT"

    celltype_group_frame <- read.table(
        text = gsub("        ","",celltype),
        sep = "\t", 
        col.names = c("V1")
    )

    immune_celltype_array <- celltype_group_frame$V1

    cell_stat_frame <- as.data.frame.array(table(
        as.character(tnk_sce$sample_id),
        as.character(tnk_sce$clusters)
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

}








