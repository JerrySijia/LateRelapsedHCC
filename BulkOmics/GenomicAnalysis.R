
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(dplyr)

wes_outdir <- "/WES"
wgs_outdir <- "/WGS"


group_cols <- c(
    PHN = "#ef9020",
    PHT = "#00af3e", 
    LRN = "#0081b4", 
    LRT = "#CD1E24"
)


#VIP Figure2A TMB Comparison BoxPlot between ZPT and ZRT
if (F) {
    fig_outdir <- paste0(wes_outdir, "/", "Figure1")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    tmb_frame <- read.table(
        file = "/WES/P01_WES_Basic/P00_WES_SomaticMutationProcess/TMB_result.tsv",
        sep = "\t", 
        header = T, 
        stringsAsFactors = F
    )

    tmb_frame["group"] <- paste0(substr(tmb_frame$SampleId, 0, 2), substr(tmb_frame$SampleId,5,5))
    tmb_frame$group <- gsub("ZP","PH",tmb_frame$group)
    tmb_frame$group <- gsub("ZR","LR",tmb_frame$group)
    
    tmb_frame$group <- factor(tmb_frame$group,levels = c("PHT","LRT"))

    compaired <- list(c("PHT", "LRT"))

    gp <- ggplot(tmb_frame, aes(x = group, y = TMB, fill = group)) +
        geom_boxplot(width = 0.2) +
        scale_fill_manual(values = group_cols[compaired[[1]]]) +
        xlab("") +
        ylab("Tumor Mutation Burden") +
        theme(panel.background = element_blank(),
                panel.grid = element_blank(),
                panel.border = element_rect(color = "black",fill = NA)) +
        geom_signif(comparisons = compaired,
                    step_increase = 0.3,
                    map_signif_level = F,
                    test = t.test)

    ggsave(filename = "TMBComparisonBoxPlot.pdf",plot = gp, width = 4,height = 6)
    ggsave(filename = "TMBComparisonBoxPlot.png",plot = gp, width = 4,height = 6)

}


#VIP HRD Comparison BoxPlot between ZPT and ZRT

if (F) {
    fig_outdir <- paste0(wes_outdir, "/", "Figure2")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    hdr_frame <- read.table(
        file = "/WES/P01_WES_Basic/P02_WES_SequenzaByCohort/hrd_summary.txt",
        sep = "\t", 
        header = T, 
        stringsAsFactors = F
    )
    hdr_frame$SampleId <- gsub("ZP","PH",hdr_frame$SampleId)
    hdr_frame$SampleId <- gsub("ZR","LR",hdr_frame$SampleId)
    hdr_frame["group"] <- paste0(
        substr(hdr_frame$SampleId, 0, 2),
        substr(hdr_frame$SampleId, 5, 5)
    )

    compaired <- list(c("PHT", "LRT"))

    hrdsum <- ggplot(hdr_frame,aes(x = group,y = HRDSum,fill = group)) + 
    geom_boxplot(width = 0.2) +
    scale_fill_manual(values = group_cols[compaired[[1]]]) +
    xlab("") +
    theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black",fill = NA)) +
    geom_signif(comparisons = compaired,
                step_increase = 0.3,
                map_signif_level = F,
                test = t.test)

    hrd <- ggplot(hdr_frame,aes(x = group,y = HRD,fill = group)) + 
    geom_boxplot(width = 0.2) +
    scale_fill_manual(values = group_cols[compaired[[1]]]) +
    xlab("") +
    theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black",fill = NA)) +
    geom_signif(comparisons = compaired,
                step_increase = 0.3,
                map_signif_level = F,
                test = t.test)


    lst <- ggplot(hdr_frame,aes(x = group,y = LST,fill = group)) + 
    geom_boxplot(width = 0.2) +
    scale_fill_manual(values = group_cols[compaired[[1]]]) +
    xlab("") +
    theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black",fill = NA)) +
    geom_signif(comparisons = compaired,
                step_increase = 0.3,
                map_signif_level = F,
                test = t.test)


    teloai <- ggplot(hdr_frame,aes(x = group,y = TeloAI,fill = group)) + 
    geom_boxplot(width = 0.2) +
    scale_fill_manual(values = group_cols[compaired[[1]]]) +
    xlab("") +
    theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black",fill = NA)) +
    geom_signif(comparisons = compaired,
                step_increase = 0.3,
                map_signif_level = F,
                test = t.test)

    gg_list <- list()
    gg_list[[1]] <- hrdsum
    gg_list[[2]] <- hrd
    gg_list[[3]] <- lst
    gg_list[[4]] <- teloai
    png("HrdComparisonBoxPlot.png",width = 1024,height = 1024)
    grid.arrange(grobs = gg_list,nrow = 2)
    dev.off()
    pdf("HrdComparisonBoxPlot.pdf",width = 12,height = 12)
    grid.arrange(grobs = gg_list,nrow = 2)
    dev.off()
}


# Figure 3 Genomic Landscape between ZPT and ZRT
if (F) {
    fig_outdir <- paste0(wes_outdir, "/", "Figure3")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)



    # create multi-omic schema pie chart
    # six omics
    data <- data.frame(
        group = c("WES","RNA-Seq","WGS","scRNA-Seq","scATAC-Seq","scTCR-Seq"),
        value = c(5,5,5,5,5,5)
    )

    data <- data %>% 
        arrange(desc(group)) %>%
        mutate(prop = value / sum(data$value) *100) %>%
        mutate(ypos = cumsum(prop)- 0.5*prop )

    gp <- ggplot(data, aes(x = "", y = prop, fill = group)) +
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar("y", start=0) +
        theme_void() + 
        theme(legend.position="none") +
        scale_fill_brewer(palette="Set1")

    ggsave(gp, file = "six_pie.svg", width = 8, height = 8)


    # five omics
    data <- data.frame(
        group = c("WES","RNA-Seq","WGS","scRNA-Seq","scTCR-Seq"),
        value = c(5,5,5,5,5)
    )

    data <- data %>% 
        arrange(desc(group)) %>%
        mutate(prop = value / sum(data$value) *100) %>%
        mutate(ypos = cumsum(prop)- 0.5*prop )

    gp <- ggplot(data, aes(x = "", y = prop, fill = group)) +
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar("y", start=0) +
        theme_void() + 
        theme(legend.position="none") +
        scale_fill_brewer(palette="Set1")

    ggsave(gp, file = "five_pie.svg", width = 8, height = 8)



    # genomic SNV heatmap
    #path <- "/WES/P01_WES_Basic/varscan_snv"
    all_sample_snv_frame <- data.frame()
    for (sample_id in list.files(paste0(wes_outdir,"/","P01_WES_Basic/varscan_snv"))) {
        snv_frame <- read.table(
            file = paste0(wes_outdir, "/", "P01_WES_Basic/varscan_snv", "/",sample_id,"/", sample_id, ".varscan_snv.reformated.hg19_multianno.txt"),
            sep = "\t", 
            header = T, 
            stringsAsFactors = F
        )
        snv_frame <- snv_frame %>%
            distinct(Gene.refGene, .keep_all = T) %>%
            filter(Func.refGene %in% c("exonic")) %>%
            dplyr::select(Gene.refGene, ExonicFunc.refGene)
        colnames(snv_frame) <- c("GeneName",sample_id)
        if (nrow(all_sample_snv_frame) == 0) {
            all_sample_snv_frame <- snv_frame
        } else {
            all_sample_snv_frame <- merge(all_sample_snv_frame,snv_frame,by = "GeneName",all= TRUE)
        }
    }


    all_sample_snv_frame = replace(all_sample_snv_frame,is.na(all_sample_snv_frame),"no")
    all_sample_snv_frame <- all_sample_snv_frame %>% mutate(across(everything(), ~ gsub("^synonymous SNV|unknown","no",.)))
    row.names(all_sample_snv_frame) <- all_sample_snv_frame$GeneName
    all_sample_snv_frame$GeneName <- NULL

    #all_sample_snv_frame <- all_sample_snv_frame[, grepl("ZR",colnames(all_sample_snv_frame))]
    all_sample_snv_frame <- all_sample_snv_frame[order(rowSums(all_sample_snv_frame != "no"), decreasing = T), ]
    all_sample_snv_frame <- all_sample_snv_frame[1:10, ]
    all_sample_snv_frame <- all_sample_snv_frame[,order(
        all_sample_snv_frame[1, ],
        all_sample_snv_frame[2, ], 
        all_sample_snv_frame[3, ], 
        all_sample_snv_frame[4, ], 
        all_sample_snv_frame[5, ],
        decreasing = T
    )]

    # right mutation frequency barplot
    mutation_freq_array <- rowSums(all_sample_snv_frame != "no") / ncol(all_sample_snv_frame)
    mut_freq_ha <- HeatmapAnnotation(
        snv = anno_barplot(mutation_freq_array,
            gp = gpar(
                fill = "#174B89",
                col = "#174B89",
                fontfamily = "Arial"
            ),
            ylim = c(0,0.4),
            axis_param = list(side = "bottom", at = c(0, 0.3)),
        ), 
        width = unit(4, "cm"),
        which = "row",
        annotation_label = "% mutated"
    )

    # top TMB barplot
    tmb_frame <- read.table(file = "/WES/P01_WES_Basic/P00_WES_SomaticMutationProcess/TMB_result.tsv", sep = "\t", header = T, stringsAsFactors = F)
    tmb_array <- tmb_frame$TMB
    names(tmb_array) <- tmb_frame$SampleId
    tmb_array <- as.numeric(tmb_array[colnames(all_sample_snv_frame)])

    tmb_ha = HeatmapAnnotation(TMB = anno_barplot(tmb_array,
        bar_width = 0.6,
        height = unit(3, "cm"),
        gp = gpar(fill = "#DF536B", col = "white")
    ), annotation_name_gp = gpar(fontsize = 16))


    ht <- Heatmap(as.matrix(all_sample_snv_frame),
        name = "landscape",
        rect_gp = gpar(col = "black", lwd = 1),
        show_heatmap_legend = FALSE,
        col = structure(c("#327BB7", 
                    "#CD1E24", 
                    "#1E843F", 
                    "#EFEFEF", 
                    "#B968A5"),
            names = c(
                "nonsynonymous SNV",
                "stopgain",
                "stoploss",
                "no",
                "frameshift substitution"
            )
        ),
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 12),
        show_column_names = TRUE,
        height = unit(7, "cm"),
        bottom_annotation = NULL,
        top_annotation = tmb_ha,
        right_annotation = mut_freq_ha,
        cluster_rows = FALSE,
        cluster_columns = FALSE
    )
    
    library(grImport)

    sample_omics <- "
        ZP01T	/WES/Figure3/five_pie.svg
        ZP02T	/WES/Figure3/six_pie.svg
        ZP03T	/WES/Figure3/six_pie.svg
        ZP04T	/WES/Figure3/six_pie.svg
        ZP05T	/WES/Figure3/six_pie.svg
        ZP06T	/WES/Figure3/six_pie.svg
        ZP07T	/WES/Figure3/six_pie.svg
        ZP08T	/WES/Figure3/six_pie.svg
        ZP09T	/WES/Figure3/six_pie.svg
        ZP10T	/WES/Figure3/six_pie.svg
        ZP11T	/WES/Figure3/six_pie.svg
        ZR01T	/WES/Figure3/six_pie.svg
        ZR02T	/WES/Figure3/six_pie.svg
        ZR03T	/WES/Figure3/six_pie.svg
        ZR04T	/WES/Figure3/six_pie.svg
        ZR05T	/WES/Figure3/six_pie.svg
        ZR06T	/WES/Figure3/six_pie.svg
        ZR07T	/WES/Figure3/six_pie.svg"

	sample_omics_frame <- read.table(
		text = gsub("        ","",sample_omics),
		sep = "\t",
		row.names = 1, 
		col.names = c("sample_id","svg")
	)

	omics_image_pdf_array <- sample_omics_frame[colnames(all_sample_snv_frame),]

    #omics_image_pdf_array = rep("/WES/Figure3/pie.svg", 18)
    omics_ha = HeatmapAnnotation(
        Datatype = anno_image(omics_image_pdf_array, border = FALSE, space = unit(0.1, "mm")),
        annotation_name_gp = gpar(fontsize = 12),
		annotation_name_side = "left"
    )

    # Clinical information Heatmap
    clinical_info_frame <- read.table(
        file = "/WES/Figure3/data/clinical_info.tsv",
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE
    )

    clinical_info_frame <- subset(clinical_info_frame,
        select = c("Patient", "Sex", "Age", "Status", "Stage")
    )

    # age divided into ranges
    clinical_info_frame$Age <- ifelse(clinical_info_frame$Age < 60,
        ifelse(clinical_info_frame$Age > 70, ">70", "60-70"), "<60"
    )

    

    rownames(clinical_info_frame) <- clinical_info_frame$Patient
    clinical_info_frame <- clinical_info_frame[substr(colnames(all_sample_snv_frame),0,4), ]
    clinical_info_frame <- subset(clinical_info_frame,select = c(-Patient))
    clinical_info_frame <- clinical_info_frame[]
    clinical_info_frame <- t(clinical_info_frame)

    color_array = c()
    label_array = c()

    grid_NA_color = "grey"
    grid_border_color = "white"


    # Sex
    color_array <- c(color_array,c("#9A3334","#625D9E"))
    label_array <- c(label_array,c("Male","Female"))

    # Age
    color_array <- c(color_array,c("#E95C59","#68A180","orange"))
    label_array <- c(label_array,c("<60","60-70",">70"))

    # Status
    color_array <- c(color_array,c("#BD956A","#6778AE"))
    label_array <- c(label_array,c("primary","relapsed"))

    # Stage
    color_array <- c(color_array,c("#ef9020","#00af3e","#0081b4"))
    label_array <- c(label_array,c("T1b","T2","T4"))

    # Clinical information Heatmap
    clinical_info_ht <- Heatmap(clinical_info_frame,
                                col = structure(color_array,names = label_array),
                                rect_gp = gpar(col = grid_border_color,lwd = 3),
                                show_column_dend = FALSE,
                                show_column_names = T,
                                column_names_rot = 60,
                                column_names_gp = gpar(fontsize = 12),
                                show_row_names = T,
                                row_names_gp = gpar(fontsize = 12),
                                row_names_side = "left",
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                top_annotation = omics_ha,
                                show_heatmap_legend = FALSE,
                                height = unit(4,"cm"),
                                na_col = "grey")


    png("genomic_landscape.png",width = 1024,height = 1024)
    draw(ht %v% clinical_info_ht)
    dev.off()

    pdf("genomic_landscape.pdf",width = 8,height = 12)
    draw(ht %v% clinical_info_ht)
    dev.off()

    # legend plot
    library(ComplexHeatmap)
    mutation_lgd = Legend(at = c("nonsynonymous SNV","stopgain","stoploss","no","frameshift substitution"), 
        title = "Mutation classification", 
        legend_gp = gpar(fill = c("#327BB7", "#CD1E24", "#1E843F", "#EFEFEF", "#B968A5")))
    omics_lgd = Legend(at = c("WES","RNA-Seq","WGS","scRNA-Seq","scATAC-Seq","scTCR-Seq"), 
        title = "Omics", legend_gp = gpar(fill = c("#FC822A","#E12626","#FFFF5A","#52AE55","#3B7DB4","#964F9F")))
    sex_lgd = Legend(at = c("Male","Female"), title = "Sex", legend_gp = gpar(fill = c("#9A3334","#625D9E")))
    age_lgd = Legend(at = c("<60","60-70",">70"), title = "Age", legend_gp = gpar(fill = c("#E95C59","#68A180","orange")))
    gender_lgd = Legend(at = c("Male","Female"), title = "Gender", legend_gp = gpar(fill = c("#BD956A","#6778AE")))
    type_lgd = Legend(at = c("primary","late relapsed"), title = "Sample type", legend_gp = gpar(fill = c("#9A3334","#625D9E")))
    stage_lgd = Legend(at = c("T1b","T2","T4"), title = "Stage", legend_gp = gpar(fill = c("#ef9020","#00af3e","#0081b4")))

    pdf("legend.pdf",w = 8,h = 8)
    draw(packLegend(mutation_lgd,omics_lgd,sex_lgd, age_lgd, gender_lgd,type_lgd,stage_lgd, 
        max_height = unit(10, "cm"),
        column_gap = unit(1, "cm")))
    dev.off()

}





# Figure 4 Purity comparison between ZPT and ZRT
if (F) {
    fig_outdir <- paste0(outdir, "/", "Figure4")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    purity_ploidy_frame <- read.table(file = "/WES/P01_WES_Basic/P02_WES_SequenzaByCohort/purity_ploidy_summary.txt",sep = "\t",header = T,stringsAsFactors = F)

    purity_ploidy_frame["group"] <- paste0(substr(purity_ploidy_frame$SampleId, 0, 2),"T")

    compaired <- list(c("ZPT", "ZRT"))

    purity <- ggplot(purity_ploidy_frame,aes(x = group,y = Purity,fill = group)) + 
    geom_boxplot(width = 0.4,outlier.alpha = 0) +
    scale_fill_manual(values = group_cols) +
    xlab("") +
    theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black",fill = NA),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 20)) +
    geom_signif(
        comparisons = compaired,
        step_increase = 0.3,
        map_signif_level = T,
        test = t.test
    ) + guides(fill=FALSE) 

    ploidy <- ggplot(purity_ploidy_frame,aes(x = group,y = Ploidy,fill = group)) + 
    geom_boxplot(width = 0.4,outlier.alpha = 0) +
    scale_fill_manual(values = group_cols) +
    xlab("") +
    theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black",fill = NA),
            axis.text=element_text(size=16),
            axis.title = element_text(size = 20)) +
    geom_signif(
        comparisons = compaired,
        step_increase = 0.3,
        map_signif_level = T,
        test = t.test
    ) +
        guides(fill=FALSE) 



    gg_list <- list()
    gg_list[[1]] <- purity
    gg_list[[2]] <- ploidy
    png("PurityPloidyComparisonBoxPlot.png",width = 1024 * 1.5,height = 1024)
    grid.arrange(grobs = gg_list,nrow = 1)
    dev.off()
    pdf("PurityPloidyComparisonBoxPlot.pdf",width = 9,height = 6)
    grid.arrange(grobs = gg_list,nrow = 1)
    dev.off()

}



#VIP SFigure4C-Figure5 Neoantigen comparison between HPT and LRT
if (F) {
    fig_outdir <- paste0(wes_outdir, "/", "Figure5")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    sample_neoantigen_stat <- "
        PH01T	40
        PH02T	96
        PH03T	128
        PH04T	20
        PH05T	60
        PH06T	61
        PH07T	195
        PH08T	93
        PH09T	81
        PH10T	185
        PH11T	146
        LR01T	71
        LR02T	62
        LR03T	71
        LR04T	2
        LR05T	42
        LR06T	95
        LR07T	33"

    tmb_frame <- read.table(
        text = gsub("        ","",sample_neoantigen_stat),
        sep = "\t", 
        col.names = c("sample_id", "neoantigen_count")
    )

    tmb_frame$TMB <- (tmb_frame$neoantigen_count * 1000000) / 60456963

    tmb_frame["group"] <- substr(tmb_frame$sample_id,0,2)

    compaired <- list(c("PHT", "LRT"))

    gp <- ggplot(tmb_frame, aes(x = group, y = TMB, fill = group)) +
        geom_boxplot(width = 0.2) +
        scale_fill_manual(values = group_cols) +
        xlab("") +
        ylab("Neoantigen Burden") +
    theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black",fill = NA)) +
    geom_signif(comparisons = compaired,
                step_increase = 0.3,
                map_signif_level = F,
                test = t.test)

    ggsave(filename = "NeoantigenTMBComparisonBoxPlot.pdf",plot = gp, width = 6,height = 6)
    ggsave(filename = "NeoantigenTMBComparisonBoxPlot.png",plot = gp, width = 6,height = 6)

}




############################# WGS Analysis ######################################


if (F) {
    fig_outdir <- paste0(wes_outdir, "/", "FigureG1")
    if (!file.exists(fig_outdir)) {
        dir.create(fig_outdir)
    }
    setwd(fig_outdir)

    all_seg_frame <- data.frame()
    for (s in list.files(paste0(wgs_outdir,"/HMMcopy"))) {
        seg_frame <- read.table(paste0(wgs_outdir, "/HMMcopy/", s, "/", s, ".segments.txt"),
            header = T,
            sep = "\t",
            stringsAsFactors = F
        )
        seg_frame <- seg_frame %>%
            filter(!(chr %in% c("X", "Y"))) %>%
            mutate(sample_id = s)
            all_seg_frame <- rbind(all_seg_frame, seg_frame)
    }
cnv <- all_seg_frame

pdf("CNV_heatmap.pdf",width = 18, height = 6)

leng<-23
ploidy <- 2

patient_margin = 0.2
sample_margin = 0.1
plot_cnv_band_width = 0.6

### chromosome length (1-22), centromere position
chrLen<-c(249250621,
          243199373,
          198022430,
          191154276,
          180915260,
          171115067,
          159138663,
          146364022,
          141213431,
          135534747,
          135006516,
          133851895,
          115169878,
          107349540,
          102531392,
          90354753,
          81195210,
          78077248,
          59128983,
          63025520,
          48129895,
          51304566,
          0)
names(chrLen)<-c(1:22)
cumChrLen<-c(0,cumsum(chrLen))[-length(chrLen)-1]
names(cumChrLen)<-c(1:22)
chrpos <- cumChrLen[1:leng-1]+chrLen[1:leng-1]/2


#cnv part
cnv$sample_id <- gsub("ZP","PH",cnv$sample_id)
cnv$sample_id <- gsub("ZR","LR",cnv$sample_id)
cnv$patient_id <- cnv$sample_id

patient_id_array <- as.array(unique(as.character(cnv$patient_id)))

sample_id_array <- as.array(unique(as.character(cnv$sample)))

#plot
plot_height = length(sample_id_array) * 0.5 + (length(patient_id_array) - 1) * 0.2
plot_width = 22 * 2 * 0.5

sample_number <- length(sample_id_array)
patient_number <- length(patient_id_array)

cols<-c("#73C375",
        "#00008B",
        "#F0F0F0",
        "#A50F15",
        "#DE2D26",
        "#FB6A4A")

# 4: "amplification"; 3: "gain"; 2: "neutral"; 1: "loss";-1: "upde"
names(cols)<-c('1','2','3','4','5','6')

par(mar=c(2,1,2,3))
plot(c(cumChrLen[leng],cumChrLen[leng]),
     c(0,(sample_number - patient_number) * sample_margin + (patient_number - 1) * patient_margin + sample_number),
     col="white",
     xlim=c(0,cumChrLen[leng]),
     xaxt="n",
     yaxt="n",
     xlab="",
     ylab="",
     bty="n")

plotted_sample_num = 0
plotted_patient_num = 0

for(i in 1:patient_number){
    patient_frame <- subset(cnv,patient_id == patient_id_array[i])
    patient_sample_id_array <- as.array(unique(as.character(patient_frame$sample_id)))
    rect(0 - 1000000,
         plotted_sample_num + patient_margin * plotted_patient_num,
         cumChrLen[leng] + 100000,
         plotted_sample_num + patient_margin * plotted_patient_num + length(patient_sample_id_array),lwd = 1.5,border = "grey")
    
    text(cumChrLen[leng],
         plotted_sample_num + patient_margin * plotted_patient_num + length(patient_sample_id_array) * 0.5,
         labels = patient_id_array[i],
         pos=4,
         xpd=TRUE,
         cex=1.6)
    
    for (s_id in patient_sample_id_array) {
        rect(0,
             plotted_sample_num + patient_margin * plotted_patient_num + 0.2,
             cumChrLen[leng],
             plotted_sample_num + patient_margin * plotted_patient_num + plot_cnv_band_width + 0.2,
             border = NA,
             col = "#F4F0EF")
        
        patient_sample_frame <- subset(patient_frame,sample_id == s_id)
        patient_sample_frame$start <- patient_sample_frame$start+cumChrLen[as.character(patient_sample_frame$chr)]
        patient_sample_frame$end <- patient_sample_frame$end+cumChrLen[as.character(patient_sample_frame$chr)]
        patient_sample_frame$col <- cols[as.character(patient_sample_frame$state)]
        
        for( j in nrow(patient_sample_frame):1){
            rect(patient_sample_frame$start[j],
                 plotted_sample_num + patient_margin * plotted_patient_num + 0.2,
                 patient_sample_frame$end[j],
                 plotted_sample_num + patient_margin * plotted_patient_num + plot_cnv_band_width + 0.2,
                 col=patient_sample_frame$col[j],
                 border=NA)
        }
        
        
        
        plotted_sample_num = plotted_sample_num + 1
    }
    
    
    for( t in cumChrLen[1:leng] ){
        segments(t,
                 plotted_sample_num - length(patient_sample_id_array) + patient_margin * plotted_patient_num,
                 t,
                 plotted_sample_num - length(patient_sample_id_array) + patient_margin * plotted_patient_num + length(patient_sample_id_array),
                 xpd=TRUE,
                 lwd=2,
                 col = "grey")
    }
    
    
    plotted_patient_num = plotted_patient_num + 1
}

segments(0,
         (sample_number - patient_number) * sample_margin + (patient_number - 1) * patient_margin + sample_number,
         cumChrLen[leng],
         (sample_number - patient_number) * sample_margin + (patient_number - 1) * patient_margin + sample_number,
         lwd=1.5,xpd=TRUE)


segments(0,
         (sample_number - patient_number) * sample_margin + (patient_number - 1) * patient_margin + sample_number + 3,
         cumChrLen[leng],
         (sample_number - patient_number) * sample_margin + (patient_number - 1) * patient_margin + sample_number + 3,
         lwd=1.5,xpd=TRUE)

labels=names(chrpos)
text(chrpos,rep((sample_number - patient_number) * sample_margin + (patient_number - 1) * patient_margin + sample_number,
                length(chrpos)),
                labels=labels,
                family = "",
                pos=3,cex=1.8,xpd=TRUE)

dev.off()

}














































