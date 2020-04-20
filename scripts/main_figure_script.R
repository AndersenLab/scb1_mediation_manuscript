# scb-1 paper figure draft

library(tidyverse)
library(linkagemapping)
library(viridis)

###############################
# IMPORTANT!!!!!
# Set working directory
setwd("~/Dropbox/AndersenLab/LabFolders/Katie/projects/scb1_mediation_manuscript/")
###############################


###############################
# Figure 1
# Linkage mapping results
###############################

# load data
mappingdf <- read.csv("data/FileS3_mappingdf.csv")

# add condtrt
mappingdf <- mappingdf %>%
    dplyr::mutate(condtrt = paste0(condition, "_", trait))

# choose one trait for each drug (highest LOD) that overlaps with scb-1 and show lod plot for all (only chrV?)
chrVtraits <- mappingdf %>%
    dplyr::filter(chr == "V", 
                  ci_l_pos < 11119730,
                  ci_r_pos > 11118721,
                  trait != "median.EXT") %>%
    na.omit() %>%
    dplyr::group_by(condition) %>%
    dplyr::arrange(desc(lod)) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::pull(condtrt)

# linkage mapping function - facet by all conditions
plot_lods <- function(map, tsize = 12) {
    map1 <- map %>% 
        dplyr::group_by(marker, condtrt) %>% 
        dplyr::filter(lod == max(lod))
    cis <- map %>% 
        dplyr::group_by(marker, condtrt) %>% 
        dplyr::mutate(maxlod = max(lod)) %>%
        dplyr::group_by(iteration, condtrt) %>% 
        dplyr::filter(!is.na(var_exp)) %>%
        dplyr::do(head(., n = 1))
    
    totalmap <- NULL
    for(i in unique(cis$condtrt)) {
        drugci <- cis %>%
            dplyr::filter(condtrt == i)
        drugmap <- map1 %>%
            dplyr::filter(condtrt == i)
        map2 <- linkagemapping:::cidefiner(drugci, drugmap)
        totalmap <- rbind(totalmap, map2)
    }
    
    totalmap$condtrt <- gsub("_", "\n", totalmap$condtrt)
    cis$condtrt <- gsub("_", "\n", cis$condtrt)
    
    ggplot2::ggplot(totalmap) + 
        ggplot2::aes(x = pos/1e+06,y = lod) + 
        ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,ymin = 0, ymax = ci_lod), fill = "skyblue", alpha = 0.5) +
        ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06,y = (1.05 * maxlod)), fill = "red", shape = 25,
                                size = tsize/7, show.legend = FALSE, color = "red") + 
        ggplot2::geom_text(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.2 * maxlod), label = paste0(100 *round(var_exp, digits = 3), "%")),
                               size = tsize/5, colour = "black", hjust = "inward") +
        ggplot2::geom_line(size = tsize/25, alpha = 0.85) +
        ggplot2::facet_grid(condition ~ chr, scales = "free") +
        ggplot2::labs(x = "Genomic position (Mb)", y = "LOD") +
        ggplot2::scale_colour_discrete(name = "Mapping\nIteration") +
        ggplot2::theme_bw(tsize) +
        ggplot2::theme(
            axis.text = ggplot2::element_text(color = "black", face = "bold"),
            axis.title = ggplot2::element_text(face = "bold", color = "black"),
            strip.text = ggplot2::element_text(face = "bold", color = "black"),
            plot.title = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(color = NA, size = 0.6))
}

# plot each of 8 traits
plot_lods(mappingdf %>%
              dplyr::filter(condtrt %in% chrVtraits),
          tsize = 10) + # 8
    ggplot2::geom_vline(data = data.frame(chr = "V", condtrt = gsub("_", "\n", chrVtraits)), aes(xintercept = 11.11), linetype = "dashed", color = "navy", size = 0.5)

# save figure
ggsave("figures/figure1_linkageplots.png", height = 9, width = 7.5)


###############################
# Figure 2
# NIL results chrV
###############################

# load data
nil_genotypes <- read.csv("data/FileS5_nil_genotypes.csv")
nildose <- read.csv("data/FileS6_nildose.csv")
HTA_stats <- read.csv("data/FileS7_HTA_stats.csv")

# plot all drugs for median.EXT
trt <- "median.EXT"
tsize <- 10 # 5

# select conditions
nil_selected <- nildose %>%
    dplyr::filter(selected == T) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(condition = ifelse(grepl("0", drug), control, condition),
                  control = ifelse(condition == control, "None", control))

# regress control condition
nil_regressed <- easysorter::regress(nil_selected) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(relative_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = relative_pheno) %>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain %in% c("N2", "ECA1132", "ECA1135") ~ "N2",
                                                 strain %in% c("CB4856", "ECA1133", "ECA1134") ~ "CB",
                                                 TRUE ~ "NIL"),
                  groups = dplyr::case_when(strain %in% c("N2", "ECA1132", "ECA1135", "ECA232") ~ "N2",
                                            strain %in% c("CB4856", "ECA1133", "ECA1134", "ECA1114") ~ "CB",
                                            TRUE ~ "NIL")) 

maxpheno <- nil_regressed %>%
    dplyr::group_by(groups, condition) %>%
    dplyr::summarise(maxphen = max(phenotype))

# stats2
stats2 <- HTA_stats %>%
    dplyr::filter(comparison %in% c("N2-ECA232", "ECA1114-CB4856"),
                  experiment == "NIL",
                  trait == trt) %>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(groups =  ifelse(grepl("N2", comparison), "N2", "CB")) %>%
    dplyr::full_join(maxpheno, by = c("condition", "groups")) %>%
    tidyr::separate(comparison, c("strain1", "strain2"), "-", remove = F) %>%
    dplyr::mutate(yval = maxphen + 0.15,
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"),
                  movex = dplyr::case_when(pval < 0.0001 ~ 0.5-1,
                                           pval < 0.001 ~ 0.5-1.5,
                                           pval < 0.01 ~ 0.5-2,
                                           pval < 0.05 ~ 0.5-3.5,
                                           TRUE ~ 0.5-1.5),
                  cols = ifelse(grepl("N2", comparison), "N2", "CB"))

# plot phenotype
pheno <-  nil_regressed %>%
    # dplyr::full_join(stats, by = c("strain", "condition", "trait")) %>%
    dplyr::mutate(strain = factor(strain, 
                                  levels = rev(c("N2", "ECA232", "ECA1132", "ECA1135", "ECA1114", "ECA1133", "ECA1134", "CB4856")))) %>%
    # labels = c("N2", "CB4856", "ECA232\nN2[V,CB>N2]", "ECA1114\nCB[V,N2>CB]", "ECA1132\nN2[scb-1∆]", "ECA1135\nN2[scb-1∆]", "ECA1133\nCB[scb-1∆]", "ECA1134\nCB[scb-1∆]"))) %>%
    dplyr::group_by(strain, condition) %>%
    dplyr::mutate(phen = max(phenotype) + 0.2) %>%
    ggplot(.) +
    aes(x = strain, y = phenotype, fill = strain_fill) +
    geom_jitter(width = 0.1, size = 0.05) +
    geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.2) +
    ggplot2::geom_segment(data = stats2, aes(x = strain1, xend = strain2, y = yval, yend = yval, color = cols), inherit.aes = F, size = 0.5) +
    ggplot2::geom_text(data = stats2, aes(x = strain1, label = sig, y = yval + .15, hjust = movex, color = cols), size = tsize/4, inherit.aes = F, angle = -90) +
    scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "NIL" = "grey")) +
    scale_color_manual(values = c("N2" = "orange", "CB" = "blue")) +
    scale_y_continuous(limits = c(0, 1.36),
                       breaks = seq(0,1, 0.5),
                       expand=expand_scale(mult=c(0.2,0.2))) +
    theme_bw(tsize) +
    theme(axis.text.x = element_text(face="bold", color="black"),
          axis.title.x = element_text(face="bold", color="black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text = element_text(face = "bold", color = "black"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    coord_flip() +
    labs(x = " ", y = "Relative median optical density") +
    facet_grid(.~condition)

# plot chrV genotype
chrVgeno <- nil_genotypes %>%
    dplyr::filter(chrom == "V") %>%
    dplyr::mutate(chrom = "chrV",
                  sample = factor(sample, levels = rev(c("N2", "ECA232", "ECA1132", "ECA1135", "ECA1114", "ECA1133", "ECA1134", "CB4856")))) %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(~chrom, scales = "free",  space = "free")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(axis.text.x = element_text(face="bold", color="black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(face="bold", color="black"),
          axis.title.y = element_text(face="bold", color="black"),
          strip.text = element_text(face = "bold", color = "black"),
          plot.title = element_text(face="bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    labs(x = "Position (Mb)", y = "") +
    geom_vline(xintercept = 11.11, linetype = "dashed")

cowplot::plot_grid(chrVgeno, pheno, nrow = 1, ncol = 2, rel_widths = c(1, 4.5), align = "h", axis = "b", labels = c("A", "B"))

# save plot
ggsave("figures/figure2_nilplots.png", width = 7.5, height = 2)


###############################
# Figure 3
# eqtl reanalysis
###############################

# get gene locations
eqtlmap <- read.csv("data/FileS9_eqtlmap.csv")

# chr lengths
df_chr_length <- data.frame(
    chr = c("I", 'II', "III", "IV", "V", "X"),
    stop = c(15072434, 15279421, 13783801, 17493829, 20924180, 17718942), 
    start = 0
)

# factor to order chr
eqtlmap$chr <- factor(eqtlmap$chr, levels = c("I", "II", "III", "IV", "V", "X"))
eqtlmap$probe_chr <- factor(eqtlmap$probe_chr, levels = c("X", "V", "IV", "III", "II", "I"))

eqtlmap <- eqtlmap %>%
    dplyr::filter(!is.na(probe_chr))

toplod <- eqtlmap %>%
    dplyr::filter(lod > 20)
lowlod <- eqtlmap %>%
    dplyr::filter(lod < 20)

# plot eQTL peaks
tsize <- 10
eqtlplot <- ggplot(eqtlmap) +
    aes(x = pos / 1e6, y = probe_start / 1e6, color = lod, alpha = 0.8, size = var_exp) +
    geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 1, ymax=1.01), color='transparent', fill='transparent', size =0.1, inherit.aes =  F) +
    geom_point(data = lowlod) +
    geom_point(data = toplod) +
    facet_grid(probe_chr~chr, scales = "free", space = "free") +
    theme_bw(tsize) +
    labs(x = "QTL position (Mb)", y = "Gene position (Mb)") +
    viridis::scale_color_viridis(option = "C", name = "LOD") +
    scale_size_continuous(range = c(0.1,2), guide = "none") +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold", color = "black"),
          strip.text = element_text(face = "bold")) +
    scale_alpha(guide = "none") +
    ggplot2::geom_vline(data = data.frame(chr = "V"), aes(xintercept = 11.11), linetype = "dashed", color = "navy", size = 0.5)

# plot hotspots
# CIS vs. TRANS
cis <- eqtlmap %>%
    dplyr::mutate(dist = abs(pos - probe_start)) %>%
    dplyr::filter(chr == probe_chr,
                  dist < 1e6)

# same chr trans
samechr_trans <- eqtlmap %>%
    dplyr::mutate(dist = abs(pos - probe_start)) %>%
    dplyr::filter(chr == probe_chr, 
                  dist > 1e6)

# diff chr trans
diffchr_trans <- eqtlmap %>%
    dplyr::mutate(dist = abs(pos - probe_start)) %>%
    dplyr::filter(chr != probe_chr)

# all trans
alltrans <- rbind(diffchr_trans, samechr_trans)

# only use trans eQTL for identifying hotspots
linkagemapping::load_cross_obj("N2xCB4856cross_full")

# riail marker conversion
mappos <- qtl::pull.map(N2xCB4856cross_full2, as.table = TRUE) %>%
    dplyr::mutate(id = rownames(.),
                  cM = pos) %>%
    tidyr::separate(id, into = c("chr", "pos"), sep = "_", remove = FALSE) %>%
    dplyr::group_by(chr) %>%
    dplyr::arrange(chr, cM) %>%
    dplyr::mutate(first = first(cM))%>%
    dplyr::mutate(cM = cM - first + 1,
                  pos = as.numeric(pos)) # changing to start at 1 to make more sense

# function to find hotspots
QTLperbin <- function(mappos, uniqueQTL, numcm = 30, pval = .01, bonferroni = "FALSE") {
    
    #round to the nearest bin (first bin is bin "0")
    rounded <- mappos %>%
        dplyr::group_by(chr) %>%
        dplyr::mutate(bin_cM = numcm*(floor(cM/numcm)))
    
    #find cutoffs for each bin
    grouped <- rounded %>%
        dplyr::group_by(chr, bin_cM) %>%
        dplyr::summarise(start = min(pos), end = max(pos)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(all_bin = rownames(.)) %>%
        dplyr::group_by(chr, bin_cM)
    
    
    #Count how many "unique" QTL we found
    allunique <- uniqueQTL %>%
        na.omit() %>%
        dplyr::group_by(chr, drugtrait) %>%
        dplyr::distinct()
    
    #count number of bins in each chr (keeping in mind that bin order starts at "0" not "1")
    binsperchr <- grouped %>%
        dplyr::group_by(chr) %>%
        dplyr::summarize(totalbin = n()) # change to fix chrX problems
    # dplyr::summarize(totalbin = max(bin_cM)/numcm + 1)
    
    totalbins <- sum(binsperchr$totalbin)
    
    #number of QTL per bin based on poisson distribution with mean of...:
    lambda <- nrow(allunique)/totalbins
    
    pval <- ifelse(bonferroni == "TRUE", 1-(pval/totalbins), 1-pval)
    
    
    #plot by bin
    mappedunique <- uniqueQTL %>%
        na.omit()
    
    #add bin information to map
    mappedwithbin <- mappedunique %>%
        dplyr::mutate(bin = 0)
    
    
    for (i in 1:nrow(mappedwithbin)) {
        for(j in 1:nrow(grouped)) {
            mappedwithbin$bin[i] <- 
                ifelse(mappedwithbin$chr[i] == grouped$chr[j] &&
                           mappedwithbin$pos[i] <= grouped$end[j] &&
                           mappedwithbin$pos[i] >= grouped$start[j],
                       grouped$bin_cM[j],
                       mappedwithbin$bin[i])
        }
    }
    
    mappedwithbin <- mappedwithbin %>%
        dplyr::ungroup() %>%
        dplyr::group_by(marker, trait) %>%
        dplyr::distinct(.keep_all = TRUE)
    
    # count number of peaks in each bin
    mappedwithbin <- mappedwithbin %>%
        dplyr::group_by(chr, bin) %>%
        dplyr::mutate(peaks_bin = n()) %>%
        dplyr::select(chr, bin, peaks_bin) %>%
        dplyr::distinct()
    
    peaksperbin <- rounded %>%
        dplyr::left_join(mappedwithbin, by = c("chr", c("bin_cM" = "bin"))) %>%
        dplyr::mutate(peaks_bin = replace_na(peaks_bin, 0))
    
    peaksperbin <- peaksperbin %>%
        dplyr::ungroup() %>%
        dplyr::group_by(chr) %>%
        dplyr::arrange(chr,pos)
    
    plot <- ggplot2::ggplot(peaksperbin) +
        ggplot2::aes(x = pos/1e6, y = peaks_bin) +
        ggplot2::theme_bw() +
        ggplot2::geom_line()+
        ggplot2::geom_hline(yintercept = qpois(pval, lambda), color = "red") +
        ggplot2::facet_grid(.~chr, scales = "free", space = "free") +
        ggplot2::labs(x = "Genomic position (Mb)", y = "Number of peaks in bin") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, face="bold", color="black"),
                       axis.text.y = ggplot2::element_text(size=10, face="bold", color="black"),
                       axis.title.x = ggplot2::element_text(size=12, face="bold", color="black"),
                       axis.title.y = ggplot2::element_text(size=12, face="bold", color="black"),
                       strip.text.x = ggplot2::element_text(size=12, face="bold", color="black"),
                       strip.text.y = ggplot2::element_text(size=12, color="black", face="bold", angle = 0),
                       strip.background = ggplot2::element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
                       plot.title = ggplot2::element_blank())
    
    summary <- list(plot, c(data.frame(c(nrow(allunique)), totalbins, lambda, bonferroni, pval)), peaksperbin, grouped)
    assign("summary", summary, envir = globalenv())
    
}

# add drugtrait
peaks <- alltrans %>%
    dplyr::mutate(drugtrait = trait)

QTLperbin(mappos, peaks, numcm = 5, pval = 0.01, bonferroni = "TRUE")

# how many bins are significant?
which_bin <- summary[[4]] %>%
    dplyr::select(chr, bin_cM, all_bin)

sig <- summary[[3]] %>%
    dplyr::filter(peaks_bin > qpois(summary[[2]]$pval, summary[[2]]$lambda)) %>%
    dplyr::left_join(which_bin) %>% 
    dplyr::distinct(chr, bin_cM, peaks_bin, all_bin)

# draw star over hotspot by hand:
hs <- data.frame(chr = c("II", "IV", "IV", "IV", "V", "V", "V", "X", "X"), pos = c(5,5,10,16,5,11,14,4.5,7), 
                 peaks_bin = c(30, 70, 30, 50, 30, 30, 30, 80, 30), sig = c(F, T, F, T, F, F, F, T, F)) %>%
    dplyr::mutate(pos = pos*1e6)

yval <- -10
cm5 <- summary[[3]] %>%
    dplyr::mutate(hotspot = " ") %>%
    ggplot2::ggplot() +
    ggplot2::aes(x = pos/1e6, y = peaks_bin) +
    geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 1, ymax=1.01), color='transparent', fill='transparent', size =0.1, inherit.aes = F) +
    ggplot2::theme_bw(tsize) +
    ggplot2::geom_line()+
    ggplot2::geom_text(data = hs, label = "*", fontface = "bold", size = 5, aes(color = sig)) +
    ggplot2::geom_segment(data = hs %>% dplyr::filter(chr == "IV"), aes(x = 5, xend = 12, y = yval, yend = yval), size = 1, color = "grey30") + #chrIVL drug hotspot
    ggplot2::geom_segment(data = hs %>% dplyr::filter(chr == "IV"), aes(x = 15, xend = 17, y = yval, yend = yval), size = 1, color = "grey30") + #chrIVR drug hotspot
    ggplot2::geom_segment(data = hs %>% dplyr::filter(chr == "V"), aes(x = 7, xend = 13, y = yval, yend = yval), size = 1, color = "grey30") + #chrV drug hotspot
    ggplot2::scale_color_manual(values = c("TRUE" = "blue", "FALSE" = 'black')) +
    ggplot2::geom_hline(yintercept = qpois(summary[[2]]$pval, summary[[2]]$lambda), color = "red") +
    ggplot2::facet_grid(.~chr, scales = "free", space = "free") +
    ggplot2::labs(x = "QTL position (Mb)", y = "Number of QTL") +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold", color = "black"),
          strip.text = element_text(face = "bold")) +
    ggplot2::geom_vline(data = data.frame(chr = "V"), aes(xintercept = 11.11), linetype = "dashed", color = "navy", size = 0.5)

# plot eqtl and hotspots together
cowplot::plot_grid(eqtlplot, cm5, nrow = 2, align = "v", axis = "lr", labels = c("A", "B"), rel_heights = c(1, 0.3))

# save
ggsave("figures/figure3_eqtlplot.png", height = 8, width = 7.5)


###############################
# Figure 4
# mediation
###############################

# load data
scb1_mediation <- read.csv("data/FileS12_scb1_mediation.csv")
set1_mappingdf <- read.csv("data/FileS11_set1_mappingdf.csv")

# add condtrt
mappingdf <- set1_mappingdf %>%
    dplyr::mutate(condtrt = paste0(condition, "_", trait))

# choose one trait for each drug (highest LOD) that overlaps with scb-1 and show lod plot for all (only chrV?)
chrVtraits <- mappingdf %>%
    dplyr::filter(chr == "V",
                  ci_l_pos < 11119730,
                  ci_r_pos > 11118721) %>%
    na.omit() %>%
    dplyr::filter(drugtrait %in% c("bleomycin_median.EXT", "amsacrine_median.EXT", "cisplatin_median.TOF", "silver_median.norm.EXT"))

tsize <- 10
# plot all three drugs together
scb1_mediation %>%
    dplyr::filter(trait %in% c("bleomycin_median.EXT", "amsacrine_median.EXT", "cisplatin_median.TOF", "silver_median.norm.EXT"))%>%
    dplyr::mutate(abs_est = abs(estimate)) %>%
    tidyr::separate(marker, c("chrom", "pos"), "_") %>%
    dplyr::mutate(pos = as.numeric(pos)) %>%
    dplyr::filter(var == "prop_med") %>%
    dplyr::arrange(pos) %>%
    dplyr::mutate(scb1 = ifelse(probe == "A_12_P104350", "yes", "no"),
                  condition = stringr::str_split_fixed(trait, "_", 2)[,1]) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(q90 = quantile(abs_est, probs = 0.9)[[1]],
                  q95 = quantile(abs_est, probs = 0.95)[[1]],
                  q99 = quantile(abs_est, probs = 0.99)[[1]]) %>%
    ggplot(.) +
    aes(x = pos/1e6, y = abs_est, fill = scb1, shape = scb1, size = scb1) +
    geom_point(aes(alpha = prob), color = "black") +
    scale_alpha_continuous(range = c(1, 0.1)) +
    theme_bw(10) +
    scale_fill_manual(values = c("yes" = "red", "no" = "black")) +
    scale_shape_manual(values = c("yes" = 23, "no" = 21)) +
    scale_size_manual(values = c("yes" = 3, "no" = 1)) +
    labs(x = "Genomic position (Mb)", y = "Mediation estimate") +
    theme(panel.grid = element_blank(),
          legend.position = "none",
          axis.text = element_text(face="bold", color="black"),
          axis.title = element_text(face="bold", color="black"),
          strip.text = element_text(face = "bold", color = "black")) +
    geom_hline(aes(yintercept = q90), color = "grey") +
    facet_wrap(~condition, scales = "free", nrow = 1) +
    geom_vline(data = chrVtraits, aes(xintercept = ci_l_pos/1e6), color = "blue", linetype = "dashed") +
    geom_vline(data = chrVtraits, aes(xintercept = ci_r_pos/1e6), color = "blue", linetype = "dashed")

ggsave("figures/figure4_mediation.png", height = 2.5, width = 7.5)


###############################
# Figure 5
# scb-1 del results
###############################

# load data
scb1_pruned <- read.csv("data/FileS13_scb1_pruned.csv")
nil_genotypes <- read.csv("data/FileS5_nil_genotypes.csv")
HTA_stats <- read.csv("data/FileS7_HTA_stats.csv")

# plot all drugs for median.EXT
trt <- "median.EXT"
tsize <- 10

# regress scb-1
scb1_regressed <- easysorter::regress(scb1_pruned) %>%
    dplyr::filter(strain %in% c("N2", "CB4856", "ECA1134", "ECA1132")) %>%
    dplyr::filter(trait == trt) %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(relative_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T)))) %>%
    dplyr::mutate(phenotype = relative_pheno) %>%
    dplyr::mutate(strain_fill = dplyr::case_when(strain %in% c("N2") ~ "N2",
                                                 strain %in% c("CB4856") ~ "CB",
                                                 TRUE ~ "NIL"),
                  groups = dplyr::case_when(strain %in% c("N2", "ECA1132", "ECA1135", "ECA232") ~ "N2",
                                            strain %in% c("CB4856", "ECA1133", "ECA1134", "ECA1114") ~ "CB",
                                            TRUE ~ "NIL"))

maxpheno <- scb1_regressed %>%
    dplyr::group_by(groups, condition) %>%
    dplyr::summarise(maxphen = max(phenotype))

# stats2
stats2 <- HTA_stats %>%
    dplyr::filter(comparison %in% c("N2-ECA1132", "ECA1134-CB4856"),
                  experiment == "CRISPR",
                  trait == trt) %>%
    dplyr::select(condition, trait, comparison, pval = adj.p.value) %>%
    dplyr::mutate(groups =  ifelse(grepl("N2", comparison), "N2", "CB")) %>%
    dplyr::full_join(maxpheno, by = c("condition", "groups")) %>%
    tidyr::separate(comparison, c("strain1", "strain2"), "-", remove = F) %>%
    dplyr::mutate(yval = maxphen + 0.1,
                  sig = dplyr::case_when(pval < 0.0001 ~ "****",
                                         pval < 0.001 ~ "***",
                                         pval < 0.01 ~ "**",
                                         pval < 0.05 ~ "*",
                                         TRUE ~ "ns"),
                  movex = dplyr::case_when(pval < 0.0001 ~ 0.5-1,
                                           pval < 0.001 ~ 0.5-1.5,
                                           pval < 0.01 ~ 0.5-2,
                                           pval < 0.05 ~ 0.5-3.5,
                                           TRUE ~ 0.5-1.5),
                  cols = ifelse(grepl("N2", comparison), "N2", "CB"))

# plot
pheno <- scb1_regressed %>%
    dplyr::mutate(strain = factor(strain, 
                                  levels = rev(c("N2", "ECA1132", "ECA1134", "CB4856")))) %>%
    dplyr::group_by(strain, condition) %>%
    dplyr::filter(!is.na(strain)) %>%
    dplyr::mutate(phen = max(phenotype) + 0.2) %>%
    ggplot(.) +
    aes(x = strain, y = phenotype, fill = strain_fill) +
    geom_jitter(width = 0.1, size = 0.05) +
    geom_boxplot(outlier.color = NA, alpha = 0.5, size = 0.2) +
    ggplot2::geom_segment(data = stats2, aes(x = strain1, xend = strain2, y = yval, yend = yval, color = cols), inherit.aes = F, size = 0.5) +
    ggplot2::geom_text(data = stats2, aes(x = strain1, label = sig, y = yval+0.15, hjust = movex, color = cols), size = tsize/4, inherit.aes = F, angle = -90) +
    scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "NIL" = "grey")) +
    scale_color_manual(values = c("N2" = "orange", "CB" = "blue")) +
    scale_y_continuous(limits = c(0, 1.36),
                       breaks = seq(0,1, 0.5),
                       expand=expand_scale(mult=c(0.2,0.2))) +
    theme_bw(tsize) +
    theme(axis.text.x = element_text(face="bold", color="black"),
          axis.title.x = element_text(face="bold", color="black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text = element_text(face = "bold", color = "black"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    coord_flip() +
    labs(x = " ", y = "Relative median optical density") +
    facet_grid(.~condition)

# fake ECA1132 geno
eca1132 <- nil_genotypes %>%
    dplyr::filter(sample == "N2") %>%
    dplyr::mutate(sample = "ECA1132")

# fake ECA1133 geno
eca1134 <- nil_genotypes %>%
    dplyr::filter(sample == "CB4856") %>%
    dplyr::mutate(sample = "ECA1134")

# add to nil genotypes
nil_genotypes <- nil_genotypes %>%
    dplyr::bind_rows(eca1134, eca1132)

# plot chrV genotype
chrVgeno <- nil_genotypes %>%
    dplyr::filter(chrom == "V",
                  sample %in% c("N2", "CB4856", "ECA1132", "ECA1134")) %>%
    dplyr::mutate(chrom = "chrV",
                  sample = factor(sample, levels = rev(c("N2", "ECA1132", "ECA1134", "CB4856")))) %>%
    dplyr::distinct() %>%
    ggplot(.)+
    geom_segment(aes(x = start/1e6, y = sample, xend = end/1e6, yend = sample, color = gt_name, size = 2), alpha = 0.7)+
    facet_grid(~chrom, scales = "free",  space = "free")+
    scale_color_manual(values=c("N2"="orange","CB4856"="blue"))+
    theme_bw(tsize) +
    theme(axis.text.x = element_text(face="bold", color="black"),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(face="bold", color="black"),
          axis.title.y = element_text(face="bold", color="black"),
          strip.text = element_text(face = "bold", color = "black"),
          plot.title = element_text(face="bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    labs(x = "Position (Mb)", y = "") +
    geom_vline(xintercept = 11.11, linetype = "dashed") +
    geom_point(data = nil_genotypes %>%
                   dplyr::filter(sample %in% c("ECA1134", "ECA1132"), chrom == "V") %>%
                   dplyr::mutate(chrom = "chrV"),
               aes(x = 11.11, y = sample), shape = 24, color = "black", fill = "grey", size = 3)

cowplot::plot_grid(chrVgeno, pheno, nrow = 1, ncol = 2, rel_widths = c(1, 4.5), align = "h", axis = "b", labels = c("A", "B"))

# save plot
ggsave("figures/figure5_scb1delplot.png", width = 7.5, height = 2)

