# Plot COs in BW248 vs Barke 50K iSelect recombination data
library(tidyverse)
library(multcomp)
library(multcompView)
library(ggpubr)
library(cowplot)

chr_file <- read.table("morex_v3.txt")
load("chromosome_CO_counts.Rdata")

chr_stats <- rename(chr_file,
                    chromosome = V1,
                    chr_start = V2,
                    chr_end = V3,
                    centromere = V4)

chr_stats$centromere_percent <- chr_stats$centromere/chr_stats$chr_end

co_count_final <- merge(co_count_final,
                      chr_stats, 
                      by ="chromosome")

colnames(co_count_final)

# total cos per marker
co_count_final$wt_total_COs <- rowSums(co_count_final[,c(grep("_wt", names(co_count_final)))],
                                       na.rm = T)
co_count_final$ss_total_COs <- rowSums(co_count_final[,c(grep("_ss", names(co_count_final)))],
                                       na.rm = T)

# average crossovers per marker
co_count_final$wt_recomb_avg <- co_count_final$wt_total_COs/length(grep("_wt", names(co_count_final)))
co_count_final$ss_recomb_avg <- co_count_final$ss_total_COs/length(grep("_ss", names(co_count_final)))

# Count occurence of NA per marker (hom regions)
count_na <- function(x) sum(is.na(x)) 
co_count_final$wt_na_sum <- rowSums(is.na(co_count_final[,grep("_wt", names(co_count_final))]))
co_count_final$ss_na_sum <- rowSums(is.na(co_count_final[,grep("_ss", names(co_count_final))]))

# recombination frequency per marker
co_count_final$wt_recomb_freq <- co_count_final$wt_total_COs/(length(grep("_wt", names(co_count_final)))-co_count_final$wt_na_sum)
co_count_final$ss_recomb_freq <- co_count_final$ss_total_COs/(length(grep("_ss", names(co_count_final)))-co_count_final$ss_na_sum)

# work out the per chromosome change in total crossovers ----------------------
chr_total <- co_count_final %>% group_by(chromosome) %>%
  summarise(WT = sum(wt_total_COs),
            SS = sum(ss_total_COs))

chr_total$perc_change <- (chr_total$SS/chr_total$WT)-1

write.csv(chr_total,
          "chromosome_total_crossovers.csv",
          row.names = F)

# work out the per chromosome change in average crossovers ----------------------
chr_avg <- co_count_final %>% group_by(chromosome) %>%
  summarise(WT = sum(wt_recomb_avg),
            SS = sum(ss_recomb_avg))

chr_avg$perc_change <- (chr_avg$SS/chr_avg$WT)-1

write.csv(chr_total,
          "chromosome_average_crossovers.csv",
          row.names = F)

# create percentage bins ------------------------------------------------------
co_count_final$percent <- (co_count_final$position/co_count_final$chr_end)*100

co_count_final$perc_bin <- cut(co_count_final$percent,
                             breaks = seq(0,
                                          100,
                                          by = 2),
                             labels = paste(seq(0,
                                                98,
                                                by = 2),
                                            seq(2,
                                                100,
                                                by = 2),
                                            sep = '-' ))

# Plot a summary of recombination frequency in all chromosomes by percentage bin
# get per position bin CO number
bin_total <- co_count_final %>% group_by(perc_bin) %>%
  summarise(WT = sum(wt_recomb_freq,
                     na.rm = T),
            SS = sum(ss_recomb_freq,
                     na.rm = T))

IBM <- c("#648FFF",
         "#DC267F",
         "#FFB000")

# function to plot every nth x axis label
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}

colnames(bin_total)

colors <- c("WT" = IBM[1],
            "SS" = IBM[2])

all_chr_WT <- ggplot(bin_total, aes(perc_bin)) +
  geom_bar(aes(y=WT,
               fill = "WT"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8,)) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  ylab("Total Recombination Frequency") +
  xlab("Genomic positon percentile") +
  ylim(0,2.2)
all_chr_WT

all_chr_mt <- ggplot(bin_total, aes(perc_bin)) +
  geom_bar(aes(y=SS,
               fill = "SS"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8)) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  ylab("Total Recombination Frequency") +
  xlab("Genomic positon percentile") +
  ylim(0,2.2)
all_chr_mt

plot_grid(all_chr_WT,
          all_chr_mt,
          labels = "auto")

ggsave("total_recomb_freq_side_by_side.png",
       width = 183,
       height = 80,
       units = "mm",
       dpi = 600)

ggsave("total_recomb_freq_side_by_side.svg",
       width = 183,
       height = 80,
       units = "mm")

# Plot recombination frequency per chromosome -------------------------------
# get per position bin total CO number per chromosome
bin_total <- co_count_final %>% group_by(perc_bin, chromosome) %>%
  summarise(WT = sum(wt_recomb_freq,
                     na.rm = T),
            SS = sum(ss_recomb_freq,
                     na.rm = T))

# tidyr hates this column ID for some reason and simply will not work unless 
# we change it.
bin_total$x <- bin_total$perc_bin
colnames(bin_total)
bin_total <- bin_total[,-1]

bin_total <- complete(bin_total,
                      chromosome,
                      x,
                      fill = list(WT = 0,
                                  SS = 0))
colnames(bin_total)

long <- pivot_longer(bin_total,
                     c(3:4),
                     names_to = "Phenotype",
                     values_to = "Crossovers")

colnames(long)

one_pos <- subset(bin_total,
                  chromosome == "chr1H")
two_pos <- subset(bin_total,
                  chromosome == "chr2H")
thr_pos <- subset(bin_total,
                  chromosome == "chr3H")
fou_pos <- subset(bin_total,
                  chromosome == "chr4H")
fiv_pos <- subset(bin_total,
                  chromosome == "chr5H")
six_pos <- subset(bin_total,
                  chromosome == "chr6H")
sev_pos <- subset(bin_total,
                  chromosome == "chr7H")

one_bin_mt <- ggplot(one_pos, aes(x)) +
  geom_bar(aes(y=SS,
               fill = "SS"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "none") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  geom_vline(xintercept = "38-40",
             linetype = 2)
one_bin_mt

one_bin_WT <- ggplot(one_pos, aes(x)) +
  geom_bar(aes(y=WT,
               fill = "WT"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "none") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4))  +
  geom_vline(xintercept = "38-40",
             linetype = 2)
one_bin_WT

two_bin_mt <- ggplot(two_pos, aes(x)) +
  geom_bar(aes(y=SS,
               fill = "SS"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "ntwo") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  geom_vline(xintercept = "44-46",
             linetype = 2)
two_bin_mt

two_bin_WT <- ggplot(two_pos, aes(x)) +
  geom_bar(aes(y=WT,
               fill = "WT"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "ntwo") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4))  +
  geom_vline(xintercept = "44-46",
             linetype = 2)
two_bin_WT

thr_bin_mt <- ggplot(thr_pos, aes(x)) +
  geom_bar(aes(y=SS,
               fill = "SS"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nthree") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  geom_vline(xintercept = "42-44",
             linetype = 2)
thr_bin_mt

thr_bin_WT <- ggplot(thr_pos, aes(x)) +
  geom_bar(aes(y=WT,
               fill = "WT"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nthree") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  geom_vline(xintercept = "42-44",
             linetype = 2)
thr_bin_WT

fou_bin_mt <- ggplot(fou_pos, aes(x)) +
  geom_bar(aes(y=SS,
               fill = "SS"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nfou") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  geom_vline(xintercept = "44-46",
             linetype = 2)
fou_bin_mt

fou_bin_WT <- ggplot(fou_pos, aes(x)) +
  geom_bar(aes(y=WT,
               fill = "WT"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nfou") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4))  +
  geom_vline(xintercept = "44-46",
             linetype = 2)
fou_bin_WT

fiv_bin_mt <- ggplot(fiv_pos, aes(x)) +
  geom_bar(aes(y=SS,
               fill = "SS"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nfiv") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4))  +
  geom_vline(xintercept = "34-36",
             linetype = 2)
fiv_bin_mt

fiv_bin_WT <- ggplot(fiv_pos, aes(x)) +
  geom_bar(aes(y=WT,
               fill = "WT"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nfiv") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4))  +
  geom_vline(xintercept = "34-36",
             linetype = 2)
fiv_bin_WT

six_bin_mt <- ggplot(six_pos, aes(x)) +
  geom_bar(aes(y=SS,
               fill = "SS"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nsix") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4))  +
  geom_vline(xintercept = "44-46",
             linetype = 2)
six_bin_mt

six_bin_WT <- ggplot(six_pos, aes(x)) +
  geom_bar(aes(y=WT,
               fill = "WT"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nsix") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4)) +
  geom_vline(xintercept = "44-46",
             linetype = 2)
six_bin_WT

sev_bin_mt <- ggplot(sev_pos, aes(x)) +
  geom_bar(aes(y=SS,
               fill = "SS"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nsev") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4))   +
  geom_vline(xintercept = "50-52",
             linetype = 2)
sev_bin_mt

sev_bin_WT <- ggplot(sev_pos, aes(x)) +
  geom_bar(aes(y=WT,
               fill = "WT"),
           stat = "identity") +
  scale_fill_manual(values = colors,
                    name = "Phenotype") +
  theme(axis.text.x = element_text(size = 8,
                                   angle = 90),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.position = "nsev") +
  ylim(0,0.75) +
  scale_x_discrete(breaks = every_nth(n = 4))    +
  geom_vline(xintercept = "50-52",
             linetype = 2)
sev_bin_WT

long$chromosome <- as.factor(long$chromosome)

chroms <- levels(long$chromosome)

plot_grid(one_bin_WT,
          one_bin_mt,
          two_bin_WT,
          two_bin_mt,
          thr_bin_WT,
          thr_bin_mt,
          fou_bin_WT,
          fou_bin_mt,
          fiv_bin_WT,
          fiv_bin_mt,
          six_bin_WT,
          six_bin_mt,
          sev_bin_WT,
          sev_bin_mt,
          nrow = 7,
          ncol = 2,
          labels = c("chr1H",
                     "",
                     "chr2H",
                     "",
                     "chr3H",
                     "",
                     "chr4H",
                     "",
                     "chr5H",
                     "",
                     "chr6H",
                     "",
                     "chr7H",
                     ""),
          label_size = 8)

ggsave("total_CO_freq_each_chrom.svg",
       width = 183,
       height = 200,
       units = "mm")

ggsave("total_CO_freq_each_chrom.png",
       width = 183,
       height = 200,
       units = "mm",
       dpi = 600)

# Box plots each chromosome all individual CO counts
colnames(co_count_final)

# want chromosome totals per individual
indiv_totals <- pivot_longer(co_count_final,
                             4:163,
                             values_to = "crossovers",
                             names_to = "individuals")

colnames(indiv_totals)

indiv_totals <- indiv_totals[,c(1,18,19)]

indiv_grp <- indiv_totals %>% group_by(individuals, chromosome) %>%
  summarise(COs = sum(crossovers, na.rm = T))

indiv_grp <- separate(indiv_grp,
                        individuals,
                        c("cross","F1", "F2", "individual", "seg_phenotype"),
                        sep = "_",
                        remove = F)

colnames(indiv_grp)

indiv_grp <- unite(indiv_grp,
                   fam,
                   c(3,4,6),
                   sep = "_",
                   remove = F)

indiv_grp$chromosome <- as.factor(indiv_grp$chromosome)
indiv_grp$seg_phenotype <- as.factor(indiv_grp$seg_phenotype)

# wilcox test difference between means
wilcresult <- compare_means(COs ~ seg_phenotype,
              indiv_grp,
              method = "wilcox.test",
              p.adjust.method = "holm",
              group.by = "chromosome")

write.csv(wilcresult,
          "wilcox_means_by_chr.csv",
          row.names = F)

# t test difference between means
tresult <- compare_means(COs ~ seg_phenotype,
              indiv_grp,
              method = "t.test",
              p.adjust.method = "holm",
              group.by = "chromosome")

write.csv(tresult,
          "ttest_means_by_chr.csv",
          row.names = F)

sigdat <- data.frame(x=c(0.875,
                         1.875,
                         2.875,
                         3.875,
                         4.875,
                         5.875,
                         6.875),
                     xend=c(1.125,
                            2.125,
                            3.125,
                            4.125,
                            5.125,
                            6.125,
                            7.125),
                     y=c(7,
                         7,
                         7,
                         7,
                         7,
                         7,
                         7),
                     annotation=tresult$p.signif)

colors <- c("wt" = IBM[1],
            "ss" = IBM[2])

ggplot(indiv_grp, aes(chromosome, COs)) +
  geom_boxplot(aes(fill = phenotype,
                   colour = phenotype),
               alpha = 0.5) +
  theme(legend.text = element_text(face = "italic",
                                   size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  geom_signif(xmin = sigdat$x,
              xmax = sigdat$xend,
              y_position = sigdat$y,
              annotations = c("ns",
                              "****",
                              "ns",
                              "**",
                              "****",
                              "*",
                              "ns"),
              tip_length = 0) +
  ylim(0,8)

ggsave("crossover_boxplot.svg",
       width = 183,
       height = 120,
       units = "mm")

# sanity check sums
totals <- indiv_grp %>% group_by(phenotype, chromosome) %>%
  summarise(COs = sum(COs, na.rm = T))


length(grep("_ss", names(co_count_final)))

totals_wt <- subset(totals,
                    phenotype == "wt")
totals_ss <- subset(totals,
                    phenotype == "ss")

totals_wt$avg <- totals_wt$COs/length(grep("_wt", names(co_count_final)))
totals_ss$avg <- totals_ss$COs/length(grep("_ss", names(co_count_final)))

unique(indiv_grp$fam)
