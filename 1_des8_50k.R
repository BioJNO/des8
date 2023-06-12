# Filter 50K data for misplaced markers and F2 recombination events and plot
# the results 
library(runner)
library(tidyverse)
library(cowplot)

load("markers_sorted_post_filtering.Rdata")

# adjust the sample names leaving a standard format: 
# [cross]_[F1]_[F2family]_[individual]_[segregating_phenotype]
names(numbered_allele)

names(numbered_allele) <- str_replace(names(numbered_allele),
                                    "BW248\\(des8.l\\)xBarke ",
                                    "")
names(numbered_allele) <- str_replace(names(numbered_allele),
                                    "BW248 x Barke ",
                                    "")
names(numbered_allele) <- str_replace(names(numbered_allele),
                                    "Ear 0",
                                    "")
names(numbered_allele) <- str_replace(names(numbered_allele),
                                      "_Plant 00",
                                      "")
names(numbered_allele) <- str_replace(names(numbered_allele),
                                    "F3_",
                                    "")

# Once run for all seven chromosomes combine the exclusion criteria -----------
excl_lines <- unique(c(chr1H_removed_lines,
                       chr2H_removed_lines,
                       chr3H_removed_lines,
                       chr4H_removed_lines,
                       chr5H_removed_lines,
                       chr6H_removed_lines,
                       chr7H_removed_lines))

co_count_all <- rbind(chr1H_CO,
                      chr2H_CO,
                      chr3H_CO,
                      chr4H_CO,
                      chr5H_CO,
                      chr6H_CO,
                      chr7H_CO)

co_count_final <- co_count_all[, !names(co_count_all) %in% excl_lines]

length(grep("_wt", names(co_count_final)))
# 82
length(grep("_ss", names(co_count_final)))
# 78

save(co_count_final,
     file = "chromosome_CO_counts.Rdata")

# pick a chromosome -----------------------------------------------------------
chr_sel <- "chr7H"

chroma <- subset(numbered_allele,
                 Chromosome_v3 == chr_sel)

# plot raw marker calls in each group in order for this chromosome ------------
colnames(chroma)
chroma[,4] <- 0
chroma[,5] <- 2
chroma_long <- pivot_longer(chroma,
                            6:186,
                            names_to = "line",
                            values_to = "call")

# keep markers in order of position on chromosome
chroma_long$marker <- factor(chroma_long$marker,
                             levels = chroma$marker)

# create columns so that the families (group by F1 and F2) can be plotted in facet grid
chroma_long <- separate(chroma_long,
                        line,
                        c("cross","F1", "F2", "individual", "seg_phenotype"),
                        sep = "_",
                        remove = F)

# Plot the F3 marker calls in order along the chromosome (x axis) by line (y axis)
# grouped by F2 parent.
raw_calls_wt <- ggplot(chroma_long[c(grep("wt", chroma_long$line)),],
                       aes(marker,
                           line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank()) +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
raw_calls_wt

raw_calls_ss <- ggplot(chroma_long[c(grep("ss", chroma_long$line)),],
                       aes(marker,
                           line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank()) +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
raw_calls_ss

# Deal with poorly mapped, or multi-mapped out of place markers ---------------
adjusted_call <- chroma

# remove parental genotypes from lines (4 = Barke, 5 = Bw248)
colnames(adjusted_call)
adjusted_call <- adjusted_call[,-c(1:5)]

# Using a sliding window function set the value of each call to the median 
# within a sliding window of 20 markers
adjusted_call <- apply(adjusted_call[,c(1:181)],
                       2,
                       function (x) runner(x,
                                           k=20,
                                           lag=1,
                                           f = function(x) {median(x)}))

adjusted_call <- as.data.frame(adjusted_call)
adjusted_call[1,] <- adjusted_call[2,]

write.csv(adjusted_call,
          paste(chr_sel,
               "sliding_window_20_calls.csv",
               sep = "_"),
               row.names = T)

# plot adjusted call markers as above -----------------------------------------
adjusted_call$marker <- chroma$marker
adjusted_call$position <- chroma$Position_v3
adjusted_call$chromosome <- chroma$Chromosome_v3

adjusted_call_long <- pivot_longer(adjusted_call,
                                   1:181,
                                   names_to = "line",
                                   values_to = "call")

# keep markers in order
adjusted_call_long$marker <- factor(adjusted_call_long$marker,
                                    levels = chroma$marker)

adjusted_call_long <- separate(adjusted_call_long,
                        line,
                        c("cross","F1", "F2", "individual", "seg_phenotype"),
                        sep = "_",
                        remove = F)

colnames(adjusted_call_long)

adjusted_calls_wt <- ggplot(adjusted_call_long[c(grep("_wt", adjusted_call_long$line)),],
                                 aes(marker,
                                     line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
adjusted_calls_wt

adjusted_calls_ss <- ggplot(adjusted_call_long[c(grep("_ss", adjusted_call_long$line)),],
                            aes(marker,
                                line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
adjusted_calls_ss

# remove outlier lines and markers --------------------------------------------
# We have some lines and some markers with absurdly variable calls which is 
# likely due to technical error. 
# 
# The following section of code identifies these based on the number of 
# transition states present for each after the sliding window function.
# These are identified as positions along the chromosome with a adjusted 
# value equal to 0.5 (a transition between 0 and 1 or vice versa) or 
# a value of 1.5 (a transition between 1 and 2 or vice versa).
colnames(adjusted_call_long)
adjusted_call_long$call <- as.factor(adjusted_call_long$call)

# first outlier markers
marker_COs <- adjusted_call_long %>%
  group_by(marker) %>%
  summarize(co_count = sum(str_count(call, '0.5|1.5')))

hist(marker_COs$co_count,
     breaks = 50,
     col = "#FFB000")

# calculate quantiles at 1% intervals
qmark <- quantile(marker_COs$co_count,
                  probs = seq(0, 1, 0.005))
qmark

# note markers with CO counts above 95% quantile for removal
qmark[["95.0%"]]

remove_mark <- subset(marker_COs,
                      co_count > qmark[["95.0%"]])

remove_mark <- remove_mark$marker 

adjusted_call_mfilt <- adjusted_call_long[!adjusted_call_long$marker %in% remove_mark,]

line_COs <- adjusted_call_mfilt %>%
  group_by(line) %>%
  summarize(co_count = sum(str_count(call, '0.5|1.5')))

hist(line_COs$co_count,
     breaks = 50,
     col = "#FFB000")

# calculate quantiles at 1% intervals
qline <- quantile(line_COs$co_count,
                  probs = seq(0, 1, 0.005))
qline

# note markers with CO counts above 98.5% quantile for removal
qline[["98.5%"]]

remove_line <- subset(line_COs,
                      co_count > qline[["98.5%"]])

remove_line <- remove_line$line

remove_line

removed_lines <- adjusted_call_mfilt[adjusted_call_mfilt$line %in% remove_line,]
removed_lines$call <- as.character(removed_lines$call)
removed_lines$call <- as.numeric(removed_lines$call)

removed_lines_plot <- ggplot(removed_lines[c(grep("_wt|_ss", removed_lines$line)),],
                            aes(marker,
                                line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
removed_lines_plot

adjusted_call_lfilt <- adjusted_call_mfilt[!adjusted_call_mfilt$line %in% remove_line,]

colnames(adjusted_call_lfilt)

adjusted_call_lfilt <- unite(adjusted_call_lfilt,
                               fam,
                               c(6:7),
                               sep = "_",
                               remove = F)

adjusted_call_lfilt$call <- as.character(adjusted_call_lfilt$call)
adjusted_call_lfilt$call <- as.numeric(adjusted_call_lfilt$call)

adjusted_calls_wt_filt <- ggplot(adjusted_call_lfilt[c(grep("_wt", adjusted_call_lfilt$line)),],
                                 aes(marker,
                                     line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
adjusted_calls_wt_filt

adjusted_calls_ss_filt <- ggplot(adjusted_call_lfilt[c(grep("_ss", adjusted_call_lfilt$line)),],
                                 aes(marker,
                                     line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
adjusted_calls_ss_filt

# if all values are homozygous, F2 parent was most likely homozygous and the 
# region is uninformative (i.e. a crossover might happen but would we couldn't
# detect it)
# 
# If *almost* all the calls are homozygous this either indicates there was a
# large number of crossovers in the same location in all F3--which is very
# unlikely--or that this individual is incorrectly assigned to this 
# family--much more likely. 
# 
# To calculate recombination frequency accurately we want to determine 
# informative and uninformative regions so we can calculate recombination
# events over an interval where we would have been able to detect them.
#
# F3 individuals incorrectly assigned to an F2 family mess this up so 
# we want to get rid of them.
#
# First let's determine the consensus position per family
fam_consensus <- adjusted_call_lfilt %>%
  # add a column with count of number of times this call occurs per F2
  add_count(fam, marker, call, name = "F2_marker_count") %>%
  # select max value for each marker and F2 to represent the majority
  group_by(fam, marker) %>%
  # keep only first
  mutate(Majority = call[F2_marker_count == max(F2_marker_count)][1])

# Add a column with the count of F2 individuals
fam_consensus <- fam_consensus %>%
  add_count(fam, marker, name = "F2_indiv_count")

# Add a column giving the ratio of the call for this individual at this marker 
# as a proportion of all calls for this marker
fam_consensus$ratio <- fam_consensus$F2_marker_count/fam_consensus$F2_indiv_count

# Subset majority homozygous positions
hom_regions <- subset(fam_consensus,
                      Majority == 0 |
                      Majority == 2)

# Make sure this is a big majority
hom_regions <- hom_regions %>%
  group_by(fam, marker) %>%
  filter(any(ratio>0.75)) 

# Plot these regions to check this worked
hom_plot_wt <- ggplot(hom_regions[c(grep("_wt", hom_regions$line)),],
                                aes(marker,
                                    line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
hom_plot_wt

hom_plot_ss <- ggplot(hom_regions[c(grep("_ss", hom_regions$line)),],
                      aes(marker,
                          line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
hom_plot_ss

hist(hom_regions$ratio,
     breaks = 100,
     col = "#FFB000")

# pull out the lines that look misplaced and remove them
remove_indiv <- subset(hom_regions,
                          ratio < 0.3)
remove_indiv <- remove_indiv$line
remove_indiv_tab <- as.data.frame(table(remove_indiv))

hist(remove_indiv_tab$Freq,
     breaks = 50,
     col = "#FFB000")

# We only want sustained outlier lines, so exclude lines with only a few markers
# out of place from filtering
remove_indiv <- subset(remove_indiv_tab,
                       Freq > 10)
remove_indiv <- remove_indiv$remove_indiv
remove_indiv <- unique(remove_indiv)

# add these to out list of filtered lines for this chromosome
remove_indiv <- as.character(remove_indiv)
remove_line <- c(remove_indiv, remove_line)

assign(paste(chr_sel,
             "removed_lines",
             sep = "_"), remove_line)

# Filter out these lines
hom_regions <- hom_regions[!hom_regions$line %in% remove_indiv,]

# plot the result of filtering
hom_plot_wt_filt <- ggplot(hom_regions[c(grep("_wt", hom_regions$line)),],
                      aes(marker,
                          line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
hom_plot_wt_filt

hom_plot_ss_filt <- ggplot(hom_regions[c(grep("_ss", hom_regions$line)),],
                      aes(marker,
                          line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
hom_plot_ss_filt

# If this looks good remove these lines from the main dataset
fam_consensus <- fam_consensus[!fam_consensus$line %in% remove_indiv,]

# Re-calculate after filtering
fam_consensus <- fam_consensus %>%
  # add a column with count of number of times this call occurs per F2
  add_count(fam, marker, call, name = "F2_marker_count") %>%
  # select max value for each marker and F2 to represent the majority
  group_by(fam, marker) %>%
  # keep only first
  mutate(Majority = call[F2_marker_count == max(F2_marker_count)][1])

# Add a column with the count of F2 individuals
fam_consensus <- fam_consensus %>%
  add_count(fam, marker, name = "F2_indiv_count")

# Add a column giving the ratio of the call for this individual at this marker 
# as a proportion of all calls for this marker
fam_consensus$ratio <- fam_consensus$F2_marker_count/fam_consensus$F2_indiv_count
# add in the 0.5 and 1.5 to exclude the transitions to F2 homologous regions from 
# crossover counting
hom_regions <- subset(fam_consensus,
                      Majority == 0 |
                      Majority == 2 |
                      Majority == 0.5 |
                      Majority == 1.5)

hom_regions <- hom_regions %>%
  group_by(fam, marker) %>%
  filter(any(ratio>0.75)) 

# set homozygous (uninformative) region calls to NA
hom_regions <- unite(hom_regions,
                     marker_line,
                     c(1,4),
                     sep = "_",
                     remove = F)

fam_consensus <- unite(fam_consensus,
                       marker_line,
                       c(1,4),
                       sep = "_",
                       remove = F)

fam_consensus[fam_consensus$marker_line %in% hom_regions$marker_line, "call"] <- NA

# plot the result of this filtering
adjusted_calls_wt_con <- ggplot(fam_consensus[c(grep("_wt", fam_consensus$line)),],
                                aes(marker,
                                    line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
adjusted_calls_wt_con

adjusted_calls_ss_con <- ggplot(fam_consensus[c(grep("_ss", fam_consensus$line)),],
                                aes(marker,
                                    line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
adjusted_calls_ss_con

plot_grid(adjusted_calls_wt_filt,
          adjusted_calls_wt_con,
          ncol = 1,
          labels = "auto")

plot_grid(adjusted_calls_ss_filt,
          adjusted_calls_ss_con,
          ncol = 1,
          labels = "auto")

# Informative regions must be heterozygous in the F2, next we want to create a
# variable that stores whether each F3 is the same or different to the F3 at 
# each marker position.
#
# First, pull out the het regions ---------------------------------------------
het_regions <- fam_consensus[!fam_consensus$marker_line %in% hom_regions$marker_line,]

het_regions <- fam_consensus %>%
  group_by(fam, marker) %>%
  filter(any(call==1)) 

# Plot them to see that this looks right
het_plot_wt <- ggplot(het_regions[c(grep("_wt", het_regions$line)),],
                           aes(marker,
                               line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
het_plot_wt

het_plot_ss <- ggplot(het_regions[c(grep("_ss", het_regions$line)),],
                      aes(marker,
                          line)) +
  geom_tile(aes(fill = call)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
het_plot_ss

# We'll use the "Majority" column to store the presumed F2 genotype at this
# position (het/1)
fam_consensus[fam_consensus$marker_line %in% het_regions$marker_line, "Majority"] <- 1

# Next we'll make a variable to store whether the F3 matches the presumed F2
# where 
fam_consensus$match <- ifelse(fam_consensus$call == fam_consensus$Majority,
                              0,
                              1)

# Plot this to see if it looks right 
wt_match <- ggplot(fam_consensus[c(grep("_wt", fam_consensus$line)),],
                                aes(marker,
                                    line)) +
  geom_tile(aes(fill = match)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free") +
  scale_fill_continuous_sequential(palette = "Heat")
wt_match

plot_grid(adjusted_calls_wt_con,
          wt_match,
          ncol = 1,
          labels = "auto")

ss_match <- ggplot(fam_consensus[c(grep("_ss", fam_consensus$line)),],
                   aes(marker,
                       line)) +
  geom_tile(aes(fill = match)) +
  theme(axis.text = element_blank())  +
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free") +
  scale_fill_continuous_sequential(palette = "Heat")
ss_match

plot_grid(adjusted_calls_ss_con,
          ss_match,
          ncol = 1,
          labels = "auto")

# Now, we can create crossover counts from this data 
#
# configure the filtered marker data so we can calculate crossovers
co_count <- pivot_wider(fam_consensus,
                        id_cols = c("marker", "chromosome", "position"),
                        values_from = "match",
                        names_from = "line")
colnames(co_count)

# To avoid the loop referring to the changes it makes in real time store a
# static version of the data frame
Co_count_static <- co_count

# With the markers in order of position on the chromosome, count a crossover if
# in a het region there is a change from a match to the F2 (1) to a mismatch (0)
#
# Only count COs where the state change is sustained over three markers in each
# direction.
for (x in 4:nrow(co_count)) {
  co_count[x,4:ncol(co_count)] <- ifelse(is.na(Co_count_static[x,4:ncol(co_count)]),
                                         NA,
                                         ifelse(Co_count_static[x-3,4:ncol(co_count)] == 1 &
                                                  Co_count_static[x-2,4:ncol(co_count)] == 1 &
                                                  Co_count_static[x-1,4:ncol(co_count)] == 1 &
                                                  Co_count_static[x,4:ncol(co_count)] == 0 &
                                                  Co_count_static[x+1,4:ncol(co_count)] == 0 & 
                                                  Co_count_static[x+2,4:ncol(co_count)] == 0 &
                                                  Co_count_static[x+3,4:ncol(co_count)] == 0 |
                                                  Co_count_static[x-3,4:ncol(co_count)] == 0 &
                                                  Co_count_static[x-2,4:ncol(co_count)] == 0 &
                                                  Co_count_static[x-1,4:ncol(co_count)] == 0 &
                                                  Co_count_static[x,4:ncol(co_count)] == 1 &
                                                  Co_count_static[x+1,4:ncol(co_count)] == 1 & 
                                                  Co_count_static[x+2,4:ncol(co_count)] == 1 &
                                                  Co_count_static[x+3,4:ncol(co_count)] == 1,
                                                1,
                                                0))
}

# set the value in the first row to 0
co_count[1:3,4:ncol(co_count)] <- 0

# write out the CO count for this chromosome
write.csv(co_count,
          paste(chr_sel,
                "cos_pre_filter_sw20.csv",
                sep = "_"),
                row.names = T)

# plot crossover counts to see that they make sense
colnames(co_count)
co_count$marker
co_count_long <- pivot_longer(co_count,
                              4:ncol(co_count),
                              names_to = "line",
                              values_to = "COs")

co_count_long <- separate(co_count_long,
                               line,
                               c("cross","F1", "F2", "individual", "seg_phenotype"),
                               sep = "_",
                               remove = F)

co_count_long$marker <- factor(co_count_long$marker,
                               levels = chroma$marker)

co_plot_wt <- ggplot(co_count_long[c(grep("_wt", co_count_long$line)), ],
                         aes(marker,
                             line)) +
  geom_tile(aes(fill = COs)) +
  theme(axis.text = element_blank()) +
  scale_fill_continuous_diverging(palette = "Berlin") + 
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
co_plot_wt

plot_grid(adjusted_calls_wt_con,
          wt_match,
          co_plot_wt,
          ncol = 1,
          labels = "auto")

co_plot_ss <- ggplot(co_count_long[c(grep("_ss", co_count_long$line)), ],
                         aes(marker,
                             line)) +
  geom_tile(aes(fill = COs)) +
  theme(axis.text = element_blank()) +
  scale_fill_continuous_diverging(palette = "Berlin") + 
  facet_grid(F1 + F2 ~ .,
             scales = "free",
             space = "free")
co_plot_ss

plot_grid(adjusted_calls_ss_con,
          ss_match,
          co_plot_ss,
          ncol = 1,
          labels = "auto")

# Plot WT CO distribution by line for this chromosome
hist(colSums(co_count[,grep("_wt",
                            names(co_count))],
             na.rm = T),
     breaks = 5,
     col = "#FFB000")

# Plot bw248 co distribution by line for this chromosome
hist(colSums(co_count[,grep("_ss",
                            names(co_count))],
             na.rm = T),
     breaks = 5,
     col = "#FFB000")

# Make and save filtering plots -----------------------------------------------
plot_grid(raw_calls_wt,
          adjusted_calls_wt_filt,
          wt_match,
          co_plot_wt,
          ncol = 1,
          labels = "auto")

ggsave(paste(chr_sel,
             "_filtering_wt_sw20.png",
             sep = "_"),
       height = 183,
       width = 240,
       dpi = 1200,
       units = "mm")

plot_grid(raw_calls_ss,
          adjusted_calls_ss_filt,
          ss_match,
          co_plot_ss,
          ncol = 1,
          labels = "auto")

ggsave(paste(chr_sel,
             "_filtering_ss_sw20.png",
             sep = "_"),
       height = 183,
       width = 240,
       dpi = 1200,
       units = "mm")

# Save chromosome CO counts ---------------------------------------------------
assign(paste(chr_sel,
             "CO",
             sep = "_"), co_count)
