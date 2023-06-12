# starting from raw data 
library(tidyverse)

# join unique ids and sample metadata from 50k batch 48
fiftyk <- read.table("48_FINAL_TRANSPOSED.txt",
                     header = TRUE,
                     na.strings = "--",
                     sep = "\t")

metadata <- read.csv("sample_metadata.csv")

colnames(metadata)
colnames(fiftyk)

fiftyk <- rename(fiftyk,
                 Unique.ID = Line.Marker)

merged <- merge(metadata,
                fiftyk,
                by = "Unique.ID")

colnames(merged[,1:20])

fiftyk <- merged[,c(6,14:43812)]

fiftyk <- rename(fiftyk,
                 line = External.ID)

# get the Barke and BW248 marker data from batches 9 and 23 
b23 <- read.table("23_FINAL.txt",
                   header = T,
                   sep = "\t")

# unique id from http://narwhal.hutton.ac.uk:9123/50k/index.pl?batch=23
bw248 <- b23[grep("23A094",
             b23$Unique.ID),]

bw248$line <- "BW248"

bw248 <- bw248[,c(44042, 2:44041)]

colnames(bw248)

b9 <- read.table("9.txt",
                header = T,
                sep = "\t")

barke <- b9[grep("Barke",
                 b9$line),]

# for some reason there are more markers in batch 23. 
# make sure all datasets have the same number of markers and that they occur in
# the same order.
matchmarkers <- colnames(barke)
bw248 <- bw248[,matchmarkers]
fiftyk <- fiftyk[,matchmarkers]

names(bw248[,700:800])
names(fiftyk[,700:800])
names(barke[,700:800])

fiftyk <- rbind(fiftyk,
                barke,
                bw248)

# Morex V3 marker positions
snp_postions <- read.table("SNPPositions_50k_on_MorexV3.txt",
                           header = T)
rownames(snp_postions) <- snp_postions$marker

# invert 50K data table
rownames(fiftyk) <- fiftyk$line
tfiftyk <- t(fiftyk)
colnames(tfiftyk)
rownames(tfiftyk)
tfiftyk <- tfiftyk[-1,]

# replace the dot(.) in imported marker names with a dash(-)
rownames(tfiftyk) <- gsub("\\.",
                          "-",
                          rownames(tfiftyk))

# merge morex v3 marker positions with 50K results 
merged <- merge(tfiftyk,
                snp_postions,
                by = 0)

row.names(merged)

colnames(merged)

# 42912 markers map to Morex v3

# get a list of the markers that couldn't be found in the V3 marker mapping
no_positions <- setdiff(rownames(tfiftyk),
                        rownames(snp_postions))

# remove the lines we don't need
colnames(merged)

want <- c("BW248",
          "Plant",
          "Barke")

merged <- merged[,c(1,651:653,grep(paste(want, collapse="|"), colnames(merged)))]

colnames(merged)

merged <- merged[,c(1:4,186:187,5:185)]

# remove monomorphic markers between parental genotypes
rownames(merged) <- merged$Row.names
merged <- merged[,-1]
colnames(merged)

# save pre-filtered table
write.csv(merged,
          "des8_50K_pre-filter.csv",
          row.names = F)

prefilter <- merged

save(prefilter,
     file = "50k-pre-filter.RData")

nomono_parents <- merged[merged[,4] != merged[,5],]

# leaves 13316 het markers markers

# remove any markers with no call for any line 
no_blank <- nomono_parents[rowSums(is.na(nomono_parents)) != ncol(nomono_parents), ]

# leaves 13316 markers

# remove markers with any missing value
nafilter <- no_blank

nafilter <- na.omit(nafilter)
# leaves 12522 markers

# save after basic filtering
write.csv(nafilter,
          "des8_50K_basic_filtering.csv",
          row.names = F)

save(nafilter,
     file = "50k-basic-filter.RData")

# -----------------------------------------------------------------------------

# change calls to indication of match to parental genotype
# 0 == Barke,
# 2 == Bw248,
# 1 == Het
numbered_allele <- nafilter

colnames(numbered_allele)

for (x in 6:186) {
  numbered_allele[,x] <- ifelse(numbered_allele[,x] == numbered_allele[,4],
                                0,
                                ifelse(numbered_allele[,x] == numbered_allele[,5],
                                       2,
                                       1))
}

# put markers in order of chromosome and then physical position 
numbered_allele <- numbered_allele[order(numbered_allele$Chromosome_v3,
                                         numbered_allele$Position_v3),]

save(numbered_allele,
     file = "markers_sorted_post_filtering.Rdata")
