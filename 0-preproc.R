# Load libraries
library(readxl)
library(stringr)
library(pcaMethods)

# Load replicated data
GC <- data.frame(read_xlsx("Supplementary Table S1.xlsx",
                           sheet = 2))
LC <- data.frame(read_xlsx("Supplementary Table S2.xlsx",
                           sheet = 2))
vol <- data.frame(read_xlsx("Supplementary Table S3.xlsx",
                            sheet = 2))

# Log10 transform LC and volatiles, set null values to missing
LC[LC == 0] <- NA
GC[GC == 0] <- NA
vol[vol == 0] <- NA

LC[,-1] <- log10(LC[,-1])
GC[,-1] <- log10(GC[,-1])
vol[,-1] <- log10(vol[,-1])

# Make strings operations so that we can extract replicates
LC$ID <- str_replace_all(LC$ID,
                         pattern = "_B2", replacement = "")
LC$ID <- str_replace_all(LC$ID,
                         pattern = "T0", replacement = "NA_0") 

GC$ID <- str_replace_all(GC$ID,
                         pattern = "_B2", replacement = "")
GC$ID <- str_replace_all(GC$ID,
                         pattern = "T0", replacement = "NA_0") 

vol$ID <- str_replace_all(vol$ID,
                          pattern = " ", replacement = "_") 

vol$ID <- str_replace_all(vol$ID,
                          pattern = "_B2", replacement = "") 
vol$ID <- str_replace_all(vol$ID,
                          pattern = "CTRL", replacement = "ctr") 
vol$ID <- str_replace_all(vol$ID,
                          pattern = "_r", replacement = "_") 

# The volatile sample IDs require a bit of extra work
vol$ID <- gsub(pattern = "T0", replacement = "0D_NA", vol$ID)
vol$ID <- gsub(pattern = "O3", replacement = "03", vol$ID)
vol$ID <- gsub(pattern = "D_", replacement = "_", vol$ID)
struc <- strsplit(vol$ID, split = "_", fixed = T)
vol$ID <- paste(sep = "_", sapply(struc, "[[", 1),
                sapply(struc, "[[", 3),
                sapply(struc, "[[", 2),
                sapply(struc, "[[", 4))

# Verify that sample names are now common across the three platforms
allNames <- Reduce(intersect, list(LC$ID, GC$ID, vol$ID))

## Standardise, remove ID, add row names

# Write function to normalize to sample intensities
normSample <- function(dataset){
      # Median centering per sample
      return(t(apply(dataset, 1, function(x){x - median(x, na.rm = T)})))
}

LC <- data.frame(normSample(LC[match(allNames, LC$ID),-1]), row.names = allNames)
GC <- data.frame(normSample(GC[match(allNames, GC$ID),-1]), row.names = allNames)
vol <- data.frame(normSample(vol[match(allNames, vol$ID),-1]), row.names = allNames)

# Create list w/ all platforms, standardise
dat <- list(LC = LC, GC = GC, vol = vol)

# Filter (NA% <= 20) and impute
filtNAs <- function(dat, cutoff = 0.2){
      idx <- apply(dat, 2, function(x){mean(is.na(x))}) > cutoff
      pc <- pca(dat[,!idx], method = "bpca", nPcs = 5, scale = "uv", center = T)
      return(pc@completeObs)
}
metab <- lapply(dat, filtNAs)

save(metab, file = "data/preprocessed.RData")

