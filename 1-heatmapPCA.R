# Load libraries and data
library(pheatmap)
load("data/preprocessed.RData")

rownames(metab$LC) <- gsub(rownames(metab$LC), pattern = "_NA_",
                           replacement = "_ctr_")
rownames(metab$GC) <- gsub(rownames(metab$GC), pattern = "_NA_",
                           replacement = "_ctr_")
rownames(metab$vol) <- gsub(rownames(metab$vol), pattern = "_NA_",
                           replacement = "_ctr_")

# Identify replicates per sample
descriptors <- strsplit(rownames(metab$LC), split = "_")
uniSamples <- unique(paste(sapply(descriptors, "[[", 1),
                    sapply(descriptors, "[[", 2),
                    sapply(descriptors, "[[", 3), sep = "_"))
# Create function to compute median values per sample
medianVal <- function(x, samples){
      res <- matrix(NA, length(samples), ncol(x))
      dimnames(res) <- list(samples, colnames(x))
      for(i in 1:length(samples)){
            idx <- grep(rownames(x), pattern = paste0("^", samples[i]))
            res[i,] <- apply(x[idx,,drop=F], 2, median, na.rm = T)
      }
      return(scale(res)) # add scaling
}

# Create list with median values per sample
medianMet <- lapply(metab, medianVal, samples = uniSamples)

# Convert experimental factors to factor type
descriptorsMed <- strsplit(rownames(medianMet$LC), split = "_")
geno <- factor(sapply(descriptorsMed, "[[", 1))
treat <- factor(sapply(descriptorsMed, "[[", 2))
time <- as.numeric(sapply(descriptorsMed, "[[", 3))

### HEATMAP ###
# Combine datasets
allDat <- t(cbind(medianMet$LC, medianMet$GC, medianMet$vol))
type <- data.frame("Type" = rep(c("LC", "GC", "Vol"),
                                sapply(medianMet, ncol)),
                   row.names = unlist(sapply(medianMet, colnames)))

# Visualize
pheatmap(allDat, border_color = NA,
         annotation_col = data.frame("Genotype" = geno,
                                     "Treatment" = treat,
                                     "Time" = time,
                                     row.names = colnames(allDat)), scale = "none",
         show_rownames = F, main = "All platforms", annotation_row = type,
         show_colnames = F, breaks = seq(-4, 4, length.out = 100))
# Idem, for the three separate platforms
pheatmap(t(medianMet$LC), border_color = NA,
         annotation_col = data.frame("Genotype" = geno,
                                     "Treatment" = treat,
                                     "Time" = time,
                                     row.names = colnames(allDat)), scale = "none",
         show_colnames = F, main = "LC",
         breaks = seq(-4, 4, length.out = 100))
pheatmap(t(medianMet$GC), border_color = NA,
         annotation_col = data.frame("Genotype" = geno,
                                     "Treatment" = treat,
                                     "Time" = time,
                                     row.names = colnames(allDat)), scale = "none",
         show_colnames = F, main = "GC",
         breaks = seq(-4, 4, length.out = 100))
pheatmap(t(medianMet$vol), border_color = NA,
         annotation_col = data.frame("Genotype" = geno,
                                     "Treatment" = treat,
                                     "Time" = time,
                                     row.names = colnames(allDat)), scale = "none",
         show_colnames = F, main = "Volatiles",
         breaks = seq(-4, 4, length.out = 100))

### PCA ###
# Display genotypes with symbols
keyGen <- c("A" = "orange", "C" = "red", "CG" = "green", "F" = "blue", "SC" = "magenta")
col <- as.data.frame(t(col2rgb(keyGen[geno], alpha = F)) / 255)

# Display treatments with different colors
keyTreat <- c("03" = 15, "CO2" = 16, "ctr" = 17)
pch <- keyTreat[treat]

# Add transparency based on time points
alpha <- (time + 2)/12

set.seed(100)
myPCA <- lapply(medianMet, prcomp)
pairs(myPCA$GC$x[,1:3], col = rgb(col$red, col$green, col$blue, alpha = alpha),
     pch = pch, main = "GC", upper.panel = NULL,
     labels = paste0("PC", 1:3, "\n",
                     round(summary(myPCA$GC)$importance[2,1:3] * 100, 2),
                     "% var exp"))
legend("top", col = c(as.character(unlist(keyGen))),#, rep("black", 3)),
       pch = rep(16, 5), legend = names(keyGen),
       bty = "n", cex = .6, xpd = T)
legend("topright", col = rep("black", 3), pch = c(15, 16, 17),
       legend = c("O3", "CO2", "ctr"), bty = "n", cex = .6, xpd = T)
pairs(myPCA$LC$x[,1:3], col = rgb(col$red, col$green, col$blue, alpha = alpha),
      pch = pch, main = "LC", upper.panel = NULL,
      labels = paste0("PC", 1:3, "\n",
                      round(summary(myPCA$LC)$importance[2,1:3] * 100, 2),
                      "% var exp"))
pairs(myPCA$vol$x[,1:3], col = rgb(col$red, col$green, col$blue, alpha = alpha),
      pch = pch, main = "Volatiles", upper.panel = NULL,
      labels = paste0("PC", 1:3, "\n",
                      round(summary(myPCA$vol)$importance[2,1:3] * 100, 2),
                      "% var exp"))