# Sat Mar  2 11:34:47 2019 ------------------------------
# limma model, Venn diagrams, permutation and std genotypic coefficients
library(limma)
library(pheatmap)
library(xlsx)
library(ellipse)
library(RColorBrewer)

load("data/preprocessed.RData")

rownames(metab$LC) <- gsub(rownames(metab$LC), pattern = "_NA_",
                           replacement = "_ctr_")
rownames(metab$GC) <- gsub(rownames(metab$GC), pattern = "_NA_",
                           replacement = "_ctr_")
rownames(metab$vol) <- gsub(rownames(metab$vol), pattern = "_NA_",
                            replacement = "_ctr_")
# Standardise all
metab <- lapply(metab, scale)

# Define explanatory vars from experimental factors
key <- strsplit(rownames(metab$LC), split = "_")

annRow <- data.frame("G" = factor(sapply(key, "[[", 1)),
                     "Atm" = factor(sapply(key, "[[", 2)),
                     "Time" = as.numeric(sapply(key, "[[", 3)),
                     row.names = rownames(metab$LC))
annRow$Atm <- factor(annRow$Atm, levels = c("ctr", "03", "CO2"))

# Build design for limma - compare additive to w/ and w/o GxE and GxExT interactions
design1 <- model.matrix(~ 0 + G + Atm + Time, data = annRow)
design2 <- model.matrix(~ 0 + G*Atm + Time, data = annRow)
design3 <- model.matrix(~ 0 + G*Atm*Time, data = annRow)
modSel <- lapply(metab, function(m){
      X <- t(m)
      selectModel(X, designlist = list(GpE = design1,
                                       GxE = design2,
                                       GxExT = design3),
                  crit = "bic")
})
modSelStats <- sapply(modSel, function(x){table(x$pref)})
# additive models are overall better
apply(modSelStats, 2, function(c){c/sum(c)}) # supptable

# Contrasts to be used later
cntr <- makeContrasts(Atm03 = "Atm03",
                      AtmCO2 = "AtmCO2",
                      Time = "Time", levels = design1)

# Run limma with lapply
results <- lapply(metab, function(x){
      X <- t(x)
      fit <- lmFit(X, design = design1)
      return(fit)
})

# Venn/Euler diagrams
par(mfrow = c(3,1))
vennRes <- lapply(results, function(X){
   fit <- contrasts.fit(X, contrasts = cntr)
   condtest <- decideTests(fit, method = "global",
                           adjust.method = "BH", p.value = 0.01)
   vennDiagram(condtest, circle.col = c("red","green","blue"),
               include = c("up", "down"), counts.col = 1:2)
   return(condtest)
})

# Write out model decideTests
write.xlsx(vennRes$LC,
           file = "data/decideTests.xlsx",
           sheetName = "LC", append = F)
write.xlsx(vennRes$GC,
           file = "data/decideTests.xlsx",
           sheetName = "GC", append = T)
write.xlsx(vennRes$vol,
           file = "data/decideTests.xlsx",
           sheetName = "Vol", append = T)

# Pick the relevant coefficients and marginal F-values
varGenos <- lapply(results, function(k){
      # F-val associated with genotypic effects
      fit2 <- eBayes(k) # makes no difference using ebayes
      tmp <- topTable(fit2, coef = 1:5, number = Inf, sort.by = "none")$F
      names(tmp) <- rownames(k$coefficients)
      return(tmp)
})

# Write out model coefficients + F-vals
lcAll <- data.frame(results$LC$coefficients, Fval = varGenos$LC)
write.xlsx(lcAll,
           file = "data/limmaCoefs.xlsx",
           sheetName = "LC", append = F)
gcAll <- data.frame(results$GC$coefficients, Fval = varGenos$GC)
write.xlsx(gcAll,
           file = "data/limmaCoefs.xlsx",
           sheetName = "GC", append = T)
volAll <- data.frame(results$vol$coefficients, Fval = varGenos$vol)
write.xlsx(volAll,
           file = "data/limmaCoefs.xlsx",
           sheetName = "Vol", append = T)

# Time, O3 and CO2 effects
timeE <- lapply(results, function(x){
   fit <- contrasts.fit(x, contrasts = cntr)
   fit$coefficients[,"Time"]
})

atmO3 <- lapply(results, function(x){
      fit <- contrasts.fit(x, contrasts = cntr)
      fit$coefficients[,"Atm03"]
})
atmCO2 <- lapply(results, function(x){
      fit <- contrasts.fit(x, contrasts = cntr)
      fit$coefficients[,"AtmCO2"]
})

# Create ellipse-drawing function for 95% CI based on simulation
permFun <- function(n){
      # Simulate matrix with random vars following standard normal
      rand <- scale(matrix(rnorm(n*nrow(annRow)), nrow = nrow(annRow),
                     ncol = n, byrow = F))
      # Do same as above
      Xperm <- t(rand)
      fitperm <- contrasts.fit(lmFit(Xperm, design = design1), contrasts = cntr)
      # Extract coefficients for genotype-level F value
      # randF <- topTable(fitperm, coef = 1:5, number = Inf, sort.by = "none")$F
   
      # Atm coefficients
      o3 <- fitperm$coefficients[, "Atm03"]
      co2 <- fitperm$coefficients[, "AtmCO2"]
      
      return(data.frame(pO3 = o3, pCO2 = co2))
}

# Plot all
par(mar = c(6, 6, 3, 3), mfrow = c(3,1))
lmt <- .15
maxTime <- max(abs(unlist(timeE)))
## LC ##
idx <- atmO3$LC**2 + atmCO2$LC**2 > lmt
plot(atmO3$LC, atmCO2$LC, abline(v = 0, h = 0), main = "LC",
     ylab = expression(beta[CO[2]]),
     xlab = expression(beta[O[3]]), pch = 16,
     col = rgb(ifelse(vennRes$LC[,"Time"] == 1, 1, 0), 0,
               ifelse(vennRes$LC[,"Time"] == -1, 1, 0),
               ifelse(vennRes$LC[,"Time"] == 0, .25, .5)),
     cex = log(varGenos$LC) / 3)
# Permutation + ellipse 95% CI
permData <- permFun(1e5)
for(lev in c(.9, .95, .99)){
      lines(ellipse(cov(permData), level = lev),
            col = rgb(0,0,0,.5), lty = 2)
}
# Labels
text(atmO3$LC[idx], atmCO2$LC[idx], labels = names(varGenos$LC)[idx], pos = 3, cex = .75)

legend("bottomright", pch = 16, pt.cex = log(c(10, 100, 1000)) / 3,
       col = rgb(0,0,0,.5), bty = "n",
       legend = c(10, 10, 100), title = "Gen Var F-value")

## GC ##
idx <- atmO3$GC**2 + atmCO2$GC**2 > lmt
plot(atmO3$GC, atmCO2$GC, abline(v = 0, h = 0), main = "GC",
     ylab = expression(beta[CO[2]]),
     xlab = expression(beta[O[3]]), pch = 16,
     col = rgb(ifelse(vennRes$GC[,"Time"] == 1, 1, 0), 0,
               ifelse(vennRes$GC[,"Time"] == -1, 1, 0),
               ifelse(vennRes$GC[,"Time"] == 0, .25, .5)),
     cex = log(varGenos$GC) / 3)
# Permutation + ellipse 95% CI
for(lev in c(.9, .95, .99)){
   lines(ellipse(cov(permData), level = lev),
         col = rgb(0,0,0,.5), lty = 2)
}
# Labels
text(atmO3$GC[idx], atmCO2$GC[idx], labels = names(varGenos$GC)[idx], pos = 3, cex = .75)

## VOL ##
idx <- atmO3$vol**2 + atmCO2$vol**2 > lmt
plot(atmO3$vol, atmCO2$vol, abline(v = 0, h = 0), main = "vol",
     ylab = expression(beta[CO[2]]),
     xlab = expression(beta[O[3]]), pch = 16,
     col = rgb(ifelse(vennRes$vol[,"Time"] == 1, 1, 0), 0,
               ifelse(vennRes$vol[,"Time"] == -1, 1, 0),
               ifelse(vennRes$vol[,"Time"] == 0, .25, .5)),
     cex = log(varGenos$vol) / 3)
# Permutation + ellipse 95% CI
for(lev in c(.9, .95, .99)){
   lines(ellipse(cov(permData), level = lev),
         col = rgb(0,0,0,.5), lty = 2)
}
# Labels
text(atmO3$vol[idx], atmCO2$vol[idx], labels = names(varGenos$vol)[idx], pos = 3, cex = .75)

# Heatmap with standardised genotypic coefficients
finalD <- list(GC = gcAll,
               LC = lcAll,
               Vol = volAll)

lapply(finalD, function(m){
      fvalMat <- log(m[,"Fval",drop = F])
      colnames(fvalMat)[1] <- "logF"
      pheatmap(m[,1:5], scale = "row",
               border_color = NA, labels_col = gsub(colnames(m)[1:5],
                          pat = "^G", rep = ""),
               annotation_row = fvalMat)
      return(NULL)
})
# Write session info w/ dependencies
writeLines(capture.output(sessionInfo()), "sessionInfo")
