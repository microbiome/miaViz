data("Tengeler2020", package = "mia")
tse <- Tengeler2020

tse <- agglomerateByRank(tse, rank = "Genus")

est <- rnorm(nrow(tse), mean = 0)
conf.lower <- est - 0.1 * rpois(length(est), lambda = 3)
conf.upper <- est + 0.1 * rpois(length(est), lambda = 3)
p.val <- runif(est)
prev <- sample(100, length(est))

df <- data.frame(
  effect = est,
  lower = conf.lower,
  upper = conf.upper,
  pval = p.val,
  Prevalence = prev
)

rowData(tse) <- cbind(rowData(tse), df)

est <- rnorm(ncol(tse), mean = 0)
conf.lower <- est - 0.1 * rpois(length(est), lambda = 3)
conf.upper <- est + 0.1 * rpois(length(est), lambda = 3)
p.val <- runif(est)
prev <- sample(100, length(est))

df <- data.frame(
  effect = est,
  lower = conf.lower,
  upper = conf.upper,
  pval = p.val,
  Prevalence = prev
)

colData(tse) <- cbind(colData(tse), df)


plotForest(
  tse,
  label.by = c("CI", "P-Value", "Prevalence"),
  show.tree = TRUE,
  tip.colour.by = "Phylum"
)

tse2 <- SummarizedExperiment(assays = assays(tse),
                             rowData = rowData(tse),
                             colData = colData(tse))

plotForest(tse2, show.tree = FALSE)


dfx <- as.data.frame(rowData(tse))
plotForest(dfx)

plotForest(tse, effect.var = "effect", ci.lower.var = "hi", label.by = c())
plotForest(tse, effect.var = "effect", ci.lower.var = NULL, label.by = c())



dfx <- as.data.frame(colData(tse))
plotForest(dfx, colour.by = NULL)
plotForest(dfx, colour.by = "patient_status")

dfx2 <- rbind(dfx, dfx)
dfx2$Group <- c(rep("A", nrow(dfx2) / 2), rep("B", nrow(dfx2) / 2))

plotForest(dfx2, id.var = "sample_name", colour.by = "Group")
plotForest(dfx2, id.var = "sample_name")

