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
