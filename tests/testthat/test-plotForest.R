test_that("plotForest", {

tse <- makeTSE(nrow = 10L, ncol = 10L)

est <- rnorm(nrow(tse), mean = 0)
conf.lower <- est - 0.1 * rpois(length(est), lambda = 3)
conf.upper <- est + 0.1 * rpois(length(est), lambda = 3)
p.val <- runif(est)
extra <- sample(100, length(est))

df <- data.frame(
  effect = est,
  lower = conf.lower,
  upper = conf.upper,
  pval = p.val,
  extra = extra
)

rowData(tse) <- cbind(rowData(tse), df)
colData(tse) <- cbind(colData(tse), df)

p <- plotForest(tse, label.by = c("CI", "P-Value", "extra"))
expect_s3_class(p, "ggplot")

colTree(tse) <- NULL

expect_warning(
  p <- plotForest(tse, by = 2),
  "'show.tree' is ignored when row/colTree(x) does not exist.",
  fixed = TRUE
)

expect_s3_class(p, "ggplot")

p <- plotForest(df, colour.by = "extra")
expect_s3_class(p, "ggplot")

plotForest(
  tse,
  label.by = c("CI", "P-Value", "extra"),
  show.tree = TRUE,
  colour.by = "var1",
  edge.colour.by = "var2"
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

})