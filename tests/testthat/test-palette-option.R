test_that("plots follow the ggplot2 discrete palette options", {
    data(GlobalPatterns, package = "mia")
    tse <- GlobalPatterns[, GlobalPatterns$SampleType %in% c(
        "Feces", "Skin", "Tongue", "Ocean")]
    tse <- mia::transformAssay(tse, method = "relabundance")
    tse <- mia::agglomerateByRank(tse, rank = "Phylum")
    tse <- tse[mia::getTop(tse, 4), ]
    pal <- unname(palette.colors(palette = "Okabe-Ito"))[-1]
    fills <- function(p) unique(ggplot2::ggplot_build(p)$data[[1]]$fill)
    # Default: Tableau palette
    old <- options(
        ggplot2.discrete.colour = NULL, ggplot2.discrete.fill = NULL)
    on.exit(options(old))
    p <- plotAbundance(tse, assay.type = "relabundance")
    expect_true(all(fills(p) %in% scater:::.get_palette("tableau10medium")))
    # Character vector option
    options(ggplot2.discrete.colour = pal, ggplot2.discrete.fill = pal)
    p <- plotAbundance(tse, assay.type = "relabundance")
    expect_true(all(fills(p) %in% pal))
    # List of palettes: the first with enough colours is used
    options(ggplot2.discrete.fill = list(pal[1:2], rev(pal)))
    p <- plotAbundance(tse, assay.type = "relabundance")
    expect_true(all(fills(p) %in% pal))
    expect_false(all(fills(p) %in% pal[1:2]))
    # Too few colours: fall back to the default
    options(ggplot2.discrete.fill = pal[1:2])
    p <- plotAbundance(tse, assay.type = "relabundance")
    expect_true(all(fills(p) %in% scater:::.get_palette("tableau10medium")))
    # Helper
    options(ggplot2.discrete.colour = pal)
    expect_identical(miaViz:::.get_option_palette(FALSE, 3), pal)
    expect_null(miaViz:::.get_option_palette(FALSE, 20))
    options(ggplot2.discrete.colour = function(...) NULL)
    expect_null(miaViz:::.get_option_palette(FALSE, 3))
})

test_that("palette option counts only the groups present in the data", {
    data(GlobalPatterns, package = "mia")
    # SampleType keeps all 9 levels after subsetting, but only 3 are present
    tse <- GlobalPatterns[, GlobalPatterns$SampleType %in% c(
        "Feces", "Skin", "Tongue")]
    tse <- mia::addAlpha(tse, index = "shannon")
    pal <- unname(palette.colors(palette = "Okabe-Ito"))[-1]
    old <- options(ggplot2.discrete.colour = pal, ggplot2.discrete.fill = pal)
    on.exit(options(old))
    p <- plotBoxplot(
        tse, col.var = "shannon", x = "SampleType", fill.by = "SampleType")
    fills <- unique(ggplot2::ggplot_build(p)$data[[1]]$fill)
    expect_true(all(fills %in% pal))
})
