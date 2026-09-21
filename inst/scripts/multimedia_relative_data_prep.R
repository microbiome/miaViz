# -----------------------------------------------------------------------------
# multimedia_relative_data_prep.R
#
# Preparation script for the miaViz example datasets `se_relative` and
# `tse_relative` (see `?miaViz::multimedia-datasets`).
#
# Source: curatedMetagenomicData -> HMP_2019_ibdmdb.relative_abundance
#         (Lloyd-Price et al., 2019, Nature 569:655-662).
#
# Pipeline:
#   1. Load the iHMP IBDMDB relative-abundance TreeSummarizedExperiment.
#   2. Convert the raw assay to relative abundances in [0, 1].
#   3. Compute a per-sample dysbiosis score as the median Bray-Curtis
#      distance to the healthy reference set (excluding self).
#   4. Subset to IBD subjects with paired baseline (visit 1) and
#      post-treatment (visit 21) samples.
#   5. Define a `treatment` factor (T0 = baseline, T1 = post-treatment).
#   6. Standardize the relative-abundance assay as `scaled`.
#   7. Save `tse_relative` (TreeSummarizedExperiment) and `se_relative`
#      (SummarizedExperiment) as separate .rda files under `data/`.
#
# Run this script from the package root, or set the destination path
# in the save() calls at the end accordingly.
# -----------------------------------------------------------------------------

remove(list = ls())
invisible(gc())

##################
# Load libraries #
##################

library(curatedMetagenomicData)
library(SummarizedExperiment)
library(dplyr)
library(vegan)
library(tidyverse)
library(multimedia)
library(stringr)
library(ggplot2)
library(mia)

##################
# Load iHMP data #
##################

# Load data
tse <- curatedMetagenomicData(
    "HMP_2019_ibdmdb.relative_abundance",
    rownames = "short",
    dryrun = FALSE
)[[1]]

# Assign SampleID for matching
tse$SampleID <- colnames(tse)

# Convert relative_abundance assay to relabundance (which is in [0,1] interval)
tse <- transformAssay(
    tse, assay.type = "relative_abundance", method = "relabundance"
)

########################
# Reference nonIBD set #
########################

# Reference set: healthy
tse$disease_binary <- tse$disease == "healthy"

########################################
# Calculate Bray-Curtis dissimilarity  #
########################################

diss <- as.matrix(getDissimilarity(
    tse, method = "bray", na.rm = TRUE, assay.type = "relabundance"
))

#################################
# Calculate the dysbiosis score #
#################################

# For each sample i, compute the median distance to all healthy reference
# samples (excluding sample i itself).
sample_ids <- tse$SampleID

tse$dysbiosis <- sapply(seq_along(tse$disease_binary), function(i) {
    ref_others <- tse$disease_binary & (sample_ids != sample_ids[i])
    median(diss[i, ref_others], na.rm = TRUE)
})

#############################################################
# Subset to IBD only and only one post-baseline time point  #
#############################################################

# Keep IBD only
tse <- tse[, tse$disease != "healthy"]

# Visit 1 (baseline) or 21 (post-treatment)
tse <- tse[, tse$visit_number %in% c(1, 21)]

# Keep subjects with both visits
keep_subjects <- names(which(table(tse$subject_id) == 2))
tse <- tse[, tse$subject_id %in% keep_subjects]

# Define time point label
tse$Time_point <- ifelse(tse$visit_number == 1, "T0", "T1")

# Define treatment: Pre and Post (Baseline vs. Post-baseline)
tse$treatment <- factor(tse$Time_point, levels = c("T0", "T1"))

#############################
# Standardize the mediators #
#############################

tse <- transformAssay(
    tse, assay.type = "relabundance", method = "standardize", name = "scaled"
)

########################################
# Build SummarizedExperiment companion #
########################################

# Convert the TreeSummarizedExperiment to be SummarizedExperiment
se_relative <- as(tse, "SummarizedExperiment")
if (is.null(rownames(se_relative))) {
    rownames(se_relative) <- rownames(rowData(se_relative)) <-
        make.names(rownames(tse), unique = TRUE)
}

tse_relative <- tse

################
# Save as .rda #
################

# One object per file, matching the layout of the other miaViz datasets
# (col_graph.rda, row_graph.rda, ...).
save(se_relative,  file = "data/se_relative.rda",  compress = "xz")
save(tse_relative, file = "data/tse_relative.rda", compress = "xz")
