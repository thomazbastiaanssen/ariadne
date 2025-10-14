# library(curatedMetagenomicData)
# library(tidyverse)
# library(SummarizedExperiment)
# library(MultiFactor)
# library(Matrix)
#
#
# # prep
# # Download ko2uniref LinkMap from repo
# ko2uniref90 <- importMapping(
#     ariadne("ChocoPhlAn"), ko ~ uniref90,
#     dry.run = FALSE
# )
# map <- importModules("GBM")
# map <- c(map, ko2uniref90)
#
# id <- "2021-10-14.GopalakrishnanV_2018.gene_families"
#
# x_SE <- curatedMetagenomicData::curatedMetagenomicData(id, dryrun = FALSE)[[1L]]
# x <- assay(x_SE)[grep(".s__", row.names(x_SE), fixed = TRUE),]
#
#
# # subset for speed
# x <- x[grep("s__Bacteroides", row.names(x), fixed = TRUE),]
#
#
# feat2sub(row.names(sub2samp_mat), "|") |> str()
#
#
# colnames(feat2sub) <- c("uniref90", "subtype")
#
# LinkMap(feat2sub)
#
# dim(as.matrix(LinkMap(feat2sub)))
#
# row.names(sub2samp_mat) <- feat2sub[[2L]]
#
# sub2samp_mat %>%
#     Matrix(sparse = FALSE) %>%
#     as.matrix() %>%
#     as.data.frame(row.names = FALSE) %>%
#     mutate(subtype = feat2sub$subtype, .before = 1) %>%
#     mutate(feature = feat2sub$uniref90, .before = 2 ) %>%
#     pivot_longer(!c(subtype, feature)) %>%
#     group_by(subtype, name) %>%
#
#     reframe(nzero = sum(value == 0L),
#             nonzero = sum(value != 0L)
#             ) %>%
#     ungroup()
#
#
#
# s2s2 <- Matrix::which(sub2samp_mat != 0L, arr.ind = TRUE, useNames = TRUE)
#
# subtype <- s2s2[,1L]
# samples <- s2s2[,2L]
# attr(subtype, "levels") <- row.names(sub2samp_mat)
# class(subtype) <- "factor"
# attr(samples, "levels") <- colnames(sub2samp_mat)
# class(samples) <- "factor"
#
# sub2samp <- LinkMap(data.frame(subtype, samples))
#
#
#
# m <- MultiFactor(list(sub2samp, feat2sub, ko2uniref90))
#
#
# m
