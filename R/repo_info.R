# # library(RCurl)
# # library(XML)
# # library(tidyverse)
#
# # Web of Life
# within_links <- paste0("https://ftp.microbio.me/pub/wol2/function/",
#                 c("eggnog/", "go/", "uniref/" ), "idmaps/")
#
# wol_table <- RCurl::getURL(within_links, ftp.use.epsv = TRUE, dirlistonly = TRUE) |>
#     lapply(XML::getHTMLLinks) |>
#     stack() |>
#     filter(str_detect(values, ".map")) |>
#     filter(str_detect(values, ".md5", negate = TRUE)) |>
#     mutate(link = paste0(ind, values),
#            from = str_remove_all(ind, ".*function/|/idmaps/"),
#            to   = str_remove(values, ".map.*")) |>
#     dplyr::select(from, to, link)
#
# WoL <- .reftableToDFList(wol_table)
#
# ChocoPhlAn <- list(
#     eggnog2uniref50    = data.frame(eggnog   = factor(), uniref50 = factor()),
#     eggnog2uniref90    = data.frame(eggnog   = factor(), uniref90 = factor()),
#     go2uniref50        = data.frame(go       = factor(), uniref50 = factor()),
#     go2uniref90        = data.frame(go       = factor(), uniref90 = factor()),
#     ko2uniref50        = data.frame(ko       = factor(), uniref50 = factor()),
#     ko2uniref90        = data.frame(ko       = factor(), uniref90 = factor()),
#     level4ec2uniref50  = data.frame(ec       = factor(), uniref50 = factor()),
#     level4ec2uniref90  = data.frame(ec       = factor(), uniref90 = factor()),
#     pfam2uniref50      = data.frame(pfam     = factor(), uniref50 = factor()),
#     pfam2uniref90      = data.frame(pfam     = factor(), uniref90 = factor()),
#     uniref502uniref90  = data.frame(uniref50 = factor(), uniref90 = factor())
# )
#
#
# MultiFactor::MultiFactor(lapply(WoL, LinkMapDB, repo = "WoL"))
# MultiFactor::MultiFactor(lapply(ChocoPhlAn, LinkMapDB, repo = "ChocoPhlAn"))
#
