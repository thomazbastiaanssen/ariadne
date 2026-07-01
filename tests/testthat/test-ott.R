
test_that("ott", {

    # Query from ncbi ids
    ncbi_ids <- c(562, 1423, 1280)
    ott_out <- .queryOTT("ncbi", "ott", ncbi_ids, 1e6)
    tax_names <- .queryOTT("ncbi", "taxname", ncbi_ids, 1e6)

    # Qeury from ott ids
    ott_ids <- c(474506, 1084928, 1090496)
    ncbi_out <- .queryOTT("ott", "ncbi", ott_ids, 1e6)
    tax_names <- .queryOTT("ott", "taxname", ott_ids, 1e6)
    
    # Query silva ids
    silva_ids <- .queryOTT("ott", "silva", ott_ids, 1e6)
    silva_ids <- .queryOTT("ncbi", "silva", ncbi_ids, 1e6)
    
    # Query from taxnames
    names <- c("s__Escherichia coli", "s__Bacillus subtilis", "s__Staphylococcus aureus")
    ott_out <- .queryOTT("taxname", "ott", names)
    silva_ids <- .queryOTT("taxname", "silva", names)

})