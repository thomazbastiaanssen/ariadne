devtools::load_all("../MultiFactor")
devtools::load_all()

graph <- ariadne()

plotPath(graph, gbm ~ uniref90 ~ taxid, focus = TRUE)

linkmaps <- weavePath(graph, gbm ~ uniref50, k = 3)

mff <- MultiFactor(lapply(linkmaps, function(x) head(x, 100000)))

lmp <- MultiFactor:::.stack_by_formula(mff, gbm ~ ko + eggnog + tigr)


