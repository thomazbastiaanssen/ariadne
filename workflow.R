devtools::load_all()

graph <- ariadne(Resources)

plotPath(graph, gbm ~ uniref50, k = 3, rm.empty = TRUE)

linkmaps <- weavePath(graph, gbm ~ uniref50, k = 3)

mff <- MultiFactor(lapply(linkmaps, function(x) head(x, 100000)))

lmp <- .stack_by_formula(mff, gbm ~ ko + eggnog + tigr)


