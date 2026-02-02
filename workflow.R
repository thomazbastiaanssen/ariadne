
graph <- ariadne(Resources)

plotPath(graph)

searchPath(graph, refseq ~ rhea, k = 5)

plotPath(graph, refseq ~ rhea, k = 2)

# takes long
# linkmap <- weavePath(graph, refseq ~ rhea, k = 2)

# takes less long
linkmap <- weavePath(graph, ko ~ eggnog, k = 2)

