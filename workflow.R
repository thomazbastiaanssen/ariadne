
graph <- ariadne(Resources)

plot(graph)

search(graph, refseq ~ rhea, k = 5)

plot(graph, refseq ~ rhea, k = 2)

linkmap <- weave2(graph, refseq ~ rhea, k = 2)

linkmap <- weave2(graph, ko ~ eggnog, k = 2)


my_fun <- function(by){
  mf <- randomMultiFactor()
  linkmap <- weave(mf, by)
  return(linkmap)
}

linkmap <- my_fun(a ~ b)
