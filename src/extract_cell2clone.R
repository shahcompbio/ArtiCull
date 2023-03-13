args = commandArgs(trailingOnly=TRUE)
filename = args[1]
output = args[2]
x <- readRDS(filename)
clusters <- x$cl$clustering
write.table(clusters, output, sep="\t")
