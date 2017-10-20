network.stats=function(ifn, fdir)
{
    network = load.table(ifn)

    system(paste("mkdir -p", fdir))
    lines = NULL
    lines = c(lines, sprintf("number of anchors: %d", length(unique(network$anchor))))
    lines = c(lines, sprintf("number of elements: %d", length(unique(network$cluster))))
    lines = c(lines, sprintf("number of edges: %d", dim(network)[1]))
    ofn = paste(fdir, "/network_stats.txt", sep="")
    cat(sprintf("generating file: %s\n", ofn))
    fc = file(ofn)
    writeLines(lines, fc)
    close(fc)
}
