network.stats=function(ifn.map.select, fdir)
{
    df = load.table(ifn.map.select)
    df = df[df$type != "stable"& df$type != "unknown",]
    t.edges = table(df$type)

    system(paste("mkdir -p", fdir))

    save.table(data.frame(type=names(t.edges), count=as.vector(t.edges)), paste(fdir, "/edge_change_stats.txt", sep=""))

    lines = NULL
    lines = c(lines, sprintf("number of elements: %d", length(unique(df$cluster))))
    lines = c(lines, sprintf("number of anchors: %d", length(unique(df$anchor))))
    ofn = paste(fdir, "/network_compare_stats.txt", sep="")
    cat(sprintf("generating file: %s\n", ofn))
    fc = file(ofn)
    writeLines(lines, fc)
    close(fc)
}
