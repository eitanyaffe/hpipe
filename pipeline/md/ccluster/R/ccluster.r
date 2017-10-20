plot.marginals=function(contigs.ifn, contacts.ifn, threshold, fdir)
{
    options(stringsAsFactors=F)
    contigs = read.delim(contigs.ifn)
    contacts = read.delim(contacts.ifn)
    contacts = contacts[contacts$contig1 != contacts$contig2,]
    t = table(contacts$contig1)
    system(paste("mkdir -p", fdir))
    ofn = paste(fdir, "/contig_marginals.png", sep="")
    cat(sprintf("generating plot: %s\n", ofn))
    png(ofn, 600, 600)
    hist(t)
    abline(v=threshold, lty=2)
    dev.off()
}

plot.complete.graph=function(contigs.ifn, contacts.ifn, anchor.ifn, known.ifn, fdir)
{
    plot.graph=function(tcontacts, tcontigs, legend.colors, legend.text, ofn, contig.color.field) {
        tcontacts$col = ifelse(tcontacts$n < 10, "lightgrey", "darkgrey")

        tcontigs$col = tcontigs[,contig.color.field]

        cat(sprintf("generating plot: %s\n", ofn))
        png(ofn, 600, 600)
        plot.new()
        plot(tcontigs$x, tcontigs$y, type="p", col=NA)
        segments(tcontacts$x1, tcontacts$y1, tcontacts$x2, tcontacts$y2, col=tcontacts$col)
        points(tcontigs$x, tcontigs$y, col=tcontigs$col, pch=19, cex=0.5)
        legend("topright", fill=legend.colors, legend=legend.text)
        dev.off()
    }

    library(igraph)
    options(stringsAsFactors=F)
    contigs = read.delim(contigs.ifn)
    contacts = read.delim(contacts.ifn)
    contacts$n = contacts$contacts
    contacts = contacts[contacts$n > 0,]

    gg = graph_from_data_frame(as.matrix(contacts[,1:2]), directed = FALSE, vertices = contigs)
    cc = as.data.frame(layout_nicely(gg))

    contigs$x = cc[,1]
    contigs$y = cc[,2]

    contacts$x1 = contigs$x[match(contacts$contig1, contigs$contig)]
    contacts$y1 = contigs$y[match(contacts$contig1, contigs$contig)]

    contacts$x2 = contigs$x[match(contacts$contig2, contigs$contig)]
    contacts$y2 = contigs$y[match(contacts$contig2, contigs$contig)]

    anchors = read.delim(anchor.ifn)
    contigs$anchor = anchors$cluster[match(contigs$contig, anchors$contig)]
    contigs$anchor[is.na(contigs$anchor)] = -1
    contigs$anchor_col = ifelse(contigs$anchor == -1, "black", contigs$anchor+1)

    N = max(contigs$anchor)
    colors = 1:(N+1)
    legend = c("NA", 1:N)

    system(paste("mkdir -p", fdir))
    ofn = paste(fdir, "/contig_graph_anchors.png", sep="")
    plot.graph(tcontacts=contacts, tcontigs=contigs, legend.colors=colors, legend.text=legend, ofn=ofn, contig.color.field="anchor_col")

    if (!file.exists(known.ifn))
        return (0)

    known = read.delim(known.ifn)
    contigs$genome = known$ref[match(contigs$contig, known$contig)]
    contigs$genome[is.na(contigs$genome)] = "none"
    genomes = unique(contigs$genome)
    contigs$genome.index = match(contigs$genome, genomes)
    N = max(contigs$genome.index)
    contigs$genome_col = ifelse(contigs$genome == "none", "black", ifelse(contigs$genome == "multi", "orange", contigs$genome.index+1))
    colors = c("black", "orange", 1 + 1:N)
    legend = c("none", "multi", genomes[genomes != "none"& genomes != "multi"])

    ofn = paste(fdir, "/contig_graph_known.png", sep="")
    plot.graph(tcontacts=contacts, tcontigs=contigs, legend.colors=colors, legend.text=legend, ofn=ofn, contig.color.field="genome_col")
}

plot.cluster.tree=function(tree.ifn, contigs.ifn, cluster.table.ifn, fdir)
{
    options(stringsAsFactors=F)
    library(gplots)

    tree = read.delim(tree.ifn)
    contigs = read.delim(contigs.ifn)
    contig.cluster = read.delim(cluster.table.ifn)
    contigs$cluster = contig.cluster$cluster[match(contigs$contig, contig.cluster$contig)]
    contigs$cluster[is.na(contigs$cluster)] = -1
    contigs$col = ifelse(contigs$cluster == -1, "black", contigs$cluster+1)

    xlim = c(0, dim(contigs)[1])
    ylim = c(-1, max(tree$level))

    ncols = 1024
    panel = colorpanel(ncols, "white", "red")
    Min = min(tree$score[tree$score!=-1])
    Max = max(tree$score)
    tree$col.i = ifelse(tree$score == -1, 1, 1+floor(ncols*(tree$score - Min)/Max))
    tree$col = panel[tree$col.i]

    system(paste("mkdir -p", fdir))
    ofn = paste(fdir, "/contig_tree.png", sep="")
    cat(sprintf("generating plot: %s\n", ofn))
    png(ofn, 600, 600)
    plot.new()
    plot.window(xlim=xlim, ylim=ylim)
    rect(tree$start_contig, tree$level, tree$end_contig, tree$level+1, col=tree$col, border=NA)
    rect(contigs$order, -1, contigs$order+1, 0, col=contigs$col, border=NA)
    dev.off()
}
