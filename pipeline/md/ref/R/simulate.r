##################################################################################################
# generate read coords
##################################################################################################

generate.abundance.table=function(ifn, abundance.min, abundance.max, ofn)
{
    options(stringsAsFactors=F)
    table = read.delim(ifn)
    table$abundance = runif(dim(table)[1], min=abundance.min, max=abundance.max)

    cat(sprintf("saving table: %s\n", ofn))
    write.table(table, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

##################################################################################################
# generate read coords
##################################################################################################

generate.assembly.read.coords.table=function(genome.table.ifn, xcoverage.base, base.length, insert, insert.sd, read.length=250, ofn)
{
    options(stringsAsFactors=F)
    table = read.delim(genome.table.ifn)
    genomes = table$genome

    result = NULL

    table$xcoverage = xcoverage.base * table$abundance
    total.genome.size = sum(table$length)
    N.total = ceiling(sum(table$xcoverage*table$length) / (2*read.length))
    cat(sprintf("covering %.1fMb genome with base xcoverage of %.1f, expected number of reads: %.1fM\n",
                total.genome.size/10^6, xcoverage.base, N.total/10^6))

    cat(sprintf("going over %d genomes...\n", length(genomes)))
    for (genome in genomes) {
        genome.size = table$length[match(genome, table$genome)]
        xcoverage = table$xcoverage[match(genome, table$genome)]
        N = ceiling((xcoverage*genome.size) / (2*read.length))
        start = floor(runif(N, min=0, max=genome.size))
        end = pmax(start + base.length, start + insert + round(rnorm(N, sd=insert.sd)))
        strand = ifelse(runif(N) > 0.5, 1, -1)

        start1 = ifelse(strand ==  1, start, end-read.length)
        start2 = ifelse(strand == -1, start, end-read.length)

        df = data.frame(
            contig1=genome, start1=start1, end1=start1+read.length, strand1=strand,
            contig2=genome, start2=start2, end2=start2+read.length, strand2=-strand)

        result = rbind(result, df)
    }
    cat(sprintf("saving %d reads to file: %s\n", dim(result)[1], ofn))
    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}

generate.hic.read.coords.table=function(genome.table.ifn, nreads=10^6, exponent, p.noise=0.2, p.distal=0.5, read.length=250, ofn)
{
    cat(sprintf("reads allocated for noise: %.1f%%\n", 100*p.noise))

    p.cis = 1 - p.noise

    table = load.table(genome.table.ifn)
    table$weight = table$length * table$abundance.hic / sum(table$length * table$abundance.hic)
    table$nreads = round(table$weight * nreads * p.cis)

    genomes = table$genome
    max.genome.size = max(table$length)

    result = NULL

    # cis
    for (genome in genomes) {
        genome.size = table$length[match(genome, table$genome)]
        N = table$nreads[match(genome, table$genome)]

        distal = runif(N) <= p.distal

        # decay
        d.min = read.length * 2
        d.max = genome.size
        d = seq(from=d.min, to=d.max)
        y = d^exponent
        y = y/sum(y)
        cs = cumsum(y)
        ind = findInterval(runif(N), cs)
        ind[ind == 0] = 1
        cis.dist = d[ind]
        cis.dist[cis.dist > genome.size] = genome.size

        # random strand
        strand1 = ifelse(runif(N) > 0.5, 1, -1)
        strand2 = ifelse(runif(N) > 0.5, 1, -1)

        # we approximate distal by any random coord on the same genome
        dcoord = floor(runif(N, min=1, max=genome.size))

        # side R1
        genome1 = rep(genome, N)
        coord1 = floor(runif(N, min=0, max=genome.size))

        # side R2
        genome2 = genome1
        coord2 = ifelse(distal, dcoord, coord1 + cis.dist)

        gtable1 = data.frame(genome=genome1, coord=coord1)
        gtable2 = data.frame(genome=genome2, coord=coord2)

        df = data.frame(
            contig1=genome1, start1=coord1, end1=coord1+read.length, strand1=strand1,
            contig2=genome2, start2=coord2-read.length, end2=coord2, strand2=strand2)
        result = rbind(result, df)
    }

    # trans

    N = p.noise * nreads

    # strand
    strand1 = ifelse(runif(N) > 0.5, 1, -1)
    strand2 = ifelse(runif(N) > 0.5, 1, -1)

    # genome
    cs = c(0,cumsum(table$weight))
    cs[length(cs)] = 1
    genome1 = genomes[findInterval(runif(N), cs, all.inside=T)]
    genome2 = genomes[findInterval(runif(N), cs, all.inside=T)]

    # coord
    coord1 = floor(runif(N, min=1, max=100*max.genome.size)) %% table$length[match(genome1, table$genome)]
    coord2 = floor(runif(N, min=1, max=100*max.genome.size)) %% table$length[match(genome2, table$genome)]

    df = data.frame(
        contig1=genome1, start1=coord1, end1=coord1+read.length, strand1=strand1,
        contig2=genome2, start2=coord2, end2=coord2+read.length, strand2=strand2)

    result = rbind(result, df)

    cat(sprintf("generated trans reads: %d (%.1f%%)\n", sum(result$contig1 != result$contig2),
                100 * sum(result$contig1 != result$contig2) / dim(result)[1]))

    write.table(result, ofn, quote=F, col.names=T, row.names=F, sep="\t")
}
