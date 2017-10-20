library(Biostrings)

make.reads=function(ifn, odir)
{
    contigs = readLines(ifn)
    M = length(contigs)/2
    ids = gsub(">", "", contigs[(1:M)*2-1])
    seqs = contigs[(1:M)*2]

    # one copy of 1st contig
    v1 = seqs[1]

    # two copies of 2nd contig
    v2a = seqs[2]
    v2b = seqs[2]

    v3 = seqs[3]

    cc = c("A", "C", "G", "T")
    del=function(v, from, length) {
        paste(c(substr(v, 1, from-1), substr(v, from+length, nchar(v))), collapse="")
    }
    insert=function(v, coord, seq) {
        paste(c(substr(v, 1, coord), seq, substr(v, coord+1, nchar(v))), collapse="")
    }
    subs=function(v, coord, i=1) {
        substr(v, coord, coord) = cc[(match(substr(v, coord, coord), cc) + i - 1) %% 4 + 1]
        v
    }

    do.sub = F
    do.del = F
    do.insert = T

    # subs
    if (do.sub) {
        v1 = subs(v1, 404)
        v1 = subs(v1, 12)
        v2a = subs(v2a, 303)
        v2b = subs(v2b, 303, 2)
        v3 = subs(v3, 13)
    }

    # deletion
    if (do.del) {
        v1 = del(v1, 105, 4)
        v2a = del(v2a, 15, 6)
        v2b = del(v2b, 12, 6)
        v3 = del(v3, 20, 1)
    }

    # insertions
    if (do.insert) {
        v1 = insert(v1, 130, "GAAGCTC") # GAGCTTC
        v2a = insert(v2a, 60, "CAAAAC")
        v3 = insert(v3, 90, "GCCCCCA")
    }

    # final seq
    seq = paste(c(seqs[-3], v1, v2a, v2b, v3), collapse="")
    seq = paste(seq, seq, seq, sep="", collapse="")

    read.count = 20000
    read.len = 120

    rc = function(seq) { as.character(reverseComplement(DNAStringSet(seq))) }
    seq = paste(c(seq, rc(seq)), sep="", collapse="")
    N = nchar(seq) - read.len - 1

    # generate reads
    start = round(runif(read.count, 1, N))
    end = start + read.len - 1

    result = NULL
    for (i in 1:read.count) {
        id = paste("@r_", i, sep="", collapse="")
        read = substr(seq, start[i], end[i])
        quality = paste(rep("I", nchar(read)), collapse="")
        result = c(result, id ,read, "+", quality)
    }
    lcount = read.count*4
    r1 = result[1:(lcount/2)]
    r2 = result[(lcount/2+1):lcount]

    save=function(r, ofn) {
        fc = file(ofn)
        writeLines(r, fc)
        close(fc)
    }
    ofn1 = paste(odir, "/R1.fastq", sep="")
    ofn2 = paste(odir, "/R2.fastq", sep="")
    save(r=r1, ofn=ofn1)
    save(r=r2, ofn=ofn2)
}

