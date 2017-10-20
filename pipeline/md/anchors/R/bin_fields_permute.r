bin.fields=function(ifn, ofn, field, permute)
{
    table = load.table(ifn)
    values = unique(table[,field])

    # permute values
    if (permute) values = sample(values)
    N = length(values)
    bins = 1:N

    # save bin table
    bin.table = data.frame(bin=bins, value=values)
    names(bin.table)[2] = field
    save.table(bin.table, paste(ofn, ".", field, sep=""))

    # append bins
    table[,paste(field, "_bin", sep="")] = bin.table$bin[match(table[,field], bin.table[,field])]
    save.table(table, ofn)
}
