export.table=function(ofn, ...)
{
    vars = list(...)

    result = NULL
    for (i in 1:length(vars)) {
        name = names(vars)[i]
        var = vars[[i]]
        result = rbind(result, data.frame(id=name, value=var))
    }
    save.table(result, ofn)
}

export=function(odir, ...)
{
    vars = list(...)
    cat(sprintf("copying files to directory: %s\n", odir))
    for (i in 1:length(vars)) {
        name = names(vars)[i]
        ifn = vars[[i]]
        ofn = paste(odir, "/", name, sep="")
        command = sprintf("cp %s %s", ifn, ofn)
        cat(sprintf("copying file: %s\n", ifn))
        if (system(command) != 0) {
            stop(paste("failed command:", command))
        }
    }
}

# rename specific fields
copy.table=function(ifn, ofn, fields.in, fields.out)
{
    if (length(fields.in) != length(fields.out))
        stop(paste("in and out fields must be the same length"))

    df = load.table(ifn)
    for (field in fields.in) {
        if (!is.element(field, names(df))) {
            stop(paste("field does not exist:", field))
        }
    }

    df.out = df[,match(fields.in, names(df))]
    colnames(df.out) = fields.out

    save.table(df.out, ofn)
}
