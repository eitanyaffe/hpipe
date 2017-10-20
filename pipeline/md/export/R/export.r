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
