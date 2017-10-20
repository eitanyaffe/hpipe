args = commandArgs(T)
ifn = args[1]
prefix = args[2]

x = read.delim(ifn)
result = sprintf("-model_num %d", dim(x)[1])
for (i in 1:dim(x)[1]) {
    size = x$size[i]
    field = x$field[i]
    result = paste(result, sprintf("-f_model_field_%d %s -f_model_size_%d %d -f_model_fn_%d %s_%s.f",
        i, field,
        i, size,
        i, prefix, field))
}
cat(result)
