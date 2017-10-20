ref.table=function(id.ifn, ref.ifn, cluster.ifn, ofn)
{
    id.df = load.table(id.ifn)
    ref.df = load.table(ref.ifn)
    cluster.df = load.table(cluster.ifn)

    df = id.df
    df$genome = ref.df$genome[match(df$accession, ref.df$genome_id)]
    df$id = cluster.df$id[match(df$genome, cluster.df$set)]
    df = df[match(cluster.df$id, df$id),]

    names = c("id", "genome", "accession", "organism_name", "infraspecific_name")
    alt.names = c("id", "genome", "accession", "name")
    if (all(is.element(names, names(df)))) {
        df = df[,names]
    } else {
        df = df[,alt.names]
    }
    save.table(df, ofn)
}
