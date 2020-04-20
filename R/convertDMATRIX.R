
convertDMATRIX<-function(samples, edges, dmatrix)
{
    # reading matrix
    g <- graph_from_data_frame(edges, directed=TRUE, vertices=samples)

    # constracting tree
    samples.with.WGBS <- subset(samples, samples$meth == "Y")
    expand.grid.unique <- function(x, y, include.equals=FALSE)
    {
        x <- unique(x)
        y <- unique(y)
        k <- function(i)
        {
            z <- setdiff(y, x[seq_len(i-include.equals)])
            if(length(z)) cbind(x[i], z, deparse.level=0)
        }
        do.call(rbind, lapply(seq_along(x), k))
    }
    # making pedigres
    pairs.for.D <- expand.grid.unique(as.character(samples.with.WGBS$node), as.character(samples.with.WGBS$node))
    t1<-NULL
    t2<-NULL
    t0<-NULL
    s1<-NULL
    s2<-NULL
    for (a in 1:nrow(pairs.for.D))
    {
        # Loop the to and from over a distance matrix that matches D (only measured samples)
        path.this<-all_shortest_paths(g, from=pairs.for.D[a,1], to=pairs.for.D[a,2], mode="all")
        path.this<-names(unlist(path.this$res))
        s1[a]<-pairs.for.D[a,1]
        s2[a]<-pairs.for.D[a,2]
        t1[a]<-samples[which(samples$node == pairs.for.D[a,1]),2]$gen
        t2[a]<-samples[which(samples$node == pairs.for.D[a,2]),2]$gen
        samples.temp<-samples[which(as.character(samples$node) %in% path.this),2]
        t0[a]<-min(samples.temp)
    }
    # making data-frame
    tmpData <- data.frame(s1, s2, t0, t1, t2, stringsAsFactors=FALSE)
    colnames(tmpData)[3]<-"time0"
    colnames(tmpData)[4]<-"time1"
    colnames(tmpData)[5]<-"time2"
    # join Dmatrix with pedigree time
    tm1<-as.data.table(inner_join(dmatrix,tmpData ,  by = c("pair.1"="s1","pair.2"="s2")))
    tm2<-as.data.table(anti_join(dmatrix,tmpData ,  by = c("pair.1"="s1","pair.2"="s2")))
    tm2<-as.data.table(inner_join(tm2,tmpData ,  by = c("pair.1"="s2","pair.2"="s1")))
    dmatrixout<-as.data.table(rbind(tm1,tm2))
    # reordering table
    dmatrixout <- dmatrixout[, c(1, 2, 4, 5, 6, 3)]

    dmatrixout <- as.matrix(dmatrixout[,3:6])

    return(dmatrixout)
}
