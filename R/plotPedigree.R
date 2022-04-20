#' Plot Pedigree
#'
#' Plotting Pedigree tree
#'
#' @param nodelist input file containing information on generation times and pedigree lineages
#' "extdata" called "nodelist.fn"
#' @param edgelist input file containing edges "edgelist.fn"
#' @param sampling.design "progenitor.intermediate"; "sibling"; "progenitor.endpoint";"tree"
#' @param out.pdf output file name
#' @param output.dir output directory
#' @param plot.width  plotting with
#' @param plot.height plotting height
#' @param vertex.label label vertix
#' @param vertex.size size of vertix
#' @param aspect.ratio aspect.ration
#' @import dplyr
#' @importFrom  data.table fread
#' @importFrom  data.table fwrite
#' @importFrom  data.table as.data.table
#' @importFrom  utils combn
#' @importFrom  igraph graph_from_data_frame
#' @importFrom  igraph layout_as_tree
#' @importFrom  igraph all_shortest_paths
#' @importFrom  igraph V
#' @importFrom  igraph V<-
#' @importFrom  igraph degree
#' @importFrom  igraph plot.igraph
#' @importFrom  graphics plot
#' @return plot pedigree matrices file.
#' @export
#' @examples
#'# Get some toy data
#' file <- system.file("extdata", "dm","nodelist.fn", package="AlphaBeta")
#' file2 <- system.file("extdata","dm","edgelist.fn", package="AlphaBeta")
#' plotPedigree(nodelist = file, edgelist=file2, sampling.design="sibling",vertex.label=TRUE,
#'  out.pdf="Plot", output.dir=getwd() )

plotPedigree <- function(nodelist, edgelist, sampling.design, out.pdf=NULL, output.dir=NULL,
                         plot.width=11, plot.height=8, vertex.label=NULL, vertex.size=12, aspect.ratio=2.5)
  {
    samples <- fread(nodelist, header=TRUE, skip=0, select = c(2,3,4))
    edges   <- fread(edgelist, header=TRUE, skip=0)

    #graph layout from edges and nodes
    g <- graph_from_data_frame(edges, directed=TRUE, vertices=samples)

    #put root vertex on top with flip.y=TRUE
    lay <- layout_as_tree(g, flip.y=TRUE)

    #-----------------------------------------------------
    if ((sampling.design == "progenitor.intermediate")) {
      df <- lay
    }

    #-----------------------------------------------------
    if ((sampling.design == "progenitor.endpoint"))
    {
      #the generation difference as weights
      weight.scale <- c(1, edges$gendiff)

      #groups from user input
      grp <- c("A", edges$group)
      m <- cbind(lay, weight.scale, grp)
      m <- data.frame(m)

      #split groups as individual dataframes
      # gps <- m %>% dplyr::group_by(grp)
      # mygps <- dplyr::group_split(gps)
      mygps <- split(gps, gps$grp)

      #set gap-width for complex pedigrees. gap of 2 is good enough.
      gap = 2

      for (j in seq_along(mygps)) {
        if (j==1){
          #relevel the pedigree
          mygps[[j]]$relevel <- as.numeric(mygps[[j]]$V2)
        }
        else {
          #find minimum (last) level of previous group. Use it for re-scaling the pedigree
          mm <- min(as.numeric(mygps[[j-1]]$relevel))
          mygps[[j]]$relevel <- as.numeric(mygps[[j]]$V2)-gap-mm
        }
      }

      m <- bind_rows(mygps)
      m <- as.data.frame(m, stringsAsFactors = FALSE)
      lay[,2] <- m$relevel
      df <- lay
    }

    #-----------------------------------------------------
    if ((sampling.design == "sibling")) {

      #find root of the pedigree
      root <- which(degree(g, v=V(g), mode = "in")==0, useNames=TRUE)
      #find degree of the node
      d <- degree(g, root)

      mg <- merge(samples[2:length(samples$node),], edges, by.x="node",by.y="to", sort=FALSE)

      myg <- c("A", mg$group)
      m <- cbind(lay, samples, myg)
      m <- data.frame(m)

      #split groups as individual dataframes
      gps <- m %>% group_by(myg, meth)
      mygps <-group_split(gps)

      n = length(mygps)
      for (i in 1:n) {
        mygps[[i]]$V1 = i-0.5*(n+1)
      }

      mm <- bind_rows(mygps)
      mm  <- as.data.frame(mm, stringsAsFactors = FALSE)

      if (as.numeric(d) > 1){
        mm[1,1]=0
      }

      df <- merge(m, mm, by="node", sort=FALSE)
      df <- df[,c("V1.y", "V2.y")]
      colnames(df) <- c("V1","V2")
      df <- as.matrix(df)
    }

    #-----------------------------------------------------
    if ((sampling.design == "tree")) {

      #for trees with 1 stem. "Stem" column not needed.
      #---------------------------------------------------
      if (is.null(edges$Stem)) {
        m <- cbind(lay, samples)
        m <- data.frame(m)

        m <- m %>% dplyr::mutate(V1 = ifelse(meth=="N", 0, V1))

        i=1
        j <- length(which(m$meth=="Y"))
        k=0

        while (i <= (length(m$V1))){
          if ("Y" %in% m$meth[i])
          {
            #print(m[i,])
            m$V1[i]=j
            j=j-1
          } else {
            m$V2[i]=k
            k=k-1
          }
          i=i+1
        }
        #5 no. of meth == N
        nl <- length(which(m$meth=="N"))
        m <- m %>% mutate(V2 = ifelse(meth=="Y", -nl, V2))

        mm <- m[,c(1,2)]
        mm[,1] <- as.numeric(as.character(m$V1))
        mm[,2] <- as.numeric(as.character(m$V2))
        df  <- as.matrix(mm)
      }
      else {
        #for trees with 2 stem (rare scenario)
        #---------------------------------------------------
        tot.age <- as.numeric(samples[which(samples$Branchpoint_date==0),]$node)
        mygroup <- c(tot.age, edges$Stem)
        m <- cbind(lay, samples, mygroup)
        m <- data.frame(m)

        #for samples with methylation measurements, level branch to 0
        m <- m %>% mutate(V2 = ifelse(meth=="Y", 0, V2))

        mycount <- m %>% count(mygroup)

        i=0
        while (i < (length(m$V1)-1) ){
          i=i+2
          m$V1[i]=if (i <= mycount$n[1]) 0.5*(-i) else 0.5*(i-mycount$n[1]) # 0.5 for sample pairs
          m$V1[i+1]=if (i <= mycount$n[1]) 0.5*(-i) else 0.5*(i-mycount$n[1])
        }
        mm <- m[,c(1,2)]
        mm[,1] <- as.numeric(as.character(m$V1))
        mm[,2] <- as.numeric(as.character(m$V2))
        df  <- as.matrix(mm)
      }
    }

    #-----------------------------------------------------------

    V(g)$color <- ifelse(samples[V(g), 3] == "Y", "red", "gray")

    if (!is.null(out.pdf)) {
      pdf(file.path(output.dir, paste0(out.pdf, ".pdf")), colormodel = 'cmyk', width = plot.width, height = plot.height)
    }

    pl <- plot.igraph(g, layout = df,
                      asp=aspect.ratio,
                      #edge.width = 1,
                      vertex.size = vertex.size,
                      vertex.frame.width = 1,
                      vertex.label = if (vertex.label==FALSE) NA,
                      vertex.label.cex = 0.6,
                      vertex.label.color = "black",
                      vertex.label.dist = 1,
                      vertex.label.degree = -pi/2,
                      edge.arrow.size = 0.1,
                      vertex.color = V(g)$color
    )
    print(pl)

    if (!is.null(out.pdf)) {
      dev.off()
    }
  }


