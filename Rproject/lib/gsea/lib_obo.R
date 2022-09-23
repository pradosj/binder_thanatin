
library(IRanges)
library(igraph)


# read OBO format and return an igraph with all "is_a" and "part_of" relations
read_obo <- function(obo_file) {
  pat <- "^([^:]+) *: *(.*)$"
  obo <- readLines(obo_file)
  obo <- obo[grepl("^\\[",obo) | grepl(pat,obo)]
  obo <- DataFrame(
    section = cumsum(grepl("^\\[",obo)),
    key = sub(pat,"\\1",obo),
    value = sub(pat,"\\2",obo)
  )
  obo <- subset(obo,section %in% section[key=="[Term]"])
  obo <- subset(obo,key %in% c("is_a","name","id","is_obsolete","namespace","relationship"))
  
  get_relationship <- function(obo) {
    R <- with(obo[obo$key=="relationship",],CharacterList(split(value,section)))
    R <- sub(" *!.*$","",R) # remove comments
    pat <- "^([^ ]+) ((GO:[0-9]+)|(KW-[0-9]+))$"
    stopifnot(all(grepl(pat,unlist(R))))
    DataFrame(
      relationship_type = sub(pat,"\\1",R),
      relationship_goid = sub(pat,"\\2",R)
    )
  }
  
  obo$section <- factor(obo$section)
  obo <- DataFrame(
    goid = with(subset(obo,key=="id"),drop(CharacterList(split(value,section)))),
    term = with(subset(obo,key=="name"),drop(CharacterList(split(value,section)))),
    is_obsolete = with(subset(obo,key=="is_obsolete"),drop(CharacterList(split(value,section)))) %in% "true",
    namespace = with(subset(obo,key=="namespace"),drop(CharacterList(split(value,section)))),
    is_a = with(subset(obo,key=="is_a"),CharacterList(split(value,section))),
    get_relationship(obo)
  )
  obo$is_a <- sub(" *!.*$","",obo$is_a)
  
  is_a <- stack(setNames(obo$is_a,obo$goid))
  is_a$type <- "is_a"
  relationship <- stack(setNames(obo$relationship_goid,obo$goid))
  relationship$type <- unlist(obo$relationship_type)
  
  go <- graph_from_data_frame(rbind(relationship,is_a),vertices = as.data.frame(obo[1:4]))
  V(go)$term_and_id <- sprintf("%s: %s",V(go)$name,V(go)$term)
  
  if (!is.dag(go)) warning("Loaded graph is not a DAG")
  return(go)
}



#' Find all reachable nodes from given nodes
#' @param h an igraph
#' @param nodes query nodes
#' @param mode graph walking mode
#' @return a 2 column data.frame listing all node_id reachable from given query_nodes
reachable <- function(h,nodes=V(h),mode="out") {
  lnk <- ego(h,vcount(h),nodes = nodes,mode = mode)
  lnk <- IntegerList(lapply(lnk,as.integer))
  names(lnk) <- nodes$name
  lnk <- stack(lnk,"query_id","node_id")
  lnk$node_id <- factor(V(h)$name,V(h)$name)[lnk$node_id]
  as.data.frame(lnk)
}



#' Unfold a graph into a tree
#' @param g an igraph
#' @return 
unfold_dag <- function(g,root="root") {
  stopifnot(is_dag(g))
  g <- g %>% 
    add_vertices(1,name="__root__") %>% # Create a single root node
    add_edges(rbind(V(g)[!.from(E(g))],vcount(g)+1L))
  
  leafs <- V(g)[!.to(E(g))]
  paths <- all_simple_paths(g,vcount(g),to=leafs,mode = "in")

  n <- paths %>%
    lapply(Reduce,f=function(a,b) paste0(a,"/",b),accumulate=TRUE) %>%
    unlist %>%
    unique
  n <- sub("^[^/]*","",n)
  n <- setdiff(n,c("","/",".","./"))
  e <- cbind(dirname(n),n)
  e <- e[rowSums(!array(e %in% n,dim(e)))==0L,]
  g <- graph_from_edgelist(e)
  revmap <- as.integer(basename(V(g)$name))
  
  
  list(tree = g,vertex_index = revmap)
}



go_basic_obo_con <- function() {
  url("http://purl.obolibrary.org/obo/go/go-basic.obo")
}













