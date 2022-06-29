

library(IRanges)
library(ggplot2)

#' Test GOterm enrichement for the given list of genes
#' @param query a logical matrix of selected genes (rownames are gene symbols, colnames are names of the geneset)
#' @param genes_annotations a CharacterList that list for each gene the sets it belong to. 
#'        The list should be limited to the gene of the universe, for example genes expressed in the brain.
multi_gsea <- function(query,genes_annotations) {
    stopifnot(is.matrix(query))
    stopifnot(is.logical(query))
    
    # remove duplicates and only keep genes that are annotated
    genes_annotations <- unique(genes_annotations)
    genes_annotations <- genes_annotations[lengths(genes_annotations)>0,]
    genes_annotations <- genes_annotations[names(genes_annotations) %in% rownames(query)]
    
    query <- query[rownames(query) %in% names(genes_annotations),,drop=FALSE]
    
    gsea <- with(stack(genes_annotations),splitAsList(as.character(name),value))
    gsea <- DataFrame(
      geneset_name = names(gsea),
      geneset_genes = gsea
    )
    gsea$geneset_size <- lengths(gsea$geneset_genes)
    
    gsea$overlap_genes <- as(apply(query,2,function(v) {
      gsea$geneset_genes[gsea$geneset_genes %in% rownames(query)[v]]
    }),"DataFrame")
    gsea$overlap_size <- sapply(gsea$overlap_genes,lengths)
    gsea$query_size <- matrix(rep(colSums(query,na.rm=TRUE),each=nrow(gsea)),nrow(gsea))
    gsea$universe_size <- matrix(rep(colSums(!is.na(query)),each=nrow(gsea)),nrow(gsea))
    
    gsea$expected_overlap_size <- gsea$query_size * gsea$geneset_size / gsea$universe_size
    gsea$fold_enrichment <- gsea$overlap_size / gsea$expected_overlap_size
    gsea$enrichment_pval <- 1 - phyper(gsea$overlap_size-1L,gsea$geneset_size,gsea$universe_size-gsea$geneset_size,gsea$query_size)
    gsea$enrichment_qval <- array(p.adjust(gsea$enrichment_pval,"fdr"),dim(gsea$enrichment_pval),dimnames(gsea$enrichment_pval)) 
    
    # sort the output
    o <- gsea$enrichment_pval
    o[is.na(o)] <- +Inf
    o <- order(Biobase::rowMin(o))
    gsea[o,]
}



read_gmt <- function(con) {
  txt <- readLines(con)
  pat <- "^([^\t]*)\t([^\t]*)\t(.*)$"
  z <- DataFrame(
    STANDARD_NAME = sub(pat,"\\1",txt),
    EXTERNAL_DETAILS_URL = sub(pat,"\\2",txt),
    MEMBERS = CharacterList(strsplit(sub(pat,"\\3",txt),"\t",fixed=TRUE))
  )
}


read_gaf <- function(gaf_file) {
  txt <- readLines(gaf_file)
  stopifnot(grepl("^!gaf-version:",txt[1]))
  txt <- txt[!grepl("^!",txt)]
  gaf <- read.table(text = txt,sep="\t",quote="",comment.char = "")
  colnames(gaf) <- c("db","db_object_id","db_object_symbol","qualifier","go_id","db_ref","evidence","with_or_from","aspect","db_object_name","db_object_synonym","db_object_type","taxon","date","assigned_by","annotation_extension","gene_product_form_id")
  return(gaf)
}

gg_matrix_hcol <- function(m,fill,grp) {
  stopifnot(ncol(m)==length(fill))
  stopifnot(length(grp)==length(fill))
            
  p <- reshape2::melt(m)
  p$fill <- fill[p$Var2]
  p$grp <- grp[p$Var2]
  ggplot(p) + facet_wrap(.~grp,nrow=1) +
    geom_col(aes(x=value,y=Var1,fill=fill)) + 
    geom_vline(xintercept=0) +
    theme(panel.spacing = unit(0,"mm"),panel.grid = element_blank(),strip.text = element_text(angle=90,hjust=0)) + 
    xlab("") + ylab("")
}





gg_matrix_hcol2 <- function(m_up,m_dn,m_up_lab=m_up,m_dn_lab=m_dn) {
  stopifnot(identical(dim(m_up),dim(m_dn)))
  stopifnot(identical(dimnames(m_up),dimnames(m_dn)))
  Z <- data.frame(
    row = as.vector(row(m_up)),
    col = as.vector(col(m_up)),
    value_up = as.vector(m_up),
    value_dn = as.vector(m_dn),
    lab_up = as.vector(m_up_lab),
    lab_dn = as.vector(m_dn_lab)
  )
  p <- ggplot(Z) + facet_wrap(.~col,nrow=1,labeller = as_labeller(setNames(colnames(m_up),seq(ncol(m_up))))) +
    geom_rect(aes(xmin=0,xmax=value_up,ymin=row-0.4,ymax=row+0.4,fill="up")) +
    geom_rect(aes(xmin=0,xmax=-value_dn,ymin=row-0.4,ymax=row+0.4,fill="dn")) + 
    geom_text(aes(x=value_up/2,y=row,label=lab_up),hjust=0,size=3) +
    geom_text(aes(x=-value_dn/2,y=row,label=lab_dn),hjust=1,size=3) +
    geom_vline(xintercept=0) +
    theme(
      panel.spacing = unit(0,"mm"),
      panel.grid = element_blank(),
      strip.text = element_text(angle=90,hjust=0)
    ) + 
    xlab("") + ylab("")
}

