
library(IRanges)

# Retreive for Uniprot all GO annotations of the given organism
query_uniprot <- function(org) {
  url <- "https://www.uniprot.org/uniprot/?format=tab&columns=id,entry%20name,protein%20names,genes(PREFERRED),genes,ec,go-id,keyword-id"
  url <- paste0(url,"&query=organism:",org)
  uniprot <- read.table(url,sep="\t",comment="",quote="",stringsAsFactors=FALSE,header=TRUE,check.names=FALSE)
  uniprot <- as(uniprot,"DataFrame")
  
  uniprot$Gene_names <- CharacterList(strsplit(uniprot$"Gene names","( +)|( */ *)"))
  uniprot$Gene_ontology_IDs <- CharacterList(strsplit(uniprot$"Gene ontology IDs"," *; *"))
  uniprot$EC_number <- CharacterList(strsplit(uniprot$"EC number"," *; *"))
  uniprot$Keyword_IDs <- CharacterList(strsplit(uniprot$"Keyword ID"," *; *"))
  within(uniprot,`Gene names` <- `Gene ontology IDs` <- `EC number` <- `Keyword ID` <- NULL)
}

uniprot_keywords_obo_con <- function() {
  gzcon(url("https://www.uniprot.org/keywords/?query=*&format=obo&force=true&compress=yes"))
}


