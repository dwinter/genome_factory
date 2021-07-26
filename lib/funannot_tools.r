funanno_from_all <- function(annot_all, all_genes){
    annot <- read.table(annot_all, 
                        sep="\t", 
                        col.names=c("gene", "annotation", "value"),
                        stringsAsFactors=FALSE, 
                        quote=NULL)
    by_type <- split(annot, annot$annotation)
    all_notes <- unfold_pair(by_type[["note"]])
    sm_idx <- grepl("SMCOG\\d+", all_notes$annotation)
    sm_df <- all_notes[sm_idx,]
    sm_final <- data.frame(gene=sm_df$gene, 
                           annotation="SMCOG", 
                           value=paste0(sm_df$annotation, "|", sm_df$value))
    notes_final <- all_notes[!sm_idx,]
    dbs <- unfold_pair(by_type$db_xref)
    GO <- by_type[["go_function"]]

    GO_final <- data.frame(gene=GO$gene, 
                           annotation="go_function", 
                           value = as.character(str_match(GO$value, "\\|(\\d+)\\|")[,2]))

    very_long <- rbind.data.frame(by_type[["name"]],
                                  by_type[["product"]],
                                  notes_final, 
                                  dbs, 
                                  GO_final,
                                  sm_final)
    agg <- aggregate_by_gene(very_long)
    M <- merge(data.frame(gene=all_genes), agg, all.x=TRUE)
    res <- spread(M, annotation, value, fill="")
    #res[,c("gene", "product", "name", "antiSMASH", "BUSCO", "CAZy", "effectorP", "go_function", "InterPro", "MEROPS", "PFAM", "SMCOG"),]
    res[,c("gene", "product", "name", "antiSMASH", "BUSCO", "CAZy", "EffectorP", "go_function", "InterPro", "MEROPS", "PFAM", "SECRETED", "SMCOG"),]

}
       # final_DB <- aggregate(value ~ gene + annotation, FUN=paste0, data=dbs, collapse=":")


unfold_pair <- function(annot_tab){
    wide <- do.call(rbind.data.frame, strsplit(annot_tab[["value"]], ":"))
    names(wide) <- c("annotation", "value")
    data.frame(gene=annot_tab[["gene"]], wide)
}

aggregate_by_gene <- function(tab){
    aggregate(value ~ gene + annotation, FUN=paste0, data=tab, collapse=":")
    
}
