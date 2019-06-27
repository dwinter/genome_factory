library(stringr)

genomic_interval <- function(chrom, start, end){
    structure(list(chrom, as.integer(start),as.integer(end)), .Names=c("chrom", "start", "end"))
}

bed_to_interval <- function(row){
    genomic_interval(row[["chrom"]], row[["start"]], row[["end"]])
}

gdist <- function(a,b){
    if(is.null(a) | is.null(b)){
        return(Inf)
    }
    if(a$chrom != b$chrom){
        return(Inf)
    }
    gap <- if (a$end > b$end) a$start - b$end else b$start - a$end
    if(gap <= 0){
        #overlap
        return(0)
    }
    gap

}

read_bed <- function(fname, skip_mt = TRUE){
    res <- read.table(fname, stringsAsFactors=FALSE)
    if(skip_mt){
        res <- res[res[,1] != 'mtDNA',]
    }
    names(res)[1:3] <- c("chrom", "start", "end")
    res
}

read_gtf <- function(fname, skip_mt=TRUE, parse_attr=TRUE){
    res <- read.table(fname, sep="\t", stringsAsFactors=FALSE)
    if(skip_mt){        
        res <- res[res[,1] != 'mtDNA',]
    }
    names(res) <- c("chr", "source", "feature", "start", "end", "score", "strand"," frame", "attribute")
    if(parse_attr){
        attr_cols <- gtf_attributes(res)
        res$attribute <- NULL
        res <- cbind.data.frame(res, attr_cols)
    }
    res
}

read_RM <- function(fname, skip_mt = TRUE){
    res <- read.table(fname, sep="\t", stringsAsFactors=FALSE, skip=3)

}


gtf_attributes <- function(gtf){
    split_one <- function(x){
        pairs <- str_split(x, " ")
        structure(as.list(sapply(pairs, "[[", 2)), .Names=sapply(pairs, "[[", 1))
    }
    tokens <- str_split(gtf$attribute, ";\\W+")
    do.call(rbind.data.frame, lapply(tokens, split_one))
}




bed_len <- function(B) B$end - B$start


#' take a numeric vector an converte to ordinal catergories
#' using provided quantiles
#' x <- rpois(1000,10)
#' table(ordinalize(x, c(0.9, .95, .99)))

ordinalize <- function(x, quantiles, lowest=-Inf, highest=Inf,...){
    breaks <- c( lowest, quantile(x, quantiles, ...), highest)
    cut(x, breaks, labels=FALSE)
}


geom_bed <- function(data, y, ...){
    data$y <- y
    geom_segment(data=data, aes(x=start, xend=end, y=y, yend=y), ...)
}

gg_genome <- function(bed, geom, score_col="score"){    
   geom(data=bed, aes_string("(start+end)/2", score_col))
}


ggtrack <- function(bed, strands=FALSE, y, ...){
    bed$y <- y
    if(strands){
        cnames <- names(bed)
        bed_cols <- c("chrom", "start", "end")
        other_cols <- cnames[!cnames %in% bed_cols]
        neg_strand <-  bed$strand == "-"
        B <- bed[neg_strand,c("chrom", "end", "start", other_cols)]
        names(B)[1:3] <- bed_cols
        A <- bed[neg_strand,c("chrom", "start", "end", other_cols)]
        bed <- rbind.data.frame(A,B)
        bed <- bed[ order(bed$chrom, bed$start, bed$end),]
        l <- geom_segment(data=bed, aes(x=start, xend=end, y=y, yend=y), ...)
        return(l)
    }
    geom_segment(data=bed, aes(x=start, xend=end, y=y, yend=y), ...)

}

geom_highlight <- function(bed, ymin=0, ymax=1, ...){
    bed$ymin <- ymin
    bed$ymax <- ymax
    geom_rect(data=bed, aes(xmin=start, xmax=end, ymax=ymax, ymin=ymin), ...)

}


matching_bed_region <- function(bedA, bedB){
    chroms <- unique(bedA$chrom)
    if(length(chroms) > 1){
        stop("Only works within on chromosome")
    }
    new_start <- min(bedA$start)
    new_end <-  max(bedA$end)
    idx <- (bedB$chrom == chroms) & (bedB$start >= new_start) & (bedB$end <= new_end)
    bedB[idx,]
}


theme_genome_track <- function(base_size = 11, base_family = ""){
    theme(
     axis.line = element_line(color="black", size = 0.5),
#     axis.line.y = element_line(color="black", size = 2),
     panel.grid.major = element_blank(),
     panel.grid.minor = element_blank(),
     panel.border = element_blank(),
     panel.background = element_blank(),
     #line = element_blank(), 
     rect = element_blank(), 
     text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, lineheight = 0.9,  hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(), debug = FALSE),
     legend.text = element_text(size = rel(0.8)), 
     legend.title = element_text(hjust = 0),
     strip.text = element_text(size = rel(0.8)), 
     plot.margin = unit(c(0, 0, 0, 0), "lines"), 
     complete = TRUE
   )
}


bed_to_gff <- function(bed, feat_source="Massey", use_six_col=TRUE){
    res <- bed
    if(ncol(bed) >= 6 & use_six_col){
        res$strand <- bed[,6]
        res$score <- bed[,5]
        res$feature <- bed[,4]
    }
    else{
        res$strand <- "+"
        res$score <- "."
        res$feature <- "unknown_feature"
    }
    res$source <- feat_source
    res$group <- "."
    res$frame <- "."
    res[,c("chrom", "source", "feature", "start", "end", "score", "strand", "frame", "group")]
}

gene_centres <- function(bed, gene_name="gene", chrom=TRUE){
    starts <- aggregate( as.formula(paste0("start ~ ",gene_name)), FUN=min, data=bed)
    ends <- aggregate( as.formula(paste0("end ~ ",gene_name)), FUN=max, data=bed)
    res <- data.frame(gene=starts[,1], mid=(starts[,2] + ends[,2])/2)
    if(chrom){
       res$chroms <- aggregate( as.formula(paste0("chrom ~ ",gene_name)), FUN=unique, data=bed)[,2]
    }
    res   
}


