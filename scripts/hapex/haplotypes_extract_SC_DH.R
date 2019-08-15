setwd('~/ma/data/dh')
library(data.table)
library(tidyr)
library(dplyr)
library(reshape2)
options(max.print = 999)

## extract haplotypes in a given window and save its frequencies and sequences
filename <- '260419_step6_red_v4_SC_DH.vcf' # adjust for each population and type
#(outname <- paste0("haplocount/hapseq_55k_50kb_v4_",substr(filename, 17,21),".txt")) # 600k data
outname <- paste0("haplocount/hapseq_55k_v4_",substr(filename, 21,25),".txt") # 55k data

dat <- fread(filename, data.table=T)
map <- fread('../map/50kb_v4_bed_breakpoints.bed', data.table = T, col.names = c('chr','start','stop'))
chroms <- 1:10

genosplit <- function(v, allele){      # allele is 1 or 3
    indexer <- seq(allele,nrow(v)*3,3)
    gen <- lapply(v, function(x) unlist(sapply(x, function(y) strsplit(y, '|'))))
    li <- vector()
    for (x in 1:length(gen)){li <- append(li,paste0(gen[x][[1]][indexer], collapse = ""))}
    return(li)
}


winlst <- list()
biglst <- list()
for (chrom in chroms){
    cat(paste("Running chromosome", chrom, '...\n'))
    windows1 <- map[chr==chrom]$start
    windows2 <- map[chr==chrom]$stop
    datchr <- dat[`#CHROM`==chrom]
    pb <- txtProgressBar(min = 0, max = length(windows1), style = 3)
    for (i in 2:length(windows1)){      # skip the first window
        region <- datchr[POS>=windows1[i] & POS<windows2[i], 10:ncol(datchr)]
        winlst["chr"]  <- chrom
        winlst["pos"]  <- windows1[i]
        if (nrow(region)==0){
            winlst["snps"] <- 0
            next 
        }
        haptab <- table(append(genosplit(region,1),genosplit(region,3)))
        winlst["snps"] <- nchar(names(haptab)[1])
        winlst["haps"] <- dim(haptab)
        winlst["freq"] <- list(haptab)
        biglst[[paste0(chrom,':',windows1[i])]] <- winlst
        setTxtProgressBar(pb, i)
    }
    close(pb)
}
cat("Done.\n")

bigdf <- data.frame(t(sapply(biglst, '[', seq(max(sapply(biglst,length))))), # convert to data.frame
                    stringsAsFactors = F)

hapdf <- melt(bigdf$freq, as.is=T)
names(hapdf) <- c('seq','count', 'chrompos')
bigdf$chrompos <- row.names(bigdf)
hapdf <- full_join(hapdf, select(bigdf, -freq), by="chrompos")

hapdf$freq <- round(hapdf$count/max(hapdf$count),3) # get frequencies from counts
hapdf <- hapdf %>% group_by(chrompos) %>%             # calculate haplotype diversity
    mutate(hapdiv=round(1-sum(freq**2),3)) %>% data.frame()

fwrite(hapdf, outname, col.names=T, quote=F, sep="\t", na = "NA")  

