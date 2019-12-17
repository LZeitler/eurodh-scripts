source('source_me.R')

freqallele <- function(x,a) {
    x <- x[!is.na(x)]
    length(x[x==a])/length(x[!is.na(x)])
}

racess <- races <- c('BU','GB','RT','SC','SF')

## get datasets to compare before/after
majmin <- list()
for (race in races){
    imp <- fread(paste0('~/ma/data/dh/pipe/full4_',race,'/030718_step5_red_full4_',race,'_step2_imp.vcf'),data.table=F) # after imp
    dro <- fread(paste0('~/ma/data/dh/pipe/full4_',race,'/030718_step5_red_full4_',race,'_onlydropped.dh.vcf'),data.table=F) # before imp
    
    imp <- filter(imp,ID%in%dro$ID)

    dro <- dro %>% arrange(`#CHROM`,POS)
    
## how often does hom-alt (minor allele in dh-pop) change to hom-ref

    dm <- data.frame(matrix(nrow=nrow(dro),ncol=length(names(dro))-9)) # for imputation correct/incorrect
    names(dm) <- names(dro)[10:length(names(dro))]
    gtb <- mm <- dm                         # for before genotypes (maj/min)
    for (n in names(dm)){
        be <- be.na <- substr(dro[,n],1,1)        # has missing data
        be[be.na=='.'] <- NA                      # missing data is not ref or alt
        gtb[,n] <- be
        af <- substr(imp[,n],1,1)
        af[be.na=='.'] <- NA
        tf <- af==be
        dm[,n] <- tf
    }

    freq0 <- sapply(1:nrow(gtb), function(x) freqallele(gtb[x,],0))
    freq1 <- sapply(1:nrow(gtb), function(x) freqallele(gtb[x,],1))
    
    if (unique(freq1+freq0)==1) cat('Allele frequency sanity check ok.\n') else cat('Allele frequecies not ok.\n')
    
    maj <- freq0>=.5                        # T is major allele is 0, F if 1
    
    mm <- sapply(1:nrow(dm), function(r) {
        sapply(1:ncol(dm), function(c) {
            if(maj[r]) {                    # major allele is 0
                if(!is.na(dm[r,c])) {
                    if(dm[r,c]==T & gtb[r,c]==0) o <- 'major correct'
                    if(dm[r,c]==T & gtb[r,c]==1) o <- 'minor correct'
                    if(dm[r,c]==F & gtb[r,c]==0) o <- 'major incorrect'
                    if(dm[r,c]==F & gtb[r,c]==1) o <- 'minor incorrect'
                    o
                }
            } else {                        # major allele is 1
                if(!is.na(dm[r,c])) {
                    if(dm[r,c]==T & gtb[r,c]==0) o <- 'minor correct'
                    if(dm[r,c]==T & gtb[r,c]==1) o <- 'major correct'
                    if(dm[r,c]==F & gtb[r,c]==0) o <- 'minor incorrect'
                    if(dm[r,c]==F & gtb[r,c]==1) o <- 'major incorrect'
                    o
                }
            }
        })
    })
    mm <- t(mm)

    majmin[[race]] <- data.frame(table(unlist(mm)))

    cat(paste(race,'done.\n'))
}

mm <- do.call(rbind,majmin)
mm$race <- substr(rownames(mm),1,2)
rownames(mm) <- NULL

mm <- mm %>%
    group_by(race) %>%
    mutate(sum=sum(Freq),
           Frequency=Freq/sum) %>%
    rename(Allele=Var1) %>% data.frame

tt <- data.frame(t(matrix(unlist(strsplit(as.character(mm$Allele),"\\s")),nrow=2)))
names(tt) <- c('allele','imputation')
mm <- cbind(mm,tt)

g <- ggplot(mm)+
    geom_bar(aes(allele,Frequency*100,fill=race,alpha=imputation),stat='identity',position='dodge',color='black')+
    facet_wrap(.~race)+
    scale_fill_viridis(discrete = T)+
    scale_alpha_discrete(range = c(.4, .9))+
    labs(x='Alleles',y='Percent of SNPs',alpha='Imputation',fill='Accession')
g

saveplot(g,'imputation-error-minor')
