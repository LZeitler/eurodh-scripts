setwd('~/ma/data/dh/')
source('~/ma/r/paperplot/source_me.R')

#### this is also in fig2.R

## get probability and frequencies output from probability_farm_v4.R
freq <- fread('prob/freq_260419.txt', data.table = F)
freq <- filter(freq, !top==-Inf)
freq$in_CI <- freq$d >= freq$bot & freq$d <= freq$top
freq$outl <- !freq$in_CI


## update positions to v4
aff <- fread('../mayer/mayer_raw.txt', data.table=F, select=1:5)
freq4 <- inner_join(freq, aff[,c("V1","pos_v4")], by=c("snp"="V1")) %>%
    mutate(pos=pos_v4) %>% arrange(race,chr,pos) %>% select(-pos_v4)

source('~/ma/r/interpolate_map.R')

gerp <- fread('../map/gerpV4.bed', data.table = F, select = c(1,2,4,5,6,7))
vcf <- fread('260419_step6_v4_BU_DH.vcf', data.table = F, select = c(3:5))

freq4 <- filter(freq4) %>% mutate(chr=as.character(chr)) %>%
    mutate(test=ifelse(d>top, 'increase', ifelse(d<bot, 'decrease', 'in CI')),lr_freq=l/(2*sL)) %>%
    mutate(test1 = ifelse(test=='in CI', 'no outlier', 'outlier'))
freq4 <- inner_join(freq4,vcf,by=c('snp'='ID'))

gerpo <- full_join(freq4,select(gerp, chr=`#chr`, pos=start, id, gerp=correctedGERP,
                                gREF=B73Allele, gALT=derivedAllele)) %>%
    filter(chr %in% 1:10, gerp>=0) %>%
    arrange(race, chr, pos) %>%
    mutate(chr=as.integer(chr))

## match alleles
gerpo.f <- filter(gerpo,!is.na(REF))
gerpo.f$complement.gALT <- sapply(gerpo.f$gALT, function(x) complement(x))
gerpo.f <- gerpo.f %>% mutate(freq.gerp.lr = ifelse(ALT==gALT, lr_freq,
                                             ifelse(REF==gALT, 1-lr_freq,
                                             ifelse(ALT==complement.gALT, lr_freq,
                                             ifelse(REF==complement.gALT, 1-lr_freq, NA))))) %>%
    filter(!is.na(freq.gerp.lr))

gerpo <- bind_rows(data.frame(filter(gerpo,is.na(REF)),freq.gerp.lr=NA),
                   select(gerpo.f,names(gerpo),freq.gerp.lr)) %>%
    arrange(race,chr,pos)

start <- gerpo %>% group_by(chr) %>% summarise(start=min(pos))
gerpo <- left_join(gerpo, start, by='chr')

gerpo$cm  <- mapply(function(chr, start, pos) gendist(pos1 = start, pos2 = pos, chrom = chr),
                    gerpo$chr, gerpo$start, gerpo$pos)
gerpo <- filter(gerpo, !is.na(cm))    
gerpo$cut.gerp.lr <- cut(gerpo$freq.gerp.lr,breaks = seq(0,1,.1), include.lowest = T)

data.frame2 <- function(...){
    df <- try(data.frame(...),silent=T)
    if (class(df)!="data.frame"){
        cat('data.frame was not created!\n')
        df <- data.frame()
    }
    return(df)
}

## get the sum of gerpscores 1 cM distant from each outlier
## system('rm genespace_260419v4_bootstrapped_recessive-6.txt')
## system('rm genespace_260419v4_bootstrapped_additive-6.txt')

for (r in races){
    for (ch in 1:10){

        vcff <- fread(paste0('260419_step6_v4_',r,'.vcf'), data.table = F)
        inds <- names(vcff)[10:length(names(vcff))]
        ty <- substr(inds[1],1,2)
        ra <- substr(inds[1],4,5)        
        
        ## compute SUMS for outlier/non-outlier and track SNPs on the chip, to later match with
        ## genotype in individuals
        outls <- filter(gerpo, chr==ch, race==ra, test1=='outlier')
        noutls <- filter(gerpo, chr==ch, race==ra, test1=='no outlier')
        dsamp <- lapply(levels(gerpo$cut.gerp.lr), function(x){
            if (nrow(filter(noutls,cut.gerp.lr==x))>=nrow(filter(outls,cut.gerp.lr==x))){
                tn <- sample_n(filter(noutls,cut.gerp.lr==x),
                               size=nrow(filter(outls,cut.gerp.lr==x)),replace=F)
                to <- filter(outls,cut.gerp.lr==x)
            } else {
                to <- sample_n(filter(outls,cut.gerp.lr==x),
                               size=nrow(filter(noutls,cut.gerp.lr==x)),replace=F)
                tn <- filter(noutls,cut.gerp.lr==x)
            }
            rbind(to,tn)
        })
        dsamp <- do.call(rbind, dsamp)
        noutls <- filter(dsamp,outl==F)
        outls <- filter(dsamp,outl==T)
        all <- filter(gerpo, chr==ch, race==ra | is.na(race)) %>% arrange(pos)
        lo <- outls$cm-.5               # cm limit of lower bound
        hi <- outls$cm+.5               # cm limit of upper bound
        nlo <- noutls$cm-.5             # cm limit of lower bound (non outlier)
        nhi <- noutls$cm+.5             # cm limit of upper bound (non outlier)
        sumsl <- mapply(function(l,h){          # for non outlier
            il <- which.min(abs(all$cm-l)) # index of lower bound
            ih <- which.min(abs(all$cm-h)) # index of higher bound
            c(sum(all[il:ih,'gerp']),mean(all[il:ih,'gerp']),length(all[il:ih,'gerp']))
        }, lo, hi)
        nsumsl <- mapply(function(l,h){          # for non outlier
            il <- which.min(abs(all$cm-l)) # index of lower bound
            ih <- which.min(abs(all$cm-h)) # index of higher bound
            c(sum(all[il:ih,'gerp']),mean(all[il:ih,'gerp']),length(all[il:ih,'gerp']))
        }, nlo, nhi)
        sums <- data.frame(outls$snp,t(sumsl),stringsAsFactors=F)
        names(sums) <- c("snp","sum","mean","count")
        nsums <- data.frame(noutls$snp,t(nsumsl),stringsAsFactors=F)
        names(nsums) <- c("snp","sum","mean","count")
        ## sumdf1 <- rbind(sumdf1,         # all sums 
                        ## data.frame(sum=sums$sum, chr=ch, race=ra, type=ty, outl='outlier'),
                        ## data.frame(sum=nsums$sum, chr=ch, race=ra, type=ty, outl='no outlier'))
        gt <- filter(vcff[,c('ID',inds)],ID%in%c(sums$snp,nsums$snp)) # snps that are outlier/non-outl

        ## parse through individuals 
        for (ind in inds){

            gerp.here <- inner_join(filter(gerpo,race==ra,chr==ch), # contains all SNP for
                                        # outlier non outlier for popuation and chromosome
                                    data.frame(snp=gt$ID,
                                               hap1=as.integer(substr(gt[,ind],0,1)),
                                               hap2=as.integer(substr(gt[,ind],3,3)),
                                               stringsAsFactors=F),
                                    by = "snp") %>%
                mutate(which.allele=ifelse(gALT==ALT,1,
                                    ifelse(gALT==REF,0,
                                    ifelse(is.na(gALT),NA,2)))) %>%
                select(chr,snp,race,outl,hap1,hap2,which.allele)
            gerp.add.outl <- gerp.here %>%           # for additive per individual
                filter(which.allele!=2, outl==T) %>% # allele 1 or allele 2 match gerp allele
                filter(which.allele==hap1 | which.allele==hap2) %>% # hom-ref snps are filtered out
                mutate(modifier=ifelse(hap1==hap2,2,1)) # if homozygous for gerp allele, sum*2
                                        # modifier
            gerp.add.noutl <- gerp.here %>% filter(outl==F)
            gerp.add.noutl <- sample_n(tbl=gerp.add.noutl,size=nrow(gerp.add.outl)) %>%
                mutate(modifier=ifelse(hap1==hap2,2,1))
            gerp.rec.outl <- gerp.here %>%           # for recessive per individual
                filter(which.allele!=2, outl==T) %>% # allele 1 and allele 2 match gerp allele
                filter(which.allele==hap1 & which.allele==hap2)
            gerp.rec.noutl <- gerp.here %>% filter(outl==F)
            gerp.rec.noutl <- sample_n(tbl=gerp.rec.noutl,size=nrow(gerp.rec.outl))

            sums.a <- inner_join(gerp.add.outl,sums,by='snp') %>%
                mutate(sum.a=sum*modifier)
            nsums.a <- inner_join(gerp.add.noutl,nsums,by='snp') %>%
                mutate(sum.a=sum*modifier)
            sums.r <- inner_join(gerp.rec.outl,sums,by='snp')
            nsums.r <- inner_join(gerp.rec.noutl,nsums,by='snp')

            dfs.s <- data.frame2(
                sum=sums.r$sum,
                snp=sums.r$snp,
                mean=sums.r$mean,
                count=sums.r$count,
                ind, type=ty,
                race=ra,
                outl='outlier')
            dfs.n <- data.frame2(
                sum=nsums.r$sum,
                snp=nsums.r$snp,
                mean=nsums.r$mean,
                count=nsums.r$count,
                ind, type=ty,
                race=ra,
                outl='no outlier')
            if (nrow(dfs.s)==0 | nrow(dfs.n)==0) next # if no alleles to analyse go to next indiv

            fwrite(                     # write recessive
                rbind(dfs.s,dfs.n),
                'genespace_260419v4_recessive.txt',
                sep = '\t', append = T)

            dfs.s <- data.frame2(
                sum=sums.a$sum.a,
                snp=sums.a$snp,
                mean=sums.a$mean,
                count=sums.a$count,
                ind, type=ty,
                race=ra,
                outl='outlier')
            dfs.n <- data.frame2(
                sum=nsums.a$sum.a,
                snp=nsums.a$snp,
                mean=nsums.a$mean,
                count=nsums.a$count,
                ind, type=ty,
                race=ra,
                outl='no outlier')
            if (nrow(dfs.s)==0 | nrow(dfs.n)==0) next # if no alleles to analyse go to next indiv

            fwrite(                     # write additive
                rbind(dfs.s,dfs.n),
                'genespace_260419v4_additive.txt',
                sep = '\t', append = T)
            
            ## status report
            cat(paste('Status: Ind',round(match(ind,inds)/length(inds),2)*100,
                      '%, Chr',ch,', Pop',r,'\n'))
            
        }
    }
}
