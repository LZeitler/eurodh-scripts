source('~/ma/r/paperplot/source_me.R')
source('~/pro/eurodh-scripts/scripts/source_me.R')
source('~/ma/r/interpolate_map.R')
setwd('~/ma/data/dh/')

set.seed(123)

freq <- get.freq('~/pro/eurodh-scripts/output/freq_260419.txt') # aSFS
aff <- fread('../map/aff-v4-with-rec-1mbp.txt', data.table=F)   # annotation
vcf <- fread('260419_step6_v4_BU_DH.vcf', data.table = F,       # 600k to get REF ALT
             select = c(1:5))
gerp <- fread('../map/gerpV4.bed', data.table = F,              # V4 GERP scores
              select = c(1,2,4,5,6,7)) %>%
    filter(`#chr` %in% 1:10, start %in% vcf$POS)  %>%
    mutate(chr=as.integer(`#chr`))


freq <- inner_join(select(freq,-chr,-pos), select(aff,snp,chr,pos,rec), by="snp") %>%
    arrange(race,chr,pos)               # update positions
freq <- inner_join(freq,vcf,by=c('snp'='ID'))

gerpo <- right_join(freq,select(gerp, chr, pos=start, id, gerp=correctedGERP,
                               gREF=B73Allele, gALT=derivedAllele)) %>%
    filter(gerp>=0, snp %in% vcf$ID) %>%
    arrange(race, chr, pos)

## match alleles
gerpo$complement.gALT <- sapply(gerpo$gALT, function(x) complement(x))
gerpo <- gerpo %>% mutate(lr_freq=l/(2*sL),
                          freq.gerp.lr = ifelse(ALT==gALT, lr_freq,
                                         ifelse(REF==gALT, 1-lr_freq,
                                         ifelse(ALT==complement.gALT, lr_freq,
                                         ifelse(REF==complement.gALT, 1-lr_freq, NA))))) %>%
    group_by(chr) %>%
    mutate(start=min(pos)) %>%
    data.frame() %>% 
    filter(!is.na(freq.gerp.lr)) %>%
    arrange(race,chr,pos)

gerpo$cm  <- mapply(function(chr, start, pos) gendist(pos1 = start, pos2 = pos, chrom = chr),
                    gerpo$chr, gerpo$start, gerpo$pos)
gerpo <- filter(gerpo, !is.na(cm))    
gerpo$cut.gerp.lr <- cut(gerpo$freq.gerp.lr,breaks=seq(0,1,.1),include.lowest=T)
gerpo$cut.rec <- cut(gerpo$rec,quantile(gerpo$rec,na.rm=T),include.lowest=T)
gerpo$cut.freq.rec <- paste(gerpo$cut.gerp.lr,gerpo$cut.rec)

sample_equal <- function(dt1,dt2,variable){
    lapply(unique(rbind(dt1,dt2)[,variable]), function(x){
        dt1n <- nrow(dt1[which(dt1[,variable]==x),])
        dt2n <- nrow(dt2[which(dt2[,variable]==x),])
        ## if (dt1n==0|dt2n==0) next
        if (dt1n>=dt2n){
            tn <- sample_n(dt1[which(dt1[,variable]==x),],
                           size=dt2n,replace=F)
            to <- dt2[which(dt2[,variable]==x),]
        } else {
            to <- sample_n(dt2[which(dt2[,variable]==x),],
                           size=dt1n,replace=F)
            tn <- dt1[which(dt1[,variable]==x),]
        }
        rbind(to,tn)
    })
}

for (r in races){
    for (ch in 1:10){

        vcff <- fread(paste0('260419_step6_v4_',r,'.vcf'), data.table = F)
        inds <- names(vcff)[10:length(names(vcff))]
        ty <- substr(inds[1],1,2)
        ra <- substr(inds[1],4,5)        
        outls <- filter(gerpo, chr==ch, race==ra, outl==T)
        noutls <- filter(gerpo, chr==ch, race==ra, outl==F)

        ## downsampling (comment out if not needed)
        dsamp <- do.call(rbind,sample_equal(noutls,outls,'cut.freq.rec'))
        ## noutlss <- filter(dsamp,outl==F)
        ## outlss <- filter(dsamp,outl==T)
        means50[[paste(r,ch)]] <- dsamp # store for direct comparison later

        ## genotypes and match alleles
        gt <- filter(vcff[,c('ID',inds)],ID%in%c(noutls$snp,outls$snp)) %>%
            rename(snp=ID)
        gerp.local <- filter(gerpo,race==ra,chr==ch) %>%
            mutate(which.allele=ifelse(gALT==ALT,1,
                                ifelse(gALT==REF,0,NA))) %>%
            filter(!is.na(which.allele))
        
        ## parse through individuals 
        for (ind in inds){
            
            gerp.here <- inner_join(gerp.local,
                                    data.frame(snp=gt$snp,
                                               hap1=as.integer(substr(gt[,ind],0,1)),
                                               hap2=as.integer(substr(gt[,ind],3,3)),
                                               stringsAsFactors=F),
                                    by = "snp") %>%
                select(chr,pos,snp,race,outl,hap1,hap2,which.allele,gerp,cm,cut.freq.rec) %>%
                mutate(gerp1=ifelse(hap1==which.allele,gerp,0),
                       gerp2=ifelse(hap2==which.allele,gerp,0),
                       gerpa=gerp1+gerp2,
                       gerpr=ifelse(gerp1!=0 & gerp2!=0,gerp1,0)) %>% 
                arrange(pos)
            
            ## calculate sums (cM windows)
            all <- gerp.here
            lo <- dsamp$cm-.5               # cm limit of lower bound
            hi <- lo+1                      # cm limit of upper bound
            sumsl <- mapply(function(l,h){
                il <- which.min(abs(all$cm-l)) # index of lower bound
                ih <- which.min(abs(all$cm-h)) # index of higher bound
                c(sum(all[il:ih,'gerpa']),mean(all[il:ih,'gerpa']),length(all[il:ih,'gerpa']),
                  sum(all[il:ih,'gerpr']),mean(all[il:ih,'gerpr']),length(all[il:ih,'gerpr']))
            }, lo, hi)
            sums <- data.frame(dsamp$snp,t(sumsl),dsamp$cut.freq.rec,dsamp$outl,stringsAsFactors=F)
            names(sums) <- c("snp","sum.a","mean.a","count.a","sum.r","mean.r","count.r",'cut.freq.rec','outl')
            sums$outl <- ifelse(sums$outl,'outlier','no outlier')

            sums.r <- sums %>%
                select(-sum.a,-mean.a,-count.a,sum=sum.r,mean=mean.r,count=count.r,cut.freq.rec,outl) %>%
                mutate(ind=ind,type=ty,race=ra)
            
            fwrite(                     # write recessive
                sums.r,
                'genespace_260419v4_600_recessive_seed123.txt',
                sep = '\t', append = T)

            sums.a <- sums %>%
                select(-sum.r,-mean.r,-count.r,sum=sum.a,mean=mean.a,count=count.a,cut.freq.rec,outl) %>%
                mutate(ind=ind,type=ty,race=ra)
            
            fwrite(                     # write recessive
                sums.a,
                'genespace_260419v4_600_additive_seed123.txt',
                sep = '\t', append = T)

            ## status report
            cat(paste('Status: Ind',round(match(ind,inds)/length(inds),2)*100,
                      '%, Chr',ch,', Pop',r,'\n'))

        }
    }
}

