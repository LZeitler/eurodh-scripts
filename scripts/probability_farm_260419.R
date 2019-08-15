setwd('~/ma/data/dh/')
source('~/ma/r/paperplot/source_me.R')

## calculate counts using plink wrapper
## plink <- 'plink'                        # plink command for local
plink <- '~/ma/plink/plink'             # plink dir for remote

## races <- c("GB_LR","RT_LR","SF_LR","GB_DH","RT_DH","SF_DH","WA_LR","WA_DH")
## racess <- c("GB","RT","SF",'WA')

for (race in races){
    command <- paste0(plink, ' --vcf 260419_step6_red_', race,
                      '.vcf --freq counts --out ~/ma/data/dh/freq/260419_step6_red_', race)
    system(command)
}

## get positions of SNPs of V4
vcf <- fread('260419_step6_red_SC_DH.vcf', data.table=F, select=c(1:3))

## read in counts generated with plink
## (from original data from Melchinger et al. 2017 supplements) 
counts <- data.frame()
for (race in racess){
    dhb <- fread(paste0("freq/260419_step6_red_", race, "_DH.frq.counts"), data.table=F)
    lrb <- fread(paste0("freq/260419_step6_red_", race, "_LR.frq.counts"), data.table=F)
    temp <- rbind(data.frame(dhb, Type = "DH"), data.frame(lrb, Type = "LR"))
    temp <- mutate(temp, Size = (C1+C2)/2+G0, Race = race)
    counts <- rbind(counts, temp)
}

gt <- unique(filter(counts,Type=="DH")[,c('Race','Size')]) %>% mutate(Size=Size)   # chromosome count DH
la <- unique(filter(counts,Type=="LR")[,c('Race','Size')]) %>% mutate(Size=Size*2) # chromosome count LR

## calculate joint probabilities and ancestral frequencies
pro <- list()                           # probabilities
freq <- list()                          # frequencies and intervals
for (race in racess){
   d  <- filter(counts, Type=="DH", Race==race)$C1/2 # allele count DH
   l  <- filter(counts, Type=="LR", Race==race)$C1   # allele count LR
   sD <- filter(counts, Type=="DH", Race==race)$Size # not required for computation
   sL <- filter(counts, Type=="LR", Race==race)$Size # not required for computation
   nH <- 100                            # surviving DH set to 100
   nD <- filter(gt, Race==race)$Size    # chrom DH
   nL <- filter(la, Race==race)$Size    # chrom LR
   t <- mapply(function(d, l) {
      sum(sapply(0:nH,
                 function(i) (i/nH)^d*(1-i/nH)^(nD-d)*
                    choose(nD,d)*
                    choose(nL,l)*
                    choose(nH,i)*
                    beta(i+l+0.5,nL+nH-i-l+0.5)/beta(.5,.5)))
   }, d, l)
   pro[[race]] <- data.frame(
      prob=t,
      chr=counts$CHR,
      pos=vcf$POS,
      snp=counts$SNP
   )
   ## second part here ##
   guessp <- 1:99/100
   anc_freq <- vector()
   top <- vector()
   bot <- vector()
   for (q in 1:length(d)){
       surface <- sapply(0:nH, function(i) # list of probabilities
           sapply(guessp, function(x)
               dbinom(d[q],size=nD,prob=i/nH)*
               dbinom(i,size=nH,prob=x)*
               dbinom(l[q],size=nL,prob=x)))
       anc_freq[q] <- (which(surface ==
                             max(surface)) %% nrow(surface))/100
       anc_probs <- sapply(0:nD,function(d)
           sum(sapply(0:nH, function(i)
               dbinom(d,size=nD,prob=i/nH)*
               dbinom(i,size=nH,prob=anc_freq[q]))))
       top[q] <- max(which(cumsum(anc_probs)<0.975))-1
       bot[q] <- min(which(cumsum(anc_probs)>0.025))-1
   }
   freq[[race]] <- data.frame(
       chr=vcf$`#CHROM`,
       pos=vcf$POS,
       snp=vcf$ID,
       anc_freq, top, bot, d, l, sD, sL
   )
}

## a safe way of writing a data.frame from a list
prodf  <- data.frame()
freqdf <- data.frame()
for (i in 1:length(racess)){
   t <- pro[[i]]
   t <- t[!duplicated(t$snp),] 
   t$race <- names(pro)[i]
   prodf <- rbind(prodf,               
                  t)
   t <- freq[[i]]
   t <- t[!duplicated(t$snp),] 
   t$race <- names(freq)[i]
   freqdf <- rbind(freqdf,             
                   t)
}

fwrite(prodf,       'prob/pro_260419.txt', sep = '\t', quote = F)
fwrite(freqdf,      'prob/freq_260419.txt', sep = '\t', quote = F)

