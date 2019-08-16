## get haplotype data computed in hapex
dat <- do.call(rbind,lapply(races, function(race){
    data.frame(fread(paste0("output/hapex/hapseq_600k_v4_",race,".vcf.txt"),data.table=F),race)}))
dat$type <- sapply(as.character(dat$race), function(x) unlist(strsplit(x, "_"))[2])
dat$pop <- sapply(as.character(dat$race), function(x) unlist(strsplit(x, "_"))[1])
dat <- dat %>% mutate(hapname=paste0(chrompos, ':', seq))

## identify major haplotype in LR and look it up in DH
hapdi <- data.frame()
for (ra in racess){
    lr <- filter(dat, pop==ra, type=='LR', haps>1) # only windows with mutliple haplotypes
    lr <- lr %>% group_by(chrompos) %>% slice(which.max(freq)) %>% data.frame() # get major haplotype
    dh <- filter(dat, pop==ra, type=='DH', chrompos %in% lr$chrompos) # remove windows with only 1 hap in LR
    dh.inherit <- filter(dh, hapname %in% lr$hapname) # major haplotypes from LR that are present in DH
    dh.fixed   <- filter(dh.inherit, freq==1)
    dh.segreg  <- filter(dh.inherit, freq<1)
    dh.lost    <- filter(lr, !hapname %in% dh.inherit$hapname) %>%
        mutate(type="DH", haps=NA, count=0, freq=0, hapdiv=NA) 
    dh.all     <- rbind(mutate(dh.fixed,  majorhap='fixed'),
                        mutate(dh.segreg, majorhap='segregating'),
                        mutate(dh.lost,   majorhap='lost'))
    combined <- inner_join(lr, dh.all[,c('chrompos','freq','hapdiv','count','majorhap','haps')],
                           by='chrompos', suffix = c('.LR','.DH'))
    combined$freqdiff <- combined$freq.DH - combined$freq.LR # a LR/DH combined df
    hapdi <- rbind(hapdi, combined)                          # output
}

## add small number correction for diversity
hapdi <- hapdi %>% mutate(N.LR=ifelse(race=='BU_LR',22,
                     ifelse(race=='GB_LR',46,23)),
                N.DH=ifelse(race=='BU_LR',36,
                     ifelse(race=='GB_LR',59,
                     ifelse(race=='RT_LR',44,
                     ifelse(race=='SC_LR',58,
                     ifelse(race=='SF_LR',69,NA))))),
                hapdiv.DH=ifelse(is.na(hapdiv.DH),NA,hapdiv.DH*(N.DH/(N.DH-1))),
                hapdiv.LR=hapdiv.LR*(N.LR/(N.LR-1)))

## filter out short haplotypes
hapdi <- filter(hapdi, snps>5) 
