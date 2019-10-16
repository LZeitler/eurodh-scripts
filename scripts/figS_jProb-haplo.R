source('source_me.R')

source('scripts/make_hapdiv.R')                # load haplotype data

prodf <- get.prob('output/pro_260419.txt') %>% # load jProb data
    mutate(pop=race,
           pos.start=pos,
           pos.end=pos)

hap <- hapdi %>%
    select(chr,pos,freq.LR,race,type,pop,majorhap,freqdiff) %>%
    mutate(pos.end=pos+49999,
           pos.start=pos) %>%
    select(-pos) %>%
    filter(type=='LR')

setDT(hap)
setDT(prodf)

setkey(hap,pop,chr,pos.start,pos.end)
setkey(prodf,pop,chr,pos.start,pos.end)

ph <- foverlaps(hap,prodf,type='any') %>% # window merge
    filter(!is.na(prob)) %>%
    filter(outl==T) %>% 
    group_by(pop,chr,i.pos.start) %>%
    mutate(meanp=mean(logprob),
           meanfreqdiff=mean(freqdiff),
           meanfreqlr=mean(freq.LR),
           numo=length(logprob)) %>%
    data.frame

a <- ggplot(filter(ph))+
    geom_boxplot(aes(majorhap,logprob,fill=majorhap),alpha=.7)+
    facet_wrap(~pop,scales = 'free',labeller = as_labeller(racelabels))+
    labs(y='-log10(p) (jProb)',x='Fate')+
    scale_fill_manual(values = c('#481567FF','red','#3CBB75FF'))+
    theme(legend.position = 'none')
a





