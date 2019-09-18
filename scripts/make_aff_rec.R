source('source_me.R')
source('interpolate_map.R')


interv <- 1000000

aff <- fread('data/mayer_raw.txt', data.table=F, select=1:5)
aff <- aff %>%
    select(snp=V1,chr=chr_v4,pos=pos_v4) %>% 
    filter(chr!=0) %>%
    group_by(chr) %>%
    arrange(chr) %>% 
    mutate(start=min(pos),
           stop=max(pos)) %>%
    data.frame
windows <- do.call(rbind,lapply(1:10, function(x){
    data.frame(chr=x,
               pos=seq(filter(aff,chr==x)$start[1],filter(aff,chr==x)$stop[1],interv))
})) %>% mutate(winpos=pos)
windows$rec <- mapply(function(chr, start, pos) gendist(pos1 = start, pos2 = pos, chrom = chr),
                    windows$chr, windows$pos, windows$pos+interv)
aff <- full_join(aff,windows,by=c('chr','pos')) %>%
    arrange(chr,pos) %>%
    fill(winpos,rec,.direction='down') %>%
    filter(!is.na(snp))


fwrite(aff,'map/aff-v4-with-rec-1mbp.txt',quote=F,sep='\t')
