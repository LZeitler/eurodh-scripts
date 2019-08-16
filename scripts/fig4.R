source('source_me.R')

#### load data

freq <- get.freq('output/freq_260419.txt')                    # aSFS
prodf <- get.prob('output/pro_260419.txt')                    # jProb

#### prepare data

## calculate heterozygosity
hetcount <- function(prefix){
    plink <- 'plink'                        # plink command for local
    for (race in races){
        command <- paste0(plink, ' --vcf ', prefix,
                          race,'.vcf --freqx --out output/het/260419_step6_red_freqx_', race)
        system(command)
    }
    ## read in data
    counts <- data.frame()
    for (race in racess){
        dhb <- fread(paste0("output/het/260419_step6_red_freqx_", race, "_DH.frqx"), data.table=F)
        lrb <- fread(paste0("output/het/260419_step6_red_freqx_", race, "_LR.frqx"), data.table=F)
        temp <- rbind(data.frame(dhb, Type = "DH", race), data.frame(lrb, Type = "LR", race))
        counts <- rbind(counts, temp)
    }
    counts$total <- rowSums(counts[,5:10])
    counts$hetfreq <- counts$C.HET./counts$total
    
    counts <- filter(counts,Type=='LR')
    return(counts)
}

counts <- hetcount('data/260419_step6_red_v4_') # positions v4
fwrite(counts,'output/het/260419_step6_red_v4_freqx.txt',sep='\t')

## prepare SNP LD pruning
system('plink --vcf data/260419_step6_red_v4_LR.vcf --indep-pairphase 500 500 0.2 --out output/prune/prune')
prune <- as.vector(fread('output/prune/prune.prune.in',data.table=F,header=F)[,1])

## prep aSFS data
freq.het <- inner_join(freq,counts,by=c('race'='race','snp'='SNP'))
freq.het <- rbind(filter(freq.het,outl==F,snp%in%prune),filter(freq.het,outl==T)) # apply pruning
freq.het$lr_freq.bin <- round(freq.het$l/(2*freq.het$sL),1)

## prep probability data
prob.het <- inner_join(prodf,counts,by=c('race'='race','snp'='SNP'))
prob.het <- rbind(filter(prob.het,outl==F,snp%in%prune),filter(prob.het,outl==T)) # apply pruning
prob.het <- inner_join(prob.het,                                                  # get allele frequencies
                       select(freq.het,race,snp,lr_freq.bin),by=c('race'='race','snp'='snp'))

## downsample function
stratify <- function(data){
    out <- data.frame()
    for (ra in racess){
        for (bin in sort(unique(data$lr_freq.bin))){
            t <- filter(data,race==ra,outl==T,lr_freq.bin==bin)
            t1 <- filter(data,race==ra,outl==F,lr_freq.bin==bin)
            n <- min(nrow(t),nrow(t1))
            out <- bind_rows(out,
                             sample_n(filter(data,race==ra,outl==T,lr_freq.bin==bin),size=n),
                             sample_n(filter(data,race==ra,outl==F,lr_freq.bin==bin),size=n))
        }
    }
    return(out)
}

#### plotting

## prepare plots
fig.het <- function(data){
    c <- ggplot(data)
    c <- c+geom_violin(aes(race,hetfreq,fill=outl),alpha=.8)
    c <- c+stat_summary(aes(race,hetfreq,group=outl),alpha=.8,
                        position=position_dodge(.9),fun.y='mean',geom='point',shape=18,size=3)
    c <- c+scale_fill_viridis(discrete = T, begin = .2, end = .8, direction = -1,labels=c('non-outlier','outlier'))
    c <- c+labs(x='Accession',y='Heterozygote frequency',fill='SNP')
    c
}


## draw plots
sampled <- stratify(freq.het)
c1 <- fig.het(sampled)

sampled <- stratify(prob.het)
c2 <- fig.het(sampled)


## stitch together
f4 <- plot_grid(c1+theme(legend.position = "none"),
                c2+theme(legend.position = "none"),ncol=1,labels='AUTO')
f4 <- plot_grid(f4,get_legend(c1),nrow=1,rel_widths=c(.65,.1))
f4

saveplot(f4,'fig4-het-v') ## thesis fig 4
