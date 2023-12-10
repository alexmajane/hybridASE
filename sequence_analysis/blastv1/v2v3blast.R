blastv2 <- read.table("blast/dmel_500_dsimv2.tsv")
blastv3 <- read.table("blast/dmel_500_dsimv3.tsv")

# remove all but the top hit
blastv2.top <- blastv2[!duplicated(blastv2$V1),]
blastv3.top <- blastv3[!duplicated(blastv3$V1),]

# analyze differences in lengths and qual of top hits
blastv2.top["ref"] <- "v2"
blastv3.top["ref"] <- "v3"

blast_stats <- inner_join(blastv2.top[,c(1,3,4,12)], blastv3.top[,c(1,3,4,12)], by = "V1") %>%
                setNames(c("gene", "id.v2", "len.v2", "eval.v2", "id.v3", "len.v3", "eval.v3"))
blast_stats %<>% mutate(id.diff = id.v3 - id.v2,
                        len.diff = len.v3 - len.v2,
                        eval.diff = eval.v3 - eval.v2)

# what proportion are unchanged?
summary(blast_stats$id.diff == 0)
#    Mode   FALSE    TRUE     # 4.3%
# logical     707   16570 
summary(blast_stats$len.diff == 0)
#    Mode   FALSE    TRUE     # 3.9%
# logical     643   16634    
summary(blast_stats$eval.diff == 0)
#    Mode   FALSE    TRUE     # 4.1%
# logical     681   16596 

# distribution of positive values?
ggplot(subset(blast_stats, id.diff != 0), aes(x = id.diff)) +
  geom_density(#bins = 52, binwidth = 0.5, center = 0.75,
                 color = "mediumseagreen", fill = "transparent") + 
  geom_vline(aes(xintercept=0),
             color="grey51", linetype="solid", size=1, alpha = 0.8) +
  geom_vline(aes(xintercept=median(id.diff)),
             color="mediumseagreen", linetype="solid", size=1, alpha = 0.8) +
  theme_minimal() + xlab("v3 - v2 percent ID") 

ggplot(subset(blast_stats, len.diff != 0), aes(x = len.diff)) +
  geom_density(#bins = 52, binwidth = 0.5, center = 0.75,
    color = "mediumseagreen", fill = "transparent") + 
  geom_vline(aes(xintercept=0),
             color="grey51", linetype="solid", size=1, alpha = 0.8) +
  geom_vline(aes(xintercept=median(len.diff)),
             color="mediumseagreen", linetype="solid", size=1, alpha = 0.8) +
  theme_minimal() + xlab("v3 - v2 length")

ggplot(subset(blast_stats, eval.diff != 0), aes(x = eval.diff)) +
  geom_density(#bins = 52, binwidth = 0.5, center = 0.75,
    color = "mediumseagreen", fill = "transparent") + 
  geom_vline(aes(xintercept=0),
             color="grey51", linetype="solid", size=1, alpha = 0.8) +
  geom_vline(aes(xintercept=median(eval.diff)),
             color="mediumseagreen", linetype="solid", size=1, alpha = 0.8) +
  theme_minimal() + xlab("v3 - v2 e-value")

summary(blast_stats$id.diff)
summary(subset(blast_stats, id.diff != 0)$id.diff)
summary(blast_stats$len.diff)
summary(subset(blast_stats, len.diff != 0)$len.diff)
summary(blast_stats$eval.diff)
summary(subset(blast_stats, eval.diff != 0)$eval.diff)

###################################
## effect of loosening gap costs ##
###################################

# default gap open/gap extend of 5/2 (above) vs looser costs of 2/1
blastv2_21 <- read.table("blast/dmel_500_dsimv2_OE21.tsv")
blastv3_21 <- read.table("blast/dmel_500_dsimv3_OE21.tsv")

# summarize distribution of number of  hits per gene 
hitsv2 <- table(count(blastv2$V1)[,2])
hitsv2.21 <- table(count(blastv2_21$V1)[,2])
hitsv3 <- table(count(blastv3$V1)[,2])
hitsv3.21 <- table(count(blastv3_21$V1)[,2])

hits.df <- data.frame(do.call(cbind, lapply(list(hitsv2,
                           hitsv2.21,
                           hitsv3,
                           hitsv3.21), head, n = 10))) %>% 
                       rownames_to_column() %>%
                       setNames(c("no.hits", "v2", "v2_21", "v3", "v3_21"))
colSums(hits.df[,-1]) # almost NO change in number of things with 10 or fewer,
colSums(hits.df[1:5,-1]) # also almost no change in 5 or fewer hits.
colSums(hits.df[1:3,-1]) # modest increases in things with 3 or fewer hits
colSums(hits.df[1:3,-1]) # with the loosened gap params

hits.df %<>% mutate(diff.v2 = 100*((v2_21 - v2) / v2),
                    diff.v3 = 100*((v3_21 - v3) / v3))

## compare lengths of top hit transcripts
blastv2_21.top <- blastv2_21[!duplicated(blastv2_21$V1),]
blastv3_21.top <- blastv3_21[!duplicated(blastv3_21$V1),]

summary(blastv2.top$V4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 28.0   428.0   500.0   447.9   506.0   569.0 
summary(blastv2_21.top$V4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 28.0   495.0   501.0   476.7   512.0   621.0 
summary(blastv3.top$V4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 28.0   442.0   500.0   450.9   506.0   569.0 
summary(blastv3_21.top$V4)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 28.0   497.0   502.0   479.7   512.0   621.0 

# so we do see some improvement in candidate transcript lengths!

# other param variants
blastv3_31 <- read.table("blast/dmel_500_dsimv3_OE31.tsv")
blastv3_32 <- read.table("blast/dmel_500_dsimv3_OE32.tsv")

hitsv3.31 <- table(count(blastv3_31$V1)[,2])
hitsv3.32 <- table(count(blastv3_32$V1)[,2])

hits.df <- data.frame(do.call(cbind, lapply(list(hitsv3,
                                                 hitsv3.32,
                                                 hitsv3.31,
                                                 hitsv3.21), head, n = 10))) %>% 
                                              rownames_to_column() %>%
                                              setNames(c("no.hits", "5,2", "3,2", "3,1", "2,1"))

blastv3_31.top <- blastv3_31[!duplicated(blastv3_31$V1),]
blastv3_32.top <- blastv3_32[!duplicated(blastv3_32$V1),]

summary(blastv3_31.top$V4)
summary(blastv3_32.top$V4)

ggplot(blastv3_31.top, aes(x = V4)) + geom_density()

#### what percent of the set of orthologs have single-hit BLAST results??
hitcounts <- count(blastv3_21$V1)
singlets <- hitcounts[hitcounts$freq == 1,]$x
doublets <- hitcounts[hitcounts$freq == 2,]$x
triplets <- hitcounts[hitcounts$freq == 3,]$x

mel_sim_ortho <- fread("tximport/simV3/ortho_MS_1to1_v3filtered.tsv", header = F) %>% setNames(c("mel", "symbol", "sim"))
table(mel_sim_ortho$mel %in% singlets)
# FALSE  TRUE 
# 1247 11039 
table(mel_sim_ortho$mel %in% doublets)
# FALSE  TRUE 
# 11738   548 
table(mel_sim_ortho$mel %in% triplets)
# FALSE  TRUE 
# 12247    39 

sfps <- read.table("acp_annotation/wigby2020")
table(sfps$V3 %in% singlets)
# FALSE  TRUE 
#    37   255 

### getting set of cis genes from  hybrid_analysis_v3.R 
all_cis <- c(cis$gene, cisbytrans$gene, cisplustrans$gene)
all_cis <- synonyms$`##primary_FBid`[match(all_cis,synonyms$current_symbol)]
table(all_cis %in% singlets) # we are missing just 227 genes of 2153, acceptable filtering loss


