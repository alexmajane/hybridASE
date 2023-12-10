library(data.table)
library(magrittr)
library(dplyr)
library(plyr)
library(tibble)
library(ggplot2)
setwd("~/Documents/Begun/hybrid/bad_align_counts/")

counts_missassigned_mel <- cbind(fread("mel_1_missassigned.tab"),
                                 fread("mel_2_missassigned.tab")[,2],
                                 fread("mel_2_missassigned.tab")[,2]) %>%
                                 setNames(c("gene", "m1", "m2", "m3")) %>%
                                 mutate(m_avg = (m1 + m2 + m3)/3)
counts_missassigned_sim <- cbind(fread("sim_1_missassigned.tab"),
                                 fread("sim_2_missassigned.tab")[,2],
                                 fread("sim_2_missassigned.tab")[,2]) %>% 
                                 setNames(c("gene", "m1", "m2", "m3")) %>%
                                 mutate(m_avg = (m1 + m2 + m3)/3)
counts_missassigned_unique_mel <- cbind(fread("mel_1_sim_unique.tab"),
                                  fread("mel_2_sim_unique.tab")[,2],
                                  fread("mel_2_sim_unique.tab")[,2]) %>%
                                  setNames(c("gene", "u1", "u2", "u3")) %>%
                                  mutate(u_avg = (u1 + u2 + u3)/3)
counts_missassigned_unique_sim <- cbind(fread("sim_1_mel_unique.tab"),
                                  fread("sim_2_mel_unique.tab")[,2],
                                  fread("sim_2_mel_unique.tab")[,2]) %>% 
                                  setNames(c("gene", "r1", "r2", "r3")) %>%
                                  mutate(u_avg = (r1 + r2 + r3)/3)

### change sim names to mel and subset to the 1-to-1 orthos
simv2v3 <- fread("../tximport/simV3/Dsim_FBgn_to_PrinV3_MS_orthos", header = F) %>% setNames(c("FBgn", "symbol", "v3"))
mel_sim_ortho <- fread("../tximport/simV3/ortho_MS_1to1_v3filtered.tsv", header = F) %>% setNames(c("mel", "symbol", "sim"))

counts_missassigned_mel$gene <- mel_sim_ortho$mel[match(
  simv2v3$FBgn[match(counts_missassigned_mel$gene, simv2v3$v3)],
  mel_sim_ortho$sim
)]
counts_missassigned_unique_mel$gene <- mel_sim_ortho$mel[match(
  simv2v3$FBgn[match(counts_missassigned_unique_mel$gene, simv2v3$v3)],
  mel_sim_ortho$sim
)]
# merge with accepted overall counts so that we can look at the 
# relative size effect of missassigned counts vs correctly assigned

# get the parent counts
counts <- read.table("~/Documents/Begun/hybrid/hybrid_counts.tsv", )
counts %<>% round()
counts %<>% rownames_to_column("gene")
counts <- counts[,c(1:4,8:10)] # remove hybrid counts
counts %<>% mutate(pmel_avg = (P_mel_1 + P_mel_2 + P_mel_3) / 3)
counts %<>% mutate(psim_avg = (P_sim_1 + P_sim_2 + P_sim_3) / 3)

# join counts together
bad_counts_mel <- join_all(list(counts[,c(1:4,8)],
                                counts_missassigned_mel,
                                counts_missassigned_unique_mel), 
                           by = "gene")
bad_counts_sim <- join_all(list(counts[,c(1,5:7,9)],
                                counts_missassigned_sim,
                                counts_missassigned_unique_sim), 
                           by = "gene")
# total misassigned counts
bad_counts_mel %<>% mutate(tot_mis_avg = u_avg + m_avg) 
bad_counts_sim %<>% mutate(tot_mis_avg = u_avg + m_avg) 

# fractional change in logfold counts with misassigned reads
bad_counts_mel %<>% mutate(change = 
              (log(pmel_avg + tot_mis_avg + 1) - log(pmel_avg + 1)) / 
                      log(pmel_avg + tot_mis_avg + 1) 
              ) 
dim(subset(bad_counts_mel, change == 0)) # 7110 of the genes have no missasigned counts
bad_genes_mel <- subset(bad_counts_mel, change > 0.025)$gene # 307 genes

bad_counts_sim %<>% mutate(change = 
                             (log(psim_avg + tot_mis_avg + 1) - log(psim_avg + 1)) / 
                                        log(psim_avg + tot_mis_avg + 1) 
) 
dim(subset(bad_counts_sim, change == 0)) # 8568 of the genes have no missasigned counts
bad_genes_sim <- subset(bad_counts_sim, change > 0.025)$gene # 97 genes 

bad_genes_all <- unique(c(bad_genes_mel, bad_genes_sim)) # get all unique genes
# so we have a list of 382 unique genes that should be removed from consideration

write.table(bad_genes_all, "bad_genes.tsv", quote = F, row.names = F, col.names = F, sep = "/t")

# gene families
fam <- fread("gene_group_data_fb_2021_06.tsv")
fam[fam==""] <- NA

# join gene families to list of genes where misalignment influences
# log(counts) greater than threshold
bad_genes_all %<>% data.frame() %>% setNames("gene")
bad_genes_fam <- right_join(fam, bad_genes_all,
                           by = c("Group_member_FB_gene_id"="gene"))

# create reference terms table
fam.list <- as.list(split(fam[,1], fam$Group_member_FB_gene_id))
# create GO object
library(topGO)

expgenes <- subset(counts, rowSums(counts[,c(2:7)]) > 3)$gene
geneList <- factor(as.integer(expgenes %in% bad_genes_all$gene))
names(geneList) <- expgenes

GO <- new("topGOdata",
          allGenes = geneList,
          annot = fam.list, ontology = "MF")
        
