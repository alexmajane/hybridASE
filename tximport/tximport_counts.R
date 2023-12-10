libary(tximport)
libary(DESeq2)

# create sample reference table 
samples <- data.frame(
      c(rep("parent", 6), rep("hybrid", 6)),
      c(rep("mel", 3), rep("sim", 3), rep("mel", 3), rep("sim", 3)),
      c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3)),
      c("P_mel_1", "P_mel_2", "P_mel_3", 
        "P_sim_1", "P_sim_2", "P_sim_3",
        "F1_mel_1", "F1_mel_2", "F1_mel_3", 
        "F1_sim_1", "F1_sim_2", "F1_sim_3")
      ) %>% setNames(c("gen", "species", "cond", "sample"))


# list of salmon quant files
files_mel <- file.path("salmon", samples$sample[c(1:3,7:9)], "quant.sf")
names(files_mel) <- samples$sample[c(1:3,7:9)]
files_sim <- file.path("salmon", samples$sample[c(4:6,10:12)], "quant.sf")
names(files_sim) <- samples$sample[c(4:6,10:12)]

# import tx2gene references
tx2gene_mel <- fread("tximport/tx2gene_mel", header = F) %>% setNames(c("txname", "geneid"))
tx2gene_sim <- fread("tximport/tx2gene_sim", header = F) %>% setNames(c("txname", "geneid"))

# run tximport
library(tximport)
txi_mel <- tximport(files_mel, type = "salmon", tx2gene = tx2gene_mel)
txi_sim <- tximport(files_sim, type = "salmon", tx2gene = tx2gene_sim)

##########################################
# ortholog subsetting + name translation #
##########################################

library(data.table)
simv2v3 <- fread("tximport/simV3/Dsim_FBgn_to_PrinV3_MS_orthos", header = F) %>% setNames(c("FBgn", "symbol", "v3"))
mel_sim_ortho <- fread("tximport/simV3/ortho_MS_1to1_v3filtered.tsv", header = F) %>% setNames(c("mel", "symbol", "sim"))

txi_mel$abundance <- subset(txi_mel$abundance, rownames(txi_mel$abundance) %in% mel_sim_ortho$mel)
txi_mel$counts <-    subset(txi_mel$counts, rownames(txi_mel$counts) %in% mel_sim_ortho$mel)
txi_mel$length <-    subset(txi_mel$length, rownames(txi_mel$length) %in% mel_sim_ortho$mel)

txi_sim$abundance <- subset(txi_sim$abundance, rownames(txi_sim$abundance) %in% simv2v3$v3)
txi_sim$counts <-    subset(txi_sim$counts,    rownames(txi_sim$counts) %in%    simv2v3$v3)
txi_sim$length <-    subset(txi_sim$length,    rownames(txi_sim$length) %in%    simv2v3$v3)

## the missing genes
missing_mel <- mel_sim_ortho[!mel_sim_ortho$mel %in% rownames(txi_mel$abundance),]
write.table(missing_mel, "missing_mel.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
# list of 76 genes
simv2v3$v3[!simv2v3$v3 %in% rownames(txi_sim$abundance)]
# one gene, LOC27207718

####################################
## format data to use with DESeq2 ##
####################################

# extract estimated count matrices
counts_mel <- data.frame(txi_mel$counts)
counts_sim <- data.frame(txi_sim$counts)

# translate v3 names to sim FBgns and then to orthologous mel FBgn
row.names(counts_sim) <- mel_sim_ortho$mel[match(
            simv2v3$FBgn[match(rownames(counts_sim), simv2v3$v3)],
            mel_sim_ortho$sim
          )]

# join mel and sim, filtering to present orthologs
counts_all <- merge(counts_mel, counts_sim, by = "row.names")

# convert back to a matrix
row.names(counts_all) <- counts_all$Row.names      
counts_all <- as.matrix(counts_all[,-1])

# export
write.table(counts_all, "hybrid_counts.tsv", quote = F, sep = "\t", row.names = T, col.names = T)
