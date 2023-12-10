intersect_500 <- read.table("~/Documents/Begun/hybrid/sequence_analysis/sim.500.intersect.formatted.tsv")
intersect_750 <- read.table("~/Documents/Begun/hybrid/sequence_analysis/sim.750.intersect.formatted.tsv")
intersect_1000 <- read.table("~/Documents/Begun/hybrid/sequence_analysis/sim.1000.intersect.formatted.tsv")

intersect_500 <- intersect_500[,c(1,4,5,7,10,13,16,17,19,22,25)]
intersect_750 <- intersect_750[,c(1,4,5,7,10,13,16,17,19,22,25)]
intersect_1000 <- intersect_1000[,c(1,4,5,7,10,13,16,17,19,22,25)]

temp <- lapply(list(intersect_500,
            intersect_750,
            intersect_1000), function(x) 
                                setNames(x, c("A_chr",
                                              "A_start",
                                              "A_end",
                                              "A_strand",
                                              "A_gene",
                                              "B_chr",
                                              "B_start",
                                              "B_end",
                                              "B_strand",
                                              "B_gene",
                                              "overlap")))

intersect_500 <- temp[[1]]
intersect_750 <- temp[[2]]
intersect_1000 <- temp[[3]]

rm(temp)

intersect_500 %<>%  subset(A_gene != B_gene)
intersect_750 %<>%  subset(A_gene != B_gene)
intersect_1000 %<>% subset(A_gene != B_gene)

dim(subset(intersect_500, A_chr != B_chr)) # there are no genes on different mel chromosomes that uniquely map to overlapping spot in simulans

# the lengths of overlap as fraction of match lengths?
intersect_500 %<>% mutate(A_len = A_end - A_start + 1,
                          B_len = B_end - B_start + 1,)
intersect_500 %<>% mutate(overlap_A = overlap / A_len,
                          overlap_B = overlap / B_len,)
# there are very few things that actually have a small fraction of overlap. like 7 genes with < 30%.
# the lowest amount of overlap is still 14%. 

# SO, lets just throw all these sequences away.

table(intersect_500$A_gene %in% intersect_500$B_gene)

intersect_500_genes <- unique(intersect_500$A_gene)
# 506 out of 13014 genes (3.9%)
intersect_750_genes <- unique(intersect_750$A_gene)
# 424 out of 11906 genes (3.6%)
intersect_1000_genes <- unique(intersect_1000$A_gene)
# 393 out of 10997 genes (3.6%)

table(intersect_750_genes %in% intersect_500_genes)
# FALSE  TRUE 
# 15   409 
table(intersect_1000_genes %in% intersect_500_genes)
# FALSE  TRUE 
# 16   377 
table(intersect_1000_genes %in% intersect_750_genes)
# FALSE  TRUE 
# 2   391 

# interesting a small number of the genes identified in longer sequence sets
# are not found in the smaller set. This could be because the sequences overlap in 
# a relatively small area in the longer sequence set but a bigger proportion in the smaller sequence set

write.table(intersect_500_genes, "~/Documents/Begun/hybrid/sequence_analysis/overlapping_sim_500_blacklist", col.names = F, row.names = F, quote = F)
write.table(intersect_750_genes, "~/Documents/Begun/hybrid/sequence_analysis/overlapping_sim_750_blacklist", col.names = F, row.names = F, quote = F)
write.table(intersect_1000_genes, "~/Documents/Begun/hybrid/sequence_analysis/overlapping_sim_1000_blacklist", col.names = F, row.names = F, quote = F)
