setwd("~/Documents/Begun/hybrid/sequence_analysis/")

mel_index <- read.table("mel_ortho_CDS.tsv")
sim_index <- read.table("sim_ortho_CDS.tsv")

mel_index_sort <- mel_index[order(mel_index$V3, -(mel_index$V2)), ]
mel_longest_ORF <- mel_index_sort[!duplicated(mel_index_sort$V3), ]

sim_index_sort <- sim_index[order(sim_index$V3, -(sim_index$V2)), ]
sim_longest_ORF <- sim_index_sort[!duplicated(sim_index_sort$V3), ]

ortho <- read.table("ortho_MS_1to1_v3filtered.tsv")

unique_mel <- mel_longest_ORF$V3[!ortho$V3[match(mel_longest_ORF$V3, ortho$V1)] %in% sim_longest_ORF$V3]
unique_sim <- sim_longest_ORF$V3[!ortho$V1[match(sim_longest_ORF$V3, ortho$V3)] %in% mel_longest_ORF$V3]

# export the genes whos name in mel changed in v6.41
write.table(subset(ortho, V3 %in% unique_sim), file = "updated_FBgn_6.41",  col.names = FALSE, row.names = FALSE, sep = "\t", quote = F)

# just remve the 24 sim genes man none of them are interesting for AG biology anyways
sim_longest_ORF %<>% subset(!V3 %in% unique_sim)

write.table(mel_longest_ORF, file = "mel_longest_ORFs.tsv", col.names = FALSE, row.names = FALSE, sep = "\t", quote = F)
write.table(sim_longest_ORF, file = "sim_longest_ORFs.tsv", col.names = FALSE, row.names = FALSE, sep = "\t", quote = F)

