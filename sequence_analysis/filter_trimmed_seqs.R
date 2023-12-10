untrimmed500 <- read.table("mel.500.tsv")
trimmed500 <- read.table("mel.500.notranscriptbut3UTR.tsv")

untrimmed750 <- read.table("mel.750.tsv")
trimmed750 <- read.table("mel.750.notranscriptbut3UTR.tsv")

untrimmed1000 <- read.table("mel.1k.tsv")
trimmed1000 <- read.table("mel.1k.notranscriptbut3UTR.tsv")
trimmed500 %<>% mutate(len = V5 - V4)
trimmed750 %<>% mutate(len = V5 - V4)
trimmed1000 %<>% mutate(len = V5 - V4)

trimmed500.filtered  <- subset(trimmed500,  len > 333)
trimmed750.filtered  <- subset(trimmed750,  len > 500)
trimmed1000.filtered <- subset(trimmed1000, len > 666)

trimmed500.filtered <- left_join(trimmed500.filtered, untrimmed500[,c(4,5,9)], by = "V9") %>%
                          setNames(c("chr",
                                     "src",
                                     "feature",
                                     "start.trimmed",
                                     "end.trimmed",
                                     "score",
                                     "strand",
                                     "phase",
                                     "attr",
                                     "trimmed.len",
                                     "start.untrimmed",
                                     "end.untrimmed"))
trimmed750.filtered <- left_join(trimmed750.filtered, untrimmed750[,c(4,5,9)], by = "V9") %>%
                          setNames(c("chr",
                                     "src",
                                     "feature",
                                     "start.trimmed",
                                     "end.trimmed",
                                     "score",
                                     "strand",
                                     "phase",
                                     "attr",
                                     "trimmed.len",
                                     "start.untrimmed",
                                     "end.untrimmed"))
trimmed1000.filtered <- left_join(trimmed1000.filtered, untrimmed1000[,c(4,5,9)], by = "V9") %>%
                          setNames(c("chr",
                                     "src",
                                     "feature",
                                     "start.trimmed",
                                     "end.trimmed",
                                     "score",
                                     "strand",
                                     "phase",
                                     "attr",
                                     "trimmed.len",
                                     "start.untrimmed",
                                     "end.untrimmed"))

# now futher subset. from our filtered set we want things where the upstream
# sequence still meets the 5' boundary of the gene. So, if + stranded, the end 
# doesn't shift, and if - stranded, the start doesnt't shift
trimmed500.filtered %<>% subset((strand == "+" & end.trimmed == end.untrimmed) |
                                (strand == "-" & start.trimmed == start.untrimmed)
                                )
            # 15059 of our 15148 (99.4%) sequences fill this criteria
trimmed750.filtered %<>% subset((strand == "+" & end.trimmed == end.untrimmed) |
                                (strand == "-" & start.trimmed == start.untrimmed)
                                )
            # 14160 / 14275
trimmed1000.filtered %<>% subset((strand == "+" & end.trimmed == end.untrimmed) |
                                (strand == "-" & start.trimmed == start.untrimmed)
                                )
            # 13517 / 13665

# create bed formatted tables
trimmed500.filtered.bed <- trimmed500.filtered[,c(1,4,5,9,6,7,2,3,8)]
# bed is 0-indexed so subtract 1 from the start interval
trimmed500.filtered.bed$start.trimmed <- trimmed500.filtered.bed$start.trimmed - 1
trimmed750.filtered.bed <- trimmed750.filtered[,c(1,4,5,9,6,7,2,3,8)]
trimmed750.filtered.bed$start.trimmed <- trimmed750.filtered.bed$start.trimmed - 1
trimmed1000.filtered.bed <- trimmed1000.filtered[,c(1,4,5,9,6,7,2,3,8)]
trimmed1000.filtered.bed$start.trimmed <- trimmed1000.filtered.bed$start.trimmed - 1

# editing attr field to conform to gtf format
trimmed500.filtered$attr  <- paste0("gene_id \"", trimmed500.filtered$attr, "\";")
trimmed750.filtered$attr  <- paste0("gene_id \"", trimmed750.filtered$attr, "\";")
trimmed1000.filtered$attr <- paste0("gene_id \"", trimmed1000.filtered$attr, "\";")

## write out finalized filtered gtfs and bed files
write.table(trimmed500.filtered[,1:9],  "mel.500.trimmed.filtered.gtf",  sep = "\t", quote = F, row.names = F, col.names = F)
write.table(trimmed750.filtered[,1:9],  "mel.750.trimmed.filtered.gtf",  sep = "\t", quote = F, row.names = F, col.names = F)
write.table(trimmed1000.filtered[,1:9], "mel.1000.trimmed.filtered.gtf", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(trimmed500.filtered.bed,  "mel.500.trimmed.filtered.bed",  sep = "\t", quote = F, row.names = F, col.names = F)
write.table(trimmed750.filtered.bed,  "mel.750.trimmed.filtered.bed",  sep = "\t", quote = F, row.names = F, col.names = F)
write.table(trimmed1000.filtered.bed, "mel.1000.trimmed.filtered.bed", sep = "\t", quote = F, row.names = F, col.names = F)