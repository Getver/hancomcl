library(dplyr)
library(tidyr)



R <- read.table('/disk0/sm/microbiome/16S/1013_park/Run/RESULT_ASV_TABLE.txt', sep='\t', header=TRUE)


R0 <- read.table('/disk0/sm/microbiome/16S/1013_park/Run0/RESULT_ASV_TABLE.txt', sep='\t', header=TRUE)
R7 <- read.table('/disk0/sm/microbiome/16S/1013_park/Run7/RESULT_ASV_TABLE.txt', sep='\t', header=TRUE)

R0_col <- names(R0)
R7_col <- names(R7)

R0 <- R %>% select(R0_col)
R7 <- R %>% select(R7_col)

B3 <- read.table('/disk0/sm/microbiome/16S/1013_park/RESULT.txt', sep='\t', header=TRUE)





####
target_sample0 = 'S_C01_S1_L001'    # S_C01_S1_L001 S_S02_S12_L001 S_W03_S23_L001
target_sample7 = 'S_C71_S6_L001'    # S_C71_S6_L001 S_S74_S18_L001 S_W76_S29_L001
D0_sample = 'C01'  # C01    S12 W03
D7_sample = 'C71'  # C71    S18 W76



denom01 <- sum(B3[[D0_sample]], na.rm = TRUE)
denom71 <- sum(B3[[D7_sample]], na.rm = TRUE)

target0 <- R0 %>% filter(Identity >= 99) %>% select(Species, all_of(target_sample0))
target7 <- R7 %>% filter(Identity >= 99) %>% select(Species, all_of(target_sample7))

target3 <- B3 %>% filter(Identity >= 99) %>% mutate(S=paste(Genus, Species), D0 = .data[[D0_sample]] / denom01, D7 = .data[[D7_sample]] / denom71) %>% select(S, D0, D7)


result0 <- target0 %>% group_by(Species) %>% summarise(D0 = sum(.data[[target_sample0]], na.rm = TRUE)) %>% arrange(desc(D0))
result7 <- target7 %>% group_by(Species) %>% summarise(D7 = sum(.data[[target_sample7]], na.rm = TRUE)) %>% arrange(desc(D7))

result3 <- target3 %>% group_by(S) %>% summarise(B0 = sum(D0, na.rm = TRUE), B7 = sum(D7, na.rm = TRUE), .groups = "drop")


com0 <- merge(result0, result3, by.x = "Species", by.y = "S", all = TRUE)
com7 <- merge(result7, result3, by.x = "Species", by.y = "S", all = TRUE)

compare <- merge(com0, com7, by = c("Species", 'B0', 'B7'), all = TRUE) %>% mutate(across(where(is.numeric), ~replace_na(.x, 0)))


write.table(compare, '/disk0/sm/microbiome/16S/1013_park/RESULT_compare.txt', sep='\t', quote=FALSE, row.names=FALSE)




# target0 <- R0 %>% filter(Identity >= 99) %>% select(c(Species, all_of(target_sample0)))
# target7 <- R7 %>% filter(Identity >= 99) %>% select(c(Species, all_of(target_sample7)))


# target3 <- B3 %>% mutate(S=paste0(Genus, ' ', Species), c01=D0_sample/sum(D0_sample), c71=D7_sample/sum(D7_sample)) %>% filter(Identity >= 99) %>% select(c(S, c01, c71))



# result0 <- target0 %>% group_by(Species) %>% summarise(D0_C01 = sum(S_C01_S1_L001, na.rm = TRUE)) %>% arrange(desc(D0_C01))
# result7 <- target7 %>% group_by(Species) %>% summarise(D7_C71 = sum(S_C71_S6_L001, na.rm = TRUE)) %>% arrange(desc(D7_C71))

# result3 <- target3 %>% group_by(S) %>% summarise(abundance01 = sum(c01, na.rm = TRUE), abundance71 = sum(c71, na.rm = TRUE)) %>% arrange(desc(abundance01))




# com0 <- merge(result0, result3, by.x='Species', by.y='S', all=TRUE)
# com7 <- merge(result7, result3, by.x='Species', by.y='S', all=TRUE)



# compare <- merge(com0, com7, by=c('Species', 'abundance71', 'abundance01'), all='TRUE') # %>% replace_na(list(abundance01 = 0, abundance71 = 0))
# compare <- compare %>% mutate(across(where(is.numeric), ~replace_na(.x, 0)))

# write.table(compare, '/disk0/sm/microbiome/16S/1013_park/RESULT_compare.txt', sep='\t', quote=FALSE, row.names=FALSE)




bad_target0 <- R0 %>% filter(Identity < 99) %>% select(c(Species, S_C01_S1_L001))
bad_result0 <- bad_target0 %>% group_by(Species) %>% summarise(D0_C01 = sum(S_C01_S1_L001, na.rm = TRUE)) %>% arrange(desc(D0_C01))

