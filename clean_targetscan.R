setwd("V:/CNV-VLCT/AnnotateCNV")

ts <- read.table("TargetScan_71_Conserved_Site_Context_Scores.txt", stringsAsFactors = F, header = T, fill = TRUE)
ts <- ts[,c(2,5)]
colnames(ts) <- c("Gene Symbol", "miRNA")

ts.hsa <- unique(ts[grep("hsa-", ts$miRNA),])

# library(xlsx)
# write.xlsx(x = ts.hsa, file = "TargetScan_71_Conserved_Site_filtered.xlsx", sheetName = "TargetScan71_filtered", row.names = FALSE)
write.table(x = ts.hsa, file = "TargetScan_71_Conserved_Site_filtered.txt", sep="\t", row.names = F)


##########################################################

