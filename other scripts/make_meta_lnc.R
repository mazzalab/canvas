setwd("../resources/BED")

################## Catch MIN and MAX coordinates of RefSeq isoforms ##################
lnc <- read.table("gene.txt", stringsAsFactors = F, header = F)
colnames(lnc) <- c("chr", "start", "end", "lnc")

lnc <- cbind(lnc, gsub(":.*$", "", lnc$lnc))
lnc <- lnc[,-4]
colnames(lnc) <- c("chr", "start", "end", "lnc")

library(data.table)
DT <- data.table(lnc)
report_min <- DT[, min(start), by = list(DT$lnc,DT$chr)]
report_max <- DT[, max(end), by = list(DT$lnc,DT$chr)]
report <- data.frame(cbind(report_min, report_max))
report <- report[,c(2,3,6,1)]
colnames(report) <- c("chr", "start", "end", "lnc")
rm(report_min, report_max, DT)

report <- report[order(report$lnc),]
write.table(x=report, file = "meta_gene.txt", row.names = F, quote = F, sep = "\t")


#######################################################################################

  
  
library(xlsx)
write.xlsx(x = report, file = "meta_gene.txt", sheetName = "LongNC", row.names = FALSE)

##########################################################

