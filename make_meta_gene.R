setwd("../CNV-VLCT/")

################## Catch MIN and MAX coordinates of RefSeq isoforms ##################
# ref_seq <- read.table("refseq_noncoding_genes.txt", stringsAsFactors = F, header = F)
# ref_seq <- read.table("refseq_coding_genes.txt", stringsAsFactors = F, header = F)
ref_seq <- read.table("refseq_allgenes.txt", stringsAsFactors = F, header = F)
colnames(ref_seq) <- c("chr", "start", "end", "gene")

library(data.table)
DT <- data.table(ref_seq)
report_min <- DT[, min(start), by = list(DT$gene,DT$chr)] #list(adShown,url)
report_max <- DT[, max(end), by = list(DT$gene,DT$chr)]
report <- data.frame(cbind(report_min, report_max))

# chr_array = ref_seq[match(report[,1], ref_seq[,4]), 1]
# report <- data.frame(cbind(report, chr_array))
report <- report[,c(2,3,6,1)]
colnames(report) <- c("chr", "start", "end", "gene")

# write.table(x=report, file = "meta_noncoding_gene_refseq.txt", row.names = F, quote = F, sep = "\t")
# write.table(x=report, file = "meta_coding_gene_refseq.txt", row.names = F, quote = F, sep = "\t")
write.table(x=report, file = "meta_allgene_refseq.txt", row.names = F, quote = F, sep = "\t")
rm(report_min, report_max, DT)

#######################################################################################

############### MATCH genes in specific lists #########################################
annotate_gene <- function(bed_file, column_name, gene_list){
  # bedified <- read.table(bed_file, stringsAsFactors = F, header = T, sep = "\t", quote="", fill = T)
  bedified <- read.xlsx(bed_file, sheetIndex = 1, stringsAsFactors = F, header=TRUE, check.names = FALSE, quote="", fill = T)
  bedified[is.na(bedified)] <- "."

  genes_list <- t(read.table(gene_list, stringsAsFactors = F, header = F))
  
  genes_list_in_scope <- rep(".", nrow(bedified))
  for (i in 1:nrow(bedified)){
    genes_in_scope = bedified[i, which(colnames(bedified) == column_name)]
    if (genes_in_scope != "." & genes_in_scope != ""){
      genes_in_scope = unlist(strsplit(genes_in_scope, ","))
      LL <- list()
      for (j in 1:length(genes_in_scope)){
        if (genes_in_scope[j] %in% genes_list) {
          LL[[length(LL)+1]] <- genes_in_scope[j]
        }
      }
      
      if(length(LL) > 0){
        genes_list_in_scope[i] = paste(LL, collapse=',')
      }
    }
  }
  
  return(genes_list_in_scope)
}

count_genes_in_column <- function(column_tobe_counted){
  sapply(column_tobe_counted, function(x) {
    if(x == ".")
      return(0)
    else
      return(length(unlist(strsplit(as.character(x), ","))))
  })
}
 
  require(xlsx)
  annotation_file = "annotated_cnv_file_14022018_500Kb.xlsx" 
  
  epilessia_inside <- annotate_gene(annotation_file, "gene_inside", "epilessia_genelist.txt")
  epilessia_cross  <- annotate_gene(annotation_file, "gene_cross", "epilessia_genelist.txt")
  epilessia_distal  <- annotate_gene(annotation_file, "gene_distal", "epilessia_genelist.txt")
  
  IDa_inside <- annotate_gene(annotation_file, "gene_inside", "ID_a_genelist.txt")
  IDa_cross  <- annotate_gene(annotation_file, "gene_cross", "ID_a_genelist.txt")
  IDa_distal  <- annotate_gene(annotation_file, "gene_distal", "ID_a_genelist.txt")
  
  IDb_inside <- annotate_gene(annotation_file, "gene_inside", "ID_b_genelist.txt")
  IDb_cross  <- annotate_gene(annotation_file, "gene_cross", "ID_b_genelist.txt")
  IDb_distal  <- annotate_gene(annotation_file, "gene_distal", "ID_b_genelist.txt")
  
  malformazioni_inside <- annotate_gene(annotation_file, "gene_inside", "malformazioni_genelist.txt")
  malformazioni_cross  <- annotate_gene(annotation_file, "gene_cross", "malformazioni_genelist.txt")
  malformazioni_distal  <- annotate_gene(annotation_file, "gene_distal", "malformazioni_genelist.txt")
  
  dosage_inside <- annotate_gene(annotation_file, "gene_inside", "dosage_sensitive_genelist.txt")
  dosage_cross  <- annotate_gene(annotation_file, "gene_cross", "dosage_sensitive_genelist.txt")
  dosage_distal  <- annotate_gene(annotation_file, "gene_distal", "dosage_sensitive_genelist.txt")
  
  onologhi_inside <- annotate_gene(annotation_file, "gene_inside", "onologhi_genelist.txt")
  onologhi_cross  <- annotate_gene(annotation_file, "gene_cross", "onologhi_genelist.txt")
  onologhi_distal  <- annotate_gene(annotation_file, "gene_distal", "onologhi_genelist.txt")
  
  mendeliome_inside <- annotate_gene(annotation_file, "gene_inside", "mendeliome_genelist.txt")
  mendeliome_cross  <- annotate_gene(annotation_file, "gene_cross", "mendeliome_genelist.txt")
  mendeliome_distal  <- annotate_gene(annotation_file, "gene_distal", "mendeliome_genelist.txt")
  
  SFARI_noscore_inside <- annotate_gene(annotation_file, "gene_inside", "gene_noscore_genelist.txt")
  SFARI_noscore_cross <- annotate_gene(annotation_file, "gene_cross", "gene_noscore_genelist.txt")
  SFARI_noscore_distal <- annotate_gene(annotation_file, "gene_distal", "gene_noscore_genelist.txt")
  
  SFARI_1score_inside <- annotate_gene(annotation_file, "gene_inside", "gene-score_1_genelist.txt")
  SFARI_1score_cross <- annotate_gene(annotation_file, "gene_cross", "gene-score_1_genelist.txt")
  SFARI_1score_distal <- annotate_gene(annotation_file, "gene_distal", "gene-score_1_genelist.txt")
  
  SFARI_2score_inside <- annotate_gene(annotation_file, "gene_inside", "gene-score_2_genelist.txt")
  SFARI_2score_cross <- annotate_gene(annotation_file, "gene_cross", "gene-score_2_genelist.txt")
  SFARI_2score_distal <- annotate_gene(annotation_file, "gene_distal", "gene-score_2_genelist.txt")
  
  SFARI_3score_inside <- annotate_gene(annotation_file, "gene_inside", "gene-score_3_genelist.txt")
  SFARI_3score_cross <- annotate_gene(annotation_file, "gene_cross", "gene-score_3_genelist.txt")
  SFARI_3score_distal <- annotate_gene(annotation_file, "gene_distal", "gene-score_3_genelist.txt")
  
  SFARI_4score_inside <- annotate_gene(annotation_file, "gene_inside", "gene-score_4_genelist.txt")
  SFARI_4score_cross <- annotate_gene(annotation_file, "gene_cross", "gene-score_4_genelist.txt")
  SFARI_4score_distal <- annotate_gene(annotation_file, "gene_distal", "gene-score_4_genelist.txt")
  
  SFARI_5score_inside <- annotate_gene(annotation_file, "gene_inside", "gene-score_5_genelist.txt")
  SFARI_5score_cross <- annotate_gene(annotation_file, "gene_cross", "gene-score_5_genelist.txt")
  SFARI_5score_distal <- annotate_gene(annotation_file, "gene_distal", "gene-score_5_genelist.txt")
  
  SFARI_6score_inside <- annotate_gene(annotation_file, "gene_inside", "gene-score_6_genelist.txt")
  SFARI_6score_cross <- annotate_gene(annotation_file, "gene_cross", "gene-score_6_genelist.txt")
  SFARI_6score_distal <- annotate_gene(annotation_file, "gene_distal", "gene-score_6_genelist.txt")
  
  SFARI_Sscore_inside <- annotate_gene(annotation_file, "gene_inside", "gene-score_S_genelist.txt")
  SFARI_Sscore_cross <- annotate_gene(annotation_file, "gene_cross", "gene-score_S_genelist.txt")
  SFARI_Sscore_distal <- annotate_gene(annotation_file, "gene_distal", "gene-score_S_genelist.txt")
  
  ASD_inside <- annotate_gene(annotation_file, "gene_inside", "ASD_genelist.txt")
  ASD_cross <- annotate_gene(annotation_file, "gene_cross", "ASD_genelist.txt")
  ASD_distal <- annotate_gene(annotation_file, "gene_distal", "ASD_genelist.txt")
  
  ### Laura specific lists
  autism_pubmed_inside <- annotate_gene(annotation_file, "gene_inside", "pubmed_09-02-2018_autism_genelist.txt")
  autism_pubmed_cross <- annotate_gene(annotation_file, "gene_cross", "pubmed_09-02-2018_autism_genelist.txt")
  autism_pubmed_distal <- annotate_gene(annotation_file, "gene_distal", "pubmed_09-02-2018_autism_genelist.txt")
  
  brain_malformation_pubmed_inside <- annotate_gene(annotation_file, "gene_inside", "pubmed_09-02-2018_brain malformations_genelist.txt")
  brain_malformation_pubmed_cross <- annotate_gene(annotation_file, "gene_cross", "pubmed_09-02-2018_brain malformations_genelist.txt")
  brain_malformation_pubmed_distal <- annotate_gene(annotation_file, "gene_distal", "pubmed_09-02-2018_brain malformations_genelist.txt")
  
  epilepsy_pubmed_inside <- annotate_gene(annotation_file, "gene_inside", "pubmed_09-02-2018_epilepsy or seizures_genelist.txt")
  epilepsy_pubmed_cross <- annotate_gene(annotation_file, "gene_cross", "pubmed_09-02-2018_epilepsy or seizures_genelist.txt")
  epilepsy_pubmed_distal <- annotate_gene(annotation_file, "gene_distal", "pubmed_09-02-2018_epilepsy or seizures_genelist.txt")
  
  ID_pubmed_inside <- annotate_gene(annotation_file, "gene_inside", "pubmed_09-02-2018_intellectual disability_genelist.txt")
  ID_pubmed_cross <- annotate_gene(annotation_file, "gene_cross", "pubmed_09-02-2018_intellectual disability_genelist.txt")
  ID_pubmed_distal <- annotate_gene(annotation_file, "gene_distal", "pubmed_09-02-2018_intellectual disability_genelist.txt")
  #########################
  
  # pseudogene_inside <- annotate_gene(annotation_file, "gene_inside", "pseudogenes_genelist.txt")
  # pseudogene_cross <- annotate_gene(annotation_file, "gene_cross", "pseudogenes_genelist.txt")
  # pseudogene_distal <- annotate_gene(annotation_file, "gene_distal", "pseudogenes_genelist.txt")
  
  ## CBIND these columns to the original file
  
  require(xlsx)
  f <- read.xlsx(annotation_file, sheetIndex = 1, stringsAsFactors = F, header=TRUE, check.names = FALSE, quote="", fill = T)
  # f <- read.table(annotation_file, stringsAsFactors = F, header = T, sep = "\t", quote="", fill=T, check.names=FALSE)
  f[is.na(f)] <- "."
  
  ######### EPILESSIA #############################################
  f = cbind(f, epilessia_inside)
  epilessia_inside_num <- count_genes_in_column(f$epilessia_inside)
  f = cbind(f, epilessia_inside_num)
  
  f = cbind(f, epilessia_cross)
  epilessia_cross_num <- count_genes_in_column(f$epilessia_cross)
  f = cbind(f, epilessia_cross_num)
  
  f = cbind(f, epilessia_distal)
  epilessia_distal_num <- count_genes_in_column(f$epilessia_distal)
  f = cbind(f, epilessia_distal_num)
  #################################################################
  
  ############### IDa #############################################
  f = cbind(f, IDa_inside)
  IDa_inside_num <- count_genes_in_column(f$IDa_inside)
  f = cbind(f, IDa_inside_num)
  
  f = cbind(f, IDa_cross)
  IDa_cross_num <- count_genes_in_column(f$IDa_cross)
  f = cbind(f, IDa_cross_num)
  
  f = cbind(f, IDa_distal)
  IDa_distal_num <- count_genes_in_column(f$IDa_distal)
  f = cbind(f, IDa_distal_num)
  #################################################################
  
  ############### IDb #############################################
  f = cbind(f, IDb_inside)
  IDb_inside_num <- count_genes_in_column(f$IDb_inside)
  f = cbind(f, IDb_inside_num)
  
  f = cbind(f, IDb_cross)
  IDb_cross_num <- count_genes_in_column(f$IDb_cross)
  f = cbind(f, IDb_cross_num)
  
  f = cbind(f, IDb_distal)
  IDb_distal_num <- count_genes_in_column(f$IDb_distal)
  f = cbind(f, IDb_distal_num)
  #################################################################
  
  ############### malformazioni ###################################
  f = cbind(f, malformazioni_inside)
  malformazioni_inside_num <- count_genes_in_column(f$malformazioni_inside)
  f = cbind(f, malformazioni_inside_num)
  
  f = cbind(f, malformazioni_cross)
  malformazioni_cross_num <- count_genes_in_column(f$malformazioni_cross)
  f = cbind(f, malformazioni_cross_num)
  
  f = cbind(f, malformazioni_distal)
  malformazioni_distal_num <- count_genes_in_column(f$malformazioni_distal)
  f = cbind(f, malformazioni_distal_num)
  #################################################################
  
  ###################### dosage ###################################
  f = cbind(f, dosage_inside)
  dosage_inside_num <- count_genes_in_column(f$dosage_inside)
  f = cbind(f, dosage_inside_num)
  
  f = cbind(f, dosage_cross)
  dosage_cross_num <- count_genes_in_column(f$dosage_cross)
  f = cbind(f, dosage_cross_num)
  
  f = cbind(f, dosage_distal)
  dosage_distal_num <- count_genes_in_column(f$dosage_distal)
  f = cbind(f, dosage_distal_num)
  #################################################################
  
  ########################## onologhi #############################
  f = cbind(f, onologhi_inside)
  onologhi_inside_num <- count_genes_in_column(f$onologhi_inside)
  f = cbind(f, onologhi_inside_num)
  
  f = cbind(f, onologhi_cross)
  onologhi_cross_num <- count_genes_in_column(f$onologhi_cross)
  f = cbind(f, onologhi_cross_num)
  
  f = cbind(f, onologhi_distal)
  onologhi_distal_num <- count_genes_in_column(f$onologhi_distal)
  f = cbind(f, onologhi_distal_num)
  #################################################################
  
  ########################## mendeliome ###########################
  f = cbind(f, mendeliome_inside)
  mendeliome_inside_num <- count_genes_in_column(f$mendeliome_inside)
  f = cbind(f, mendeliome_inside_num)
  
  f = cbind(f, mendeliome_cross)
  mendeliome_cross_num <- count_genes_in_column(f$mendeliome_cross)
  f = cbind(f, mendeliome_cross_num)
  
  f = cbind(f, mendeliome_distal)
  mendeliome_distal_num <- count_genes_in_column(f$mendeliome_distal)
  f = cbind(f, mendeliome_distal_num)
  #################################################################
  
  ############################ SFARI ##############################
  f = cbind(f, SFARI_noscore_inside)
  SFARI_noscore_inside_num <- count_genes_in_column(f$SFARI_noscore_inside)
  f = cbind(f, SFARI_noscore_inside_num)
  
  f = cbind(f, SFARI_noscore_cross)
  SFARI_noscore_cross_num <- count_genes_in_column(f$SFARI_noscore_cross)
  f = cbind(f, SFARI_noscore_cross_num)
  
  f = cbind(f, SFARI_noscore_distal)
  SFARI_noscore_distal_num <- count_genes_in_column(f$SFARI_noscore_distal)
  f = cbind(f, SFARI_noscore_distal_num)
  
  ##
  
  f = cbind(f, SFARI_1score_inside)
  SFARI_1score_inside_num <- count_genes_in_column(f$SFARI_1score_inside)
  f = cbind(f, SFARI_1score_inside_num)
  
  f = cbind(f, SFARI_1score_cross)
  SFARI_1score_cross_num <- count_genes_in_column(f$SFARI_1score_cross)
  f = cbind(f, SFARI_1score_cross_num)
  
  f = cbind(f, SFARI_1score_distal)
  SFARI_1score_distal_num <- count_genes_in_column(f$SFARI_1score_distal)
  f = cbind(f, SFARI_1score_distal_num)
  
  ##
  
  f = cbind(f, SFARI_2score_inside)
  SFARI_2score_inside_num <- count_genes_in_column(f$SFARI_2score_inside)
  f = cbind(f, SFARI_2score_inside_num)
  
  f = cbind(f, SFARI_2score_cross)
  SFARI_2score_cross_num <- count_genes_in_column(f$SFARI_2score_cross)
  f = cbind(f, SFARI_2score_cross_num)
  
  f = cbind(f, SFARI_2score_distal)
  SFARI_2score_distal_num <- count_genes_in_column(f$SFARI_2score_distal)
  f = cbind(f, SFARI_2score_distal_num)
  
  ##
  
  f = cbind(f, SFARI_3score_inside)
  SFARI_3score_inside_num <- count_genes_in_column(f$SFARI_3score_inside)
  f = cbind(f, SFARI_3score_inside_num)
  
  f = cbind(f, SFARI_3score_cross)
  SFARI_3score_cross_num <- count_genes_in_column(f$SFARI_3score_cross)
  f = cbind(f, SFARI_3score_cross_num)
  
  f = cbind(f, SFARI_3score_distal)
  SFARI_3score_distal_num <- count_genes_in_column(f$SFARI_3score_distal)
  f = cbind(f, SFARI_3score_distal_num)
  
  ##
  
  f = cbind(f, SFARI_4score_inside)
  SFARI_4score_inside_num <- count_genes_in_column(f$SFARI_4score_inside)
  f = cbind(f, SFARI_4score_inside_num)
  
  f = cbind(f, SFARI_4score_cross)
  SFARI_4score_cross_num <- count_genes_in_column(f$SFARI_4score_cross)
  f = cbind(f, SFARI_4score_cross_num)
  
  f = cbind(f, SFARI_4score_distal)
  SFARI_4score_distal_num <- count_genes_in_column(f$SFARI_4score_distal)
  f = cbind(f, SFARI_4score_distal_num)
  
  ##
  
  f = cbind(f, SFARI_5score_inside)
  SFARI_5score_inside_num <- count_genes_in_column(f$SFARI_5score_inside)
  f = cbind(f, SFARI_5score_inside_num)
  
  f = cbind(f, SFARI_5score_cross)
  SFARI_5score_cross_num <- count_genes_in_column(f$SFARI_5score_cross)
  f = cbind(f, SFARI_5score_cross_num)
  
  f = cbind(f, SFARI_5score_distal)
  SFARI_5score_distal_num <- count_genes_in_column(f$SFARI_5score_distal)
  f = cbind(f, SFARI_5score_distal_num)
  
  ##
  
  f = cbind(f, SFARI_6score_inside)
  SFARI_6score_inside_num <- count_genes_in_column(f$SFARI_6score_inside)
  f = cbind(f, SFARI_6score_inside_num)
  
  f = cbind(f, SFARI_6score_cross)
  SFARI_6score_cross_num <- count_genes_in_column(f$SFARI_6score_cross)
  f = cbind(f, SFARI_6score_cross_num)
  
  f = cbind(f, SFARI_6score_distal)
  SFARI_6score_distal_num <- count_genes_in_column(f$SFARI_6score_distal)
  f = cbind(f, SFARI_6score_distal_num)
  
  ##
  
  f = cbind(f, SFARI_Sscore_inside)
  SFARI_Sscore_inside_num <- count_genes_in_column(f$SFARI_Sscore_inside)
  f = cbind(f, SFARI_Sscore_inside_num)
  
  f = cbind(f, SFARI_Sscore_cross)
  SFARI_Sscore_cross_num <- count_genes_in_column(f$SFARI_Sscore_cross)
  f = cbind(f, SFARI_Sscore_cross_num)
  
  f = cbind(f, SFARI_Sscore_distal)
  SFARI_Sscore_distal_num <- count_genes_in_column(f$SFARI_Sscore_distal)
  f = cbind(f, SFARI_Sscore_distal_num)
  #################################################################
  
  ############################# ASD ###############################
  f = cbind(f, ASD_inside)
  ASD_inside_num <- count_genes_in_column(f$ASD_inside)
  f = cbind(f, ASD_inside_num)
  
  f = cbind(f, ASD_cross)
  ASD_cross_num <- count_genes_in_column(f$ASD_cross)
  f = cbind(f, ASD_cross_num)
  
  f = cbind(f, ASD_distal)
  ASD_distal_num <- count_genes_in_column(f$ASD_distal)
  f = cbind(f, ASD_distal_num)
  #################################################################
  
  ### Laura's specific lists
  f = cbind(f, autism_pubmed_inside)
  autism_pubmed_inside_num <- count_genes_in_column(f$autism_pubmed_inside)
  f = cbind(f, autism_pubmed_inside_num)
  
  f = cbind(f, autism_pubmed_cross)
  autism_pubmed_cross_num <- count_genes_in_column(f$autism_pubmed_cross)
  f = cbind(f, autism_pubmed_cross_num)
  
  f = cbind(f, autism_pubmed_distal)
  autism_pubmed_distal_num <- count_genes_in_column(f$autism_pubmed_distal)
  f = cbind(f, autism_pubmed_distal_num)
  
  #
  
  f = cbind(f, brain_malformation_pubmed_inside)
  brain_malformation_pubmed_inside_num <- count_genes_in_column(f$brain_malformation_pubmed_inside)
  f = cbind(f, brain_malformation_pubmed_inside_num)
  
  f = cbind(f, brain_malformation_pubmed_cross)
  brain_malformation_pubmed_cross_num <- count_genes_in_column(f$brain_malformation_pubmed_cross)
  f = cbind(f, brain_malformation_pubmed_cross_num)
  
  f = cbind(f, brain_malformation_pubmed_distal)
  brain_malformation_pubmed_distal_num <- count_genes_in_column(f$brain_malformation_pubmed_distal)
  f = cbind(f, brain_malformation_pubmed_distal_num)
  
  #
  
  f = cbind(f, epilepsy_pubmed_inside)
  epilepsy_pubmed_inside_num <- count_genes_in_column(f$epilepsy_pubmed_inside)
  f = cbind(f, epilepsy_pubmed_inside_num)
  
  f = cbind(f, epilepsy_pubmed_cross)
  epilepsy_pubmed_cross_num <- count_genes_in_column(f$epilepsy_pubmed_cross)
  f = cbind(f, epilepsy_pubmed_cross_num)
  
  f = cbind(f, epilepsy_pubmed_distal)
  epilepsy_pubmed_distal_num <- count_genes_in_column(f$epilepsy_pubmed_distal)
  f = cbind(f, epilepsy_pubmed_distal_num)
  
  #
  
  f = cbind(f, ID_pubmed_inside)
  ID_pubmed_inside_num <- count_genes_in_column(f$ID_pubmed_inside)
  f = cbind(f, ID_pubmed_inside_num)
  
  f = cbind(f, ID_pubmed_cross)
  ID_pubmed_cross_num <- count_genes_in_column(f$ID_pubmed_cross)
  f = cbind(f, ID_pubmed_cross_num)
  
  f = cbind(f, ID_pubmed_distal)
  ID_pubmed_distal_num <- count_genes_in_column(f$ID_pubmed_distal)
  f = cbind(f, ID_pubmed_distal_num)
  #########################
  
  
  ############################# PSEUDOGENES########################
  # f = cbind(f, pseudogene_inside)
  # pseudogene_inside_num <- count_genes_in_column(f$pseudogene_inside)
  # f = cbind(f, pseudogene_inside_num)
  # 
  # f = cbind(f, pseudogene_cross)
  # pseudogene_cross_num <- count_genes_in_column(f$pseudogene_cross)
  # f = cbind(f, pseudogene_cross_num)
  # 
  # f = cbind(f, pseudogene_distal)
  # pseudogene_distal_num <- count_genes_in_column(f$pseudogene_distal)
  # f = cbind(f, pseudogene_distal_num)
  #################################################################
  
  
  
  library(xlsx)
  write.xlsx(x = f, file = "annotated_cnv_file_14022018_500Kb.xlsx", sheetName = "Annotated_CNV", row.names = FALSE)
#}

##########################################################

