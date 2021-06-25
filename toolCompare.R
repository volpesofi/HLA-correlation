# ----- import library ---------
suppressMessages(library(dplyr))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(stringr))
suppressMessages(library("gplots"))
suppressMessages(library(tidyr))
suppressMessages(library(viridis))
library(plyr)
library(scales)

# ----- set WorkingDirectory and create outputs folders ---------
setwd('/beegfs/scratch/ric.cosr/ric.cosr/Vo_WGS/hla_statistics/')

# ------ function --------
MAtrix_prepare = function(file, sep = "\t", digit = 2, rownames = "Barcode") {
  HLA_t = read.table(file = file , header = T, sep = sep, check.names = F)
  HLA_t$Barcode = as.character(HLA_t$Barcode)
  HLA_allele= HLA_t[,grep(x = colnames(HLA_t),pattern = "HLA")]
  fx <- function(x, N) {
    a = strsplit(x, '\\*')
    for (i in 1:length(a)) { 
      digit = unlist(str_split(string = a[[i]][2], pattern = ":"))
      a[[i]][2] = paste(digit[1:N], collapse = ":")
      if (a[[i]][1] == "n/a") {
        a[[i]] = a[[i]][1]
      } else {
        a[[i]] = paste(a[[i]], collapse = "*")
      }
    }
    return(unlist(a))
  }
  HLA_twodigit = as.data.frame(sapply(HLA_allele, fx, N = digit))
  row.names(HLA_twodigit) <- HLA_t[,rownames]
  colnames(HLA_twodigit) <- str_replace_all(string = colnames(HLA_twodigit), 
                                            pattern = "-", replacement = "")
  colnames(HLA_twodigit) <- str_replace_all(string = colnames(HLA_twodigit), 
                                            pattern = " ", replacement = "_")
  return(HLA_twodigit)
}
#file1 = file_hla_typing_1; file2= file_hla_typing_2; digit = 2; rownames = "Barcode";ALLELE_list = c("A", "B", "C")
toolCOMPARE = function(file1, file2, digit = 2, rownames = "Barcode", ALLELE_list = c("A", "B", "C")) {
  M1 = MAtrix_prepare(file = file1, digit =2, rownames = rownames)
  M2 = MAtrix_prepare(file = file2, digit =2, rownames = rownames)
  patients = as.character(unique(c(row.names(M1), row.names(M2))))
  toolComparison = list()
  for (allele in ALLELE_list) {
    df = data.frame()
    for (p in patients) {
      index1 = grep(pattern = paste("HLA",allele,sep=''), x = colnames(M1))
      index2 = grep(pattern = paste("HLA",allele,sep=''), x = colnames(M1))
      T1 = str_remove_all(as.character(M1[p,index1]), pattern = "HLA-"); 
      T2 = str_remove_all(as.character(M2[p,index2]), pattern = "HLA-"); 
      hom = length(intersect(T1,T2))
      df = rbind(df, data.frame(ID = p, homology = hom))
    }
    # save results
    toolComparison[[allele]] = df
    write.table(df, paste(tabl, "HLA-", allele, "comparison.tsv", sep=''), sep= "\t", row.names = F, quote = F)
    ## plot
    tb <- data.frame(table(df$homology))
    colnames(tb) <- c('homology','Freq')
    tb$toolCorrespondence = plyr::revalue(as.factor(tb$homology), c("2"="perfect", "1"="partial", "0"="none"))
    tb$Perc <- tb$Freq / sum(tb$Freq)*100
    write.table(tb, paste(tabl, "HLA-", allele, "perc.tsv", sep=''), sep= "\t", row.names = F, quote = F)
    bp<- ggplot(tb, aes(x=" ", y = Freq, fill=toolCorrespondence))+
      geom_bar(width = 1, stat = "identity", color = "black") + 
      scale_fill_brewer(palette="Dark2") +
      ggtitle(paste("HLA-", allele, "_", str_sub(file1, 7, -5), "_VS_" , str_sub(file2, 7, -5), sep ="")) +
      theme_minimal() +
      xlab(" ") + ylab("") +
      geom_text(aes(label = paste(round(Freq / sum(Freq) * 100, 1), "%")),
                position = position_stack(vjust = 0.5), 
                color = "white", family = "Times New Roman", 
                fontface = "bold") +
      coord_polar("y", start=0)
    ggsave(paste(plt, "HLA-", allele, "pie.png", sep=''),bp, width = 7, height = 5)
  } 
  return(toolComparison)
}

# ------ analysis --------
file_hla_typing_1 = 'input/hlavbseq_results.tsv'
file_hla_typing_2 = 'input/fished_all_tops.tsv'
file_hla_typing_1 = 'input/optitype_results_classI.tsv'
file_hla_typing_2 = 'input/fished_all_tops.tsv'

tabl = paste("tables_", str_sub(file_hla_typing_1, 7, -5), "_VS_" , str_sub(file_hla_typing_2, 7, -5), "/",  sep ='') 
plt = paste("plots_", str_sub(file_hla_typing_1, 7, -5), "_VS_" , str_sub(file_hla_typing_2, 7, -5), "/",  sep ='') 
dir.create(tabl, recursive = T, showWarnings = F)
dir.create(plt, recursive = T, showWarnings = F)

toolsCL_list = toolCOMPARE(file1 = file_hla_typing_1, 
                      file2 = file_hla_typing_2, 
                      digit = 2, 
                      rownames = "Barcode", 
                      ALLELE_list = c("A", "B", "C", "DQA1", "DQB1", "DRB1"))


toolsCL_list = toolCOMPARE(file1 = file_hla_typing_1, 
                           file2 = file_hla_typing_2, 
                           digit = 2, 
                           rownames = "Barcode", 
                           ALLELE_list = c("A", "B", "C"))

toolsCL_df = Reduce(function(x, y) merge(x, y, by = "ID"), toolsCL_list)

colnames(toolsCL_df) = c("ID","A", "B", "C", "DQA1", "DQB1", "DRB1")
colnames(toolsCL_df) = c("ID","A", "B", "C") #, "DQA1", "DQB1", "DRB1")

write.table(toolsCL_df, paste(tabl, "Summary_comparison.tsv", sep=''), sep= "\t", row.names = F, quote = F)

