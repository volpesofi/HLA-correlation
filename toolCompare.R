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
library("patchwork")

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
    tb$toolCorrespondence = plyr::revalue(as.factor(tb$homology), c("2"="2 alleles", "1"="1 alleles", "0"="none"))
    tb$Perc <- tb$Freq / sum(tb$Freq)*100
    write.table(tb, paste(tabl, "HLA-", allele, "perc.tsv", sep=''), sep= "\t", row.names = F, quote = F)
    write.table(tb[tb$toolCorrespondence == "2 alleles", "ID"], 
                paste(tabl, "HLA-", allele, "ID_2alleles.tsv", sep=''), sep= "\t", row.names = F, quote = F)
    bp<- ggplot(tb, aes(x=" ", y = Freq, fill=toolCorrespondence))+
      geom_bar(width = 1, stat = "identity", color = "black") + 
      #scale_fill_brewer(palette="YlGnBu",  direction=-1, type = "seq") +
      scale_fill_manual(values = c("#FF0000","#FFCC33", "#339900")) +
      #scale_fill_viridis(end = 0.8) +
      ggtitle(paste("HLA-", allele, "_", str_sub(file1, 7, -5), "_VS_" , str_sub(file2, 7, -5), sep ="")) +
      theme_minimal() +
      xlab(" ") + ylab("") +
      geom_text(aes(label = paste(round(Freq / sum(Freq) * 100, 1), "%")),
                position = position_stack(vjust = 0.5), 
                color = "white", 
                #family = "Times New Roman", 
                fontface = "bold") +
      coord_polar("y", start=0)
    ggsave(paste(plt, "HLA-", allele, "pie.png", sep=''),bp, width = 7, height = 5)
    ggsave(paste(plt, "HLA-", allele, "pie.pdf", sep=''),bp, device = "pdf",width = 7, height = 5)
  } 
  return(toolComparison)
}

# ------ analysis --------
file_hla_typing_1 = 'input/hlavbseq_results.tsv'
file_hla_typing_2 = 'input/fished_all_tops.tsv'
#file_hla_typing_1 = 'input/optitype_results_classI.tsv'
#file_hla_typing_2 = 'input/optitype_results_classI.tsv'

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


### -------- statistics on reliability HLA genes ----------
suppressMessages(library(venn))
file1 = file_hla_typing_1; file2 = file_hla_typing_2
rownames = "Barcode"; allele = "A"

toolCompare_HLA = function(file1, file2, rownames = "Barcode", allele = "A") {
  M1 = MAtrix_prepare(file = file1, digit =2, rownames = rownames)
  M2 = MAtrix_prepare(file = file2, digit =2, rownames = rownames)
  index1 = grep(pattern = paste("HLA",allele, sep=''), x = colnames(M1))
  index2 = grep(pattern = paste("HLA",allele, sep=''), x = colnames(M2))
  
  HLA1 =unique(c(M1[,index1[1]], M1[,index1[2]]))
  HLA2 = unique(c(M2[,index1[1]], M2[,index1[2]]))
  a = venn::venn(list(tool1 = HLA1, tool2 = HLA2), 
                 zcolor = c("#0000FF", "#FF3300"), box = FALSE, 
                 ggplot = FALSE)
  ggsave(filename = paste(plt, "HLA-", allele,"_common.png", sep = ""))
  listInt = attributes(a)$intersections
  write.xlsx(listInt, file = paste(tabl, "HLA-", allele,"_common.xlsx", sep = ""))
  
  hla.unique = sort(unique(c(HLA1, HLA2)))
  samples = unique(c(row.names(M1), row.names(M2)))
  # initialize matrices
  M_tool1 = matrix(0, nrow = length(samples), ncol = length(hla.unique), 
                   dimnames = list(sampleID = samples, hla = hla.unique))
  M_tool2 = matrix(0, nrow = length(samples), ncol = length(hla.unique), 
                   dimnames = list(sampleID = samples, hla = hla.unique))
  # fill matrices
  for (s in samples) {
    if (s %in% row.names(M1)) {
      # hla of tool 1 for patient s
      hla.1 = M1[s,index[1]] 
      hla.2 = M1[s,index[2]]
      # add count to matrix tool 1
      M_tool1[s,hla.1] = M_tool1[s,hla.1]+1
      M_tool1[s,hla.2] = M_tool1[s,hla.2]+1
    } 
    if (s %in% row.names(M2)) {
      # hla of tool 2 for patient s
      hla.1 = M2[s,index[1]] 
      hla.2 = M2[s,index[2]]
      # add count to matrix tool 2
      M_tool2[s,hla.1] = M_tool2[s,hla.1]+1
      M_tool2[s,hla.2] = M_tool2[s,hla.2]+1
    }
  }
  
  # diff matrix
  M_diff = M_tool1 - M_tool2
  myb = seq(0,2,by = 0.01)
  myc = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(myb))
  
  pheatmap(abs(M_diff), 
           border_color = 'darkgrey',
           color = myc, 
           breaks = myb,
           cluster_cols = T, cluster_rows = T,
           cellwidth = 6, cellheight = 1, 
           fontsize = 5, fontsize_row = 1, fontsize_col = 5)
  scores = colSums(abs(M_diff))/ (colSums(M_tool1)+colSums(M_tool2))
  HLA_error = data.frame(name = names(scores), 
                         error_Rate = scores, 
                         countT1 = colSums(M_tool1),
                         countT2 = colSums(M_tool2))
  HLA_error_s = HLA_error[order(HLA_error$error_Rate, decreasing = T),]
  #pheatmap(HLA_error$error_Rate)
  HLA_error_s$name = as.factor(HLA_error_s$name)
  HLA_error_s$name = factor(HLA_error_s$name, levels = HLA_error_s$name)
  
  ER = ggplot(data = HLA_error_s, aes(x="", y=name, fill = error_Rate)) + 
    #geom_raster() + 
    geom_tile(color = "darkgrey") +
    scale_fill_gradientn(colours=rev(brewer.pal(n = 7, name ="PiYG")),
                         na.value = "transparent",
                         breaks=c(0,0.25,0.5,0.75,1),
                         labels=c(0,0.25,0.5,0.75,1),
                         limits=c(0,1)) + 
    theme_minimal()
  CT1 = ggplot(data = HLA_error_s, aes(x="", y=name, fill=countT1)) + 
    geom_tile(color = "darkgrey") + 
    scale_fill_gradientn(colours=rev(brewer.pal(n = 7, name ="RdBu")),
                         na.value = "darkred",
                         breaks=c(0,5,10,15,20),
                         labels=c(0,5,10,15,20),
                         limits=c(0,20)) + theme_minimal()
  
  CT2 = ggplot(data = HLA_error_s, aes(x="", y=name, fill=countT2)) + 
    geom_tile(aes(x="", y=name, fill=countT2), color = "darkgrey") +
    #geom_raster() +
    #ylim(c(0,20)) +
   scale_fill_gradientn(colours=rev(brewer.pal(n = 7, name ="RdBu")),
                         na.value = "darkred",
                         breaks=c(0,5,10,15,20),
                         labels=c(0,5,10,15,20),
                         limits=c(0,20)) + theme_minimal()
  
  ER + CT1 + CT2
  ggsave(paste(plt, "HLA_",allele, "_error_score_plot.png", sep=''), 
         plot = ER + CT1 + CT2,
         width = 18, height = 6)
  return(HLA_error_s)
}





#rowSums(abs(M_diff))

myb = seq(0,1,by = 0.01)
myc = colorRampPalette(brewer.pal(n = 7, name ="YlGnBu"))(length(myb))
pheatmap(abs(M_diff),
         border_color = 'darkgrey',
         color = myc, 
         breaks = myb,
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         cutree_rows = 3,
         show_colnames = FALSE,
         main = paste('HLA Jaccard distance heatmap'),
         fontsize = 1,
         cellwidth = 1, cellheight = 1,
         fontsize_row = 1)

suppressMessages(library('philentropy'))
#install.packages('philentropy')
Jacard.Matrix_P <- philentropy::distance(abs(M_diff), method = "jaccard")
Jacard.Matrix_H <- philentropy::distance(t(abs(M_tool1)), method = "jaccard")
row.names(Jacard.Matrix_H) = colnames(M_diff)
colnames(Jacard.Matrix_H) = colnames(M_diff)
suppressMessages(library("RColorBrewer"))


pheatmap(Jacard.Matrix_P,
         border_color = 'darkgrey',
         color = myc, 
         breaks = myb,
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         cutree_rows = 3,
         show_colnames = FALSE,
         main = paste('HLA Jaccard distance heatmap'),
         fontsize = 1,
         cellwidth = 1, cellheight = 1,
         fontsize_row = 1) #,
         #filename = paste(dir_enrichR,file,'_JaccardDist.pdf', sep=''))