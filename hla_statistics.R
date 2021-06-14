# ----- import library ---------
suppressMessages(library(dplyr))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(stringr))
suppressMessages(library("gplots"))
suppressMessages(library(tidyr))
suppressMessages(library(viridis))

# ----- set WorkingDirectory and create outputs folders ---------
setwd('/beegfs/scratch/ric.cosr/ric.cosr/Vo_WGS/hla_statistics/')

# ------ input file ------
#file_hla_typing = 'input/fished_all_tops.tsv'
#file_hla_typing = 'input/optitype_results_classI.tsv'
file_hla_typing = 'input/hlavbseq_results.tsv'

# --- create dir for outputs ----
#tool_typing = "optitype"
tool_typing = "hlavbseq"
#tool_typing = "bwakit"
tabl = paste("tables_", tool_typing, "/",  sep ='') 
plt = paste("plots_", tool_typing, "/",sep ='') 

dir.create(tabl, recursive = T, showWarnings = F)
dir.create(plt, recursive = T, showWarnings = F)

# -----  read HLA output from file  ---------
HLA_t = read.table(file = file_hla_typing , header = T, sep = '\t', check.names = F)
HLA_t$Barcode= as.character(HLA_t$Barcode)
WGSsamples = HLA_t$Barcode
row.names(HLA_t) = HLA_t$Barcode

# modify Sample_name_WGS format for matching with Helen residual tables
HLA_t$Sample_name_WGS = str_sub(string = str_remove_all(HLA_t$Sample_name_WGS, pattern = "_"), start = 2, end = -1)

# ------ metadata ---------
# metadata_Vo.xlsx was downloaded from Dropbox on 10/05
# read metadata
metadata_Vo = read.xlsx(xlsxFile = "input/metadata_Vo.xlsx")
# modify Barcodes format to match the HLA table one
metadata_Vo[grep(pattern = "^0", x = metadata_Vo$Barcode), "Barcode"] = str_sub(metadata_Vo$Barcode[grep(pattern = "^0", x = metadata_Vo$Barcode)], 
                                                                                start = 2, end = -1)
# save Barcode as row.names
row.names(metadata_Vo) = metadata_Vo$Barcode
# filter metadata for available WGS samples
metadata_Vo_HLA = metadata_Vo[WGSsamples,]
# add info on WGS sample name to metadatas
metadata_Vo_HLA$WGS = HLA_t$Sample_name_WGS[HLA_t$Barcode %in% row.names(metadata_Vo_HLA)]
#sort metadata alphabetically
metadata_Vo_HLA = metadata_Vo_HLA[order(metadata_Vo_HLA$WGS),]

# ------ residuals -------
# Read input residual sent by Helen on 30/05
# this files are only temporary cause they were calculated on only 215 samples
# even though the analysis is implemented to use this file, 
# the linear model results are ignored for sign HLA selection 
ls_res = list.files(path = "input/", pattern=glob2rx("resid*.txt"))
# add residual to metadata
for (fl in ls_res) {
  print(str_sub(string = fl, start = 1, end = -16))
  f = read.table(paste("input/",fl, sep =''))
  metadata_Vo_HLA$resid = NA
  # correct two sampleID cause these were rerunned in a second run with different ID
  f[f$V1 == "1255WGS2D2","V1"] = "1255WGS31A4merged"
  f[f$V1 == "1255WGS1B11","V1"] = "1255WGS31E4merged"
  metadata_Vo_HLA[metadata_Vo_HLA$WGS %in% f$V1,"resid"] = f$V2 
  colnames(metadata_Vo_HLA)[dim(metadata_Vo_HLA)[2]] <- str_sub(string = fl, start = 1, end = -16)
}

# -------- my_functions -----------
# set of function to analyse the HLA matrix

#' function for importing the HLA matrix and round it to the 2nd HLA digit
#' @param file file with HLA matrix
#' @param sep file columns separator
#' @param rownames define which variable to use as rownames c("Barcode","Sample_name_WGS"),
#'  default Barcode
#' @return HLA matrix rounded at the second digit
#' @examples
#' HLA_2dig = twodigit(file = "input/all_fished_tops_pre.tsv", sep = '\t', rownames = "Sample_name_WGS")
#' @export
twodigit <- function(file, sep = "\t", rownames = "Barcode") {
  HLA_t = read.table(file = file , header = T, sep = sep, check.names = F)
  HLA_t$Barcode = as.character(HLA_t$Barcode)
  HLA_allele= HLA_t[,grep(x = colnames(HLA_t),pattern = "HLA")]
  fx <- function(x) str_sub(string = x, start = 1, end = 11)
  HLA_allele_2ndDIGIT = as.data.frame(sapply(HLA_allele, fx))
  row.names(HLA_allele_2ndDIGIT) <- HLA_t[,rownames]
  colnames(HLA_allele_2ndDIGIT) <- str_replace_all(string = colnames(HLA_allele_2ndDIGIT), pattern = "-", replacement = "")
  colnames(HLA_allele_2ndDIGIT) <- str_replace_all(string = colnames(HLA_allele_2ndDIGIT), pattern = " ", replacement = "_")
  return(HLA_allele_2ndDIGIT)
}

#' Create matrix of HLA with samplesID on the row and HLA-haplotypes on the columns
#' @param file file with the hla assessment to read
#' @param allele  one of the class I allele = c("A", "B", "C")
#' @param rownames define which variable to use as rownames, c("Sample_name_WGS", "Barcode")
#' @param tool tool used for haplotyping c("optitype", "other")
#' @return Matrix subset for one of the allele to use for further analysis
#' @examples
#' M_HLA = Create.HLAMatrix(file = "input/fished_all_tops.tsv", allele = "A",
#'                        sep = '\t', rownames = "Sample_name_WGS")
#' @export
Create.HLAMatrix <- function(file, allele, sep = '\t',rownames = "Barcode", tool) {
  HLA_t = read.table(file = file , header = T, sep = sep, check.names = F)
  HLA_t$Barcode = as.character(HLA_t$Barcode)
  HLA_allele= HLA_t[,grep(x = colnames(HLA_t),pattern = "HLA")]
  
  fx <- function(x) {
    a = strsplit(x, '\\*')
    for (i in 1:length(a)) { 
      a[[i]][2] = str_sub(string = a[[i]][2], start = 1, end = 5)
      if (a[[i]][1] == "n/a") {
        a[[i]] = a[[i]][1]
      } else {
        a[[i]] = paste(a[[i]], collapse = "*")
      }
      }
    return(unlist(a))
    ### str_sub(string = x, start = 1, end = 11)
    }
  HLA_twodigit = as.data.frame(sapply(HLA_allele, fx))
  row.names(HLA_twodigit) <- HLA_t[,rownames]
  colnames(HLA_twodigit) <- str_replace_all(string = colnames(HLA_twodigit), 
                                            pattern = "-", replacement = "")
  colnames(HLA_twodigit) <- str_replace_all(string = colnames(HLA_twodigit), 
                                            pattern = " ", replacement = "_")
  
  index = grep(pattern = paste("HLA",allele,sep=''), x = colnames(HLA_twodigit))
  
  hla.unique = sort(unique(c(HLA_twodigit[,index[1]], HLA_twodigit[,index[2]])))
  samples = row.names(HLA_twodigit)
  # initialize matrix
  M = matrix(0, nrow = length(samples), ncol = length(hla.unique), 
             dimnames = list(sampleID = samples, hla = hla.unique))
  # fill matrix
  for (s in samples) {
    hla.1 = HLA_twodigit[s,index[1]] 
    hla.2 = HLA_twodigit[s,index[2]]
    M[s,hla.1] = M[s,hla.1]+1
    M[s,hla.2] = M[s,hla.2]+1
  }
  if (tool == "optitype") {
    colnames(M) <- paste("HLA",colnames(M), sep ="-")
  }
  return(M)
}

#' function for adding variable to regress to matrix 
#' @param M the HLA count matrix
#' @param metadata metadata dataframe with samples information
#' @param info2add string with columns of metadata you wish to add to \code(M)
#' default info are c("Gender", "positive_swab", "WGS"
#' "Groundtruth_GTA", "Groundtruth_direct_contacts_GTB", "Groundtruth_indirect_contacts_GTC",
#' "Abbot_semiquantitative", "Roche_Total_ICO", "Diasorin_IgG_semiquantitative",
#' "resid_abbot", "resid_roche", "resid_diasorin", "resid_GTB", "resid_GTC")
#' @param clean_output if TRUE will set NA info for Groundtruths to 0 and add 
#' swabs columns with 0 for negative, 1 for positive, NA for missing information
#' @return Matrix \code(M) with added as many columns as the metadata information 
#' specified in \code(info2add).
#' @example 
#' M_HLA_meta = addMetadata(M = M_HLA, metadata = metadata,
#'                          clean_output = TRUE)
#' @export
addMetadata = function(M,
                       metadata,
                       info2add = c("Gender", "positive_swab", "Groundtruth_GTA",
                                    "Groundtruth_direct_contacts_GTB",
                                    "Groundtruth_indirect_contacts_GTC",
                                    "Abbot_semiquantitative", 
                                    "Roche_Total_ICO", 
                                    "Diasorin_IgG_semiquantitative",
                                    "WGS",
                                    "resid_abbot", "resid_roche", "resid_diasorin",
                                    "resid_GTB", "resid_GTC"),
                       clean_output = TRUE) {
  M = M[metadata$Barcode,] # order input matrix as metadata 
  M_hla_meta = as.data.frame(cbind(M, 
                                   metadata[, info2add]))
  if (clean_output) {
    M_hla_meta$Groundtruth_GTA[is.na(M_hla_meta$Groundtruth_GTA)] = 0
    M_hla_meta$Groundtruth_direct_contacts_GTB[is.na(M_hla_meta$Groundtruth_direct_contacts_GTB)] = 0
    M_hla_meta$Groundtruth_indirect_contacts_GTC[is.na(M_hla_meta$Groundtruth_indirect_contacts_GTC)] = 0
    M_hla_meta$swabs = NA
    M_hla_meta$swabs[M_hla_meta$positive_swab] = 1
    M_hla_meta$swabs[!M_hla_meta$positive_swab] = 0
  }
  return(M_hla_meta)
}



#' function for performing linear model on single allele with
#' binomial error distribution
#' @param M is the hla count matrix with metadata produced by CreateMatrix and addMetadata
#' @param var is the name of the column with the variable to regress
#' @param hla.as.factor logical, consider hla genotype as a factor, default FALSE
#' @param allele allele to use, c("A", "B", "C")
#' @param outfile logical, save output on a file, default FALSE
#' @return a dataframe with statistical information (Estimate and canonical pvalue) 
#' for any typed HLA of chosen \code(allele)
#' @example
#' glm_dataframe = hla_glm(M_HLA_meta , var = "swabs", hla.as.factor = FALSE, 
#'                 allele = c("A", "B", "C"), outfile = FALSE)
#' @export
#M = M_hla.complete; var = "swabs"; allele = "A"; hla.as.factor = FALSE
hla_glm = function(M, var = "swabs", hla.as.factor = FALSE, 
                   allele = c("A", "B", "C"), outfile = FALSE) {
  data = M[!is.na(M[,var]),]
  glmHLA_dataframe = data.frame()
  hla = colnames(data)[grepl(colnames(data),pattern = "HLA")]
  for (hla_i in hla) {
    if (hla.as.factor) { # consider hla value as FACTOR
      if (length(levels(factor(data[,hla_i]))) >= 2) { 
        glm_HLAi = glm(as.factor(data[,var]) ~ as.factor(data[,hla_i]), 
                       family = "binomial")
        if (dim(summary(glm_HLAi)$coefficients)[1]==3) {
          df_entry = data.frame(hla = hla_i, 
                                gen1_pvalue = summary(glm_HLAi)$coefficients[2,4],
                                gen2_pvalue = summary(glm_HLAi)$coefficients[3,4])
          
        } else if (dim(summary(glm_HLAi)$coefficients)[1]==2) {
          df_entry = data.frame(hla = hla_i, 
                                gen1_pvalue = summary(glm_HLAi)$coefficients[2,4],
                                gen2_pvalue = NA)
        } else {
          df_entry = data.frame(hla = hla_i, 
                                gen1_pvalue = NA,
                                gen2_pvalue = NA)
        }
      } else {
        df_entry = data.frame(hla = hla_i, 
                              gen1_pvalue = NA,
                              gen2_pvalue = NA)
      }
      glmHLA_dataframe = rbind(glmHLA_dataframe, df_entry)
    }
    # if hla.as.factor = FALSE
    else {
      if (length(levels(factor(data[,hla_i]))) >= 2) { 
      glm_HLAi = glm(as.factor(data[,var]) ~ data[,hla_i], 
                     family = "binomial")
      
      df_entry = data.frame(hla = hla_i, 
                            Estimate = summary(glm_HLAi)$coefficients[2,1],
                            pvalue = summary(glm_HLAi)$coefficients[2,4])
      } else {
        df_entry = data.frame(hla = hla_i, 
                              Estimate = NA,
                              pvalue = NA)
      }
      glmHLA_dataframe = rbind(glmHLA_dataframe, df_entry)
      
    }
  }
  # write results in tsv file, default FALSE will skip write.table
  if (outfile) {
    write.table(glmHLA_dataframe, file = paste(tabl,"statistics_glm_", var,
                                               "_HLAasfactor_",
                                               hla.as.factor,
                                               "_allele_", allele,
                                               ".tsv", sep = ""),
                row.names = FALSE, sep = "\t", quote = FALSE)
  }
  # return dataframe
  return(glmHLA_dataframe)
}

#' function for performing linear model on single allele
#' @param M is the hla count matrix with metadata produced by CreateMatrix and addMetadata
#' @param var is the name of the column with the variable to regress
#' @param hla.as.factor logical, consider hla genotype as a factor, default FALSE
#' @param allele allele to use, c("A", "B", "C")
#' @param outfile logical, save output on a file, default FALSE
#' @return a dataframe with statistical information (Estimate and canonical pvalue) 
#' for any typed HLA of chosen \code(allele)
#' @example
#' lm_dataframe = hla_lm(M_HLA_meta, var = "Abbot_semiquantitative", hla.as.factor = FALSE, 
#'                 allele = c("A", "B", "C"), outfile = FALSE)
#' @export
# M = M_hla.complete; var = "Abbot_semiquantitative"; allele = "A"; hla.as.factor = FALSE
hla_lm = function(M, var = "Abbot_semiquantitative", 
                  hla.as.factor = FALSE, 
                  allele = c("A", "B", "C"),
                  outfile = FALSE) {
  data = M[!is.na(M[,var]),]
  lmHLA_dataframe = data.frame()
  hla = colnames(data)[grepl(colnames(data),pattern = "HLA")]
  for (hla_i in hla) {
    #print(hla_i)
    if (hla.as.factor) {
      if (length(levels(factor(data[,hla_i]))) >= 2) { 
      lm_HLAi = lm(data[,var] ~ as.factor(data[,hla_i]))
      if (dim(summary(lm_HLAi)$coefficients)[1]==3) {
        df_entry = data.frame(hla = hla_i, 
                              gen1_pvalue = summary(lm_HLAi)$coefficients[2,4],
                              gen2_pvalue = summary(lm_HLAi)$coefficients[3,4])
        
      } else if (dim(summary(lm_HLAi)$coefficients)[1]==2) {
        df_entry = data.frame(hla = hla_i, 
                              gen1_pvalue = summary(lm_HLAi)$coefficients[2,4],
                              gen2_pvalue = NA)
      } else {
        df_entry = data.frame(hla = hla_i, 
                              gen1_pvalue = NA,
                              gen2_pvalue = NA)
      }
      } else {
        df_entry = data.frame(hla = hla_i, 
                              gen1_pvalue = NA,
                              gen2_pvalue = NA)
      }
      lmHLA_dataframe = rbind(lmHLA_dataframe, df_entry)
    }
    # if hla.as.factor = FALSE
    else {
      lm_HLAi = lm_HLAi = lm(data[,var] ~ data[,hla_i])
      if (dim(summary(lm_HLAi)$coefficients)[1]==2) {
        df_entry = data.frame(hla = hla_i, 
                              Estimate = summary(lm_HLAi)$coefficients[2,1],
                              Std.Error = summary(lm_HLAi)$coefficients[2,2],
                              pvalue = summary(lm_HLAi)$coefficients[2,4])
      } else {
        df_entry = data.frame(hla = hla_i, 
                              Estimate = NA,
                              Std.Error = NA,
                              pvalue = NA)
      }
      lmHLA_dataframe = rbind(lmHLA_dataframe, df_entry)
    }
  }
  # write results in tsv file
  if (outfile) {
    write.table(lmHLA_dataframe, file = paste(tabl, "statistics_lm_", var,
                                              "_HLAasfactor_",
                                              hla.as.factor,
                                              "_allele_", allele,
                                              ".tsv", sep = ""),
                row.names = FALSE, sep = "\t", quote = FALSE)
  }
  # return dataframe
  return(lmHLA_dataframe)
}


#' function to produce a list of dataframes with the results of glm and lm for 
#' the variable of interest: swabs, GTA, GTB, GTC, antibody level, residuals
#' @param M_hla.complete is the hla count matrix with metadata produced by CreateMatrix and addMetadata.
#' @param allele tag for output files "A", "B" "C"
#' @param countFilter number of min HLA counts use to filter data before FDR evaluation, def is 10
#' @return List of data frames with results of glm and lm for each var of interest
#' @examples 
#' results_all(M_hla.a.complete, allele = "A", countFilter = 10)
#M_hla.complete = M_hla.a.complete; allele = "A"; countFilter = 10
results_all = function(M_hla.complete, 
                       allele = c("A", "B", "C"), 
                       countFilter = 10){
  # ---- Save linear modeling with binomial err distribution ----
  glm_swabs_hla.lin = hla_glm(M = M_hla.complete, var = "swabs", 
                          allele = allele, hla.as.factor = FALSE)
  glm_GTC_hla.lin = hla_glm(M = M_hla.complete, 
                            var = "Groundtruth_indirect_contacts_GTC",
                            allele = allele, hla.as.factor = FALSE)
  glm_GTB_hla.lin = hla_glm(M = M_hla.complete, 
                            var = "Groundtruth_direct_contacts_GTB",
                            allele = allele, hla.as.factor = FALSE)
  glm_GTA_hla.lin = hla_glm(M = M_hla.complete, var = "Groundtruth_GTA",
                            allele = allele, hla.as.factor = FALSE)
  
  # --- filter HLA express in more than countFilter patients
  filter = colSums(M_hla.complete[,grepl(x = colnames(M_hla.complete), 
                                         pattern = "^HLA")]) > countFilter
  HLAexp = colnames(M_hla.complete)[grepl("^HLA", 
                                          colnames(M_hla.complete))][filter]
  dataframe2join = list(swab = glm_swabs_hla.lin, 
                        GTA = glm_GTA_hla.lin, 
                        GTB = glm_GTB_hla.lin, 
                        GTC = glm_GTC_hla.lin)
  dataframe2join_expr = list()
  i = 0
  for (df in dataframe2join) {
    i = i + 1
    df_expr = df[df$hla %in% HLAexp,]
    df_expr$FDR = p.adjust(df_expr$pvalue, method = "BH")
    dataframe2join_expr[[i]] = df_expr
  }
  names(dataframe2join_expr) = names(dataframe2join)
  RESULTS_GLM = Reduce(function(x, y) merge(x, y, by = "hla"), dataframe2join)
  # note that this function give a warning cause the columns of dataframe to merge
  # have the same name - they are therefore renames:
  conames_GLM = crossing(x = factor(names(dataframe2join), 
                                       levels = names(dataframe2join)),
                            y = factor(c("Estimate","pvalue"), 
                                       levels = c("Estimate","pvalue"))) %>% 
    unite("z", x:y, remove = TRUE)
  colnames(RESULTS_GLM) = c("hla", conames_GLM$z)
  
  RESULTS_GLM_expr = Reduce(function(x, y) merge(x, y, by = "hla"), 
                            dataframe2join_expr)
  # note that this function give a warning cause the columns of dataframe to merge
  # have the same name - they are therefore renames:
  conames_GLMexp = crossing(x = factor(names(dataframe2join_expr), 
                                       levels = names(dataframe2join_expr)),
                            y = factor(c("Estimate","pvalue", "FDR"), 
                                       levels = c("Estimate","pvalue", "FDR"))) %>% 
    unite("z", x:y, remove = TRUE)
  colnames(RESULTS_GLM_expr) = c("hla",conames_GLMexp$z)
  
  write.table(x = RESULTS_GLM, 
              file = paste(tabl,"RESULTS_HLA_", 
                           allele,"_Binomial_lm_statistics.tsv", sep = ''), 
              sep = "\t",
              row.names = F, quote = F)
  write.table(x = RESULTS_GLM_expr, 
              file = paste(tabl,"RESULTSexpr",countFilter,"_HLA_", 
                           allele,"_Binomial_lm_statistics.tsv", sep = ''), 
              sep = '\t',
              row.names = F, quote = F)
  
  # ---- Save linear modeling continuous variable ----
  lm_abbott_hla = hla_lm(M = M_hla.complete, var = "Abbot_semiquantitative", 
                           allele = allele, hla.as.factor = FALSE)
  lm_roche_hla = hla_lm(M = M_hla.complete, var = "Roche_Total_ICO", 
                          allele = allele, hla.as.factor = FALSE)
  lm_Diasorin_hla = hla_lm(M = M_hla.complete, var = "Diasorin_IgG_semiquantitative", 
                             allele = allele, hla.as.factor = FALSE)
  dataframe2join = list(Abbott = lm_abbott_hla,
                        Roche = lm_roche_hla,
                        Diasorin = lm_Diasorin_hla)
  i = 0
  dataframe2join_expr = list()
  for (df in dataframe2join) {
    i = i + 1
    df_expr = df[df$hla %in% HLAexp,]
    df_expr$FDR = p.adjust(df_expr$pvalue, method = "BH")
    dataframe2join_expr[[i]] = df_expr
  }
  names(dataframe2join_expr) = names(dataframe2join)
  
  RESULTS_LM = Reduce(function(x, y) merge(x, y, by = "hla"), dataframe2join)
  # note that this function give a warning cause the columns of dataframe to merge
  # have the same name - they are therefore renames:
  conames_LM = crossing(x = factor(names(dataframe2join), 
                                   levels = names(dataframe2join)),
                        y = factor(c("Estimate", "Std.Error","pvalue"), 
                                   levels = c("Estimate", "Std.Error","pvalue"))) %>% 
    unite("z", x:y, remove = TRUE)
  colnames(RESULTS_LM) = c("hla", conames_LM$z)
  # 
  RESULTS_LM_expr = Reduce(function(x, y) merge(x, y, by = "hla"), dataframe2join_expr)
  # note that this function give a warning cause the columns of dataframe to merge
  # have the same name - they are therefore renames:
  conames_LM = crossing(x = factor(names(dataframe2join_expr), 
                                   levels = names(dataframe2join_expr)),
                        y = factor(c("Estimate", "Std.Error","pvalue", "FDR"), 
                                   levels = c("Estimate", "Std.Error","pvalue", "FDR"))) %>% 
    unite("z", x:y, remove = TRUE)
  colnames(RESULTS_LM_expr) = c("hla", conames_LM$z)
  
  write.table(x = RESULTS_LM, 
              file = paste(tabl,"RESULTS_HLA_", allele,
                           "_antibody_linearModel_statistics.tsv", sep =''), 
              sep = "\t",
              row.names = F, quote = F)
  write.table(x = RESULTS_LM_expr, 
              file = paste(tabl,"RESULTSexpr",countFilter,"_HLA_",allele,
                           "_antibody_linearModel_statistics.tsv", sep =''), 
              sep = "\t",
              row.names = F, quote = F)
  
  # ----- Save linear modeling continuous variable - ***RESIDUALS*** ----
  lm_abbott_hla = hla_lm(M = M_hla.complete, var = "resid_abbot", 
                         allele = allele, hla.as.factor = FALSE)
  lm_roche_hla = hla_lm(M = M_hla.complete, var = "resid_roche", 
                        allele = allele, hla.as.factor = FALSE)
  lm_Diasorin_hla = hla_lm(M = M_hla.complete, var = "resid_diasorin", 
                           allele = allele, hla.as.factor = FALSE)
  lm_GTB_hla = hla_lm(M = M_hla.complete, var = "resid_GTB", 
                         allele = allele, hla.as.factor = FALSE)
  lm_GTC_hla = hla_lm(M = M_hla.complete, var = "resid_GTC", 
                        allele = allele, hla.as.factor = FALSE)

  dataframe2join = list(resAbbott = lm_abbott_hla,
                        resRoche = lm_roche_hla,
                        resDiasorin = lm_Diasorin_hla,
                        resGTB = lm_GTB_hla,
                        resGTC = lm_GTC_hla)
  i = 0
  dataframe2join_expr = list()
  for (df in dataframe2join) {
    i = i + 1
    df_expr = df[df$hla %in% HLAexp,]
    df_expr$FDR = p.adjust(df_expr$pvalue, method = "BH")
    dataframe2join_expr[[i]] = df_expr
  }
  names(dataframe2join_expr) = names(dataframe2join)
  
  RESULTS_RS = Reduce(function(x, y) merge(x, y, by = "hla"), dataframe2join)
  # note that this function give a warning cause the columns of dataframe to merge
  # have the same name - they are therefore renames:
  conames_RS = crossing(x = factor(names(dataframe2join), 
                                   levels = names(dataframe2join)),
                        y = factor(c("Estimate", "Std.Error","pvalue"), 
                                   levels = c("Estimate", "Std.Error","pvalue"))) %>% 
    unite("z", x:y, remove = TRUE)
  colnames(RESULTS_RS) = c("hla", conames_RS$z)
  # --- filter HLA express in more than countFilter patients
  RESULTS_RS_expr = Reduce(function(x, y) merge(x, y, by = "hla"), dataframe2join_expr)
  # note that this function give a warning cause the columns of dataframe to merge
  # have the same name - they are therefore renames:
  conames_RS = crossing(x = factor(names(dataframe2join_expr), 
                                   levels = names(dataframe2join_expr)),
                        y = factor(c("Estimate", "Std.Error","pvalue", "FDR"), 
                                   levels = c("Estimate", "Std.Error","pvalue", "FDR"))) %>% 
    unite("z", x:y, remove = TRUE)
  colnames(RESULTS_RS_expr) = c("hla", conames_RS$z)
  
  write.table(x = RESULTS_RS, 
              file = paste(tabl, "RESULTS_HLA_",allele,
                           "_residuals_linearModel_statistics.tsv", sep =''), 
              sep = "\t",
              row.names = F, quote = F)
  write.table(x = RESULTS_RS_expr, 
              file = paste(tabl, "RESULTSexpr",countFilter,"_HLA_",allele,
                           "_residuals_linearModel_statistics.tsv", sep =''),
              sep = "\t",
              row.names = F, quote = F)
  
  # save dataset in list
  
  res = list(binomial = RESULTS_GLM,
             binomial_expr = RESULTS_GLM_expr,
             continue = RESULTS_LM,
             continue_expr = RESULTS_LM_expr,
             residuals = RESULTS_RS,
             residuals_expr = RESULTS_RS_expr)
  return(res)
}

# -------- ** Analysis starts here ** --------

# -------- ***** Class I ***** --------
# -------- HLA-A -------- 
# create HLA matrix from HLA typing 
M_hla.a = Create.HLAMatrix(file = file_hla_typing, allele = "A", 
                           rownames = "Barcode", tool = tool_typing)
# Add metadata
M_hla.a.complete = addMetadata(M_hla.a, metadata = metadata_Vo_HLA, clean_output = TRUE)
# save the results (ignore 6 warning, columns are then renamed in the function)
RES_HLA = results_all(M_hla.complete = M_hla.a.complete, allele = "A", countFilter = 10)

# save all the result lists for all HLA-A in a dataframe
hlaa_all = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = T)]
RESULTS_A_all = Reduce(function(x, y) merge(x, y, by = "hla"), hlaa_all)
RESULTS_A_all_f = RESULTS_A_all[, grep(pattern = "Std.Error", colnames(RESULTS_A_all), invert = T)]
write.table(x = RESULTS_A_all_f, file = paste(tabl,"results_hla_a_allmodels.tsv", sep=''), 
            sep = "\t", quote = F, row.names = F)
# save all the result lists for filtered HLA-A in a dataframe 
hlaa_expr = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = F)]
RESULTS_A_expr = Reduce(function(x, y) merge(x, y, by = "hla"), hlaa_expr)
##### at the moment ignore residuals
RESULTS_A_expr = RESULTS_A_expr[,grep(pattern = "res", x = colnames(RESULTS_A_expr), invert = TRUE)]
# filter HLA with 2 variable with pvalue lowere than .1
atleast_1sign = rowSums(RESULTS_A_expr[,grep(c("pvalue"), 
                                             x = colnames(RESULTS_A_expr))] < 0.05) > 1
signTable_A = RESULTS_A_expr[atleast_1sign,]
Est = signTable_A[,grep(c("Estimate"), x = colnames(signTable_A))]
# protection
Est_p = Est < 0
protection_A = signTable_A[apply(Est_p, 1, function(x) all(x)),"hla"]
# susceptibility
Est_s = Est > 0
susceptibility_A = signTable_A[apply(Est_s, 1, function(x) all(x)),"hla"]
susceptibility_A; protection_A 

# ---- HLA-A plots --------
# ---- heatmap with pheatmap -----
allele = "A"
data = RES_HLA$continue_expr[,c("Abbott_pvalue", "Roche_pvalue", "Diasorin_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$continue_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin")
suppressMessages(library('pheatmap'))
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$continue_expr[,c("Abbott_Estimate","Roche_Estimate", "Diasorin_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$continue_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$Abbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"HLA", allele, "_antibody_pvalue_heatmap.png", sep =''))


data = RES_HLA$binomial_expr[,c("swab_pvalue", "GTA_pvalue", "GTB_pvalue", "GTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$binomial_expr$hla
colnames(data) = c("swab", "GTA", "GTB", "GTC")
ann_data = RES_HLA$binomial_expr[,c("swab_Estimate","GTA_Estimate", "GTB_Estimate", "GTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$binomial_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$swab_Estimate))
colnames(ann_data_2) = c("swab_Coeff","GTA_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  swab_Coeff = mycol,
  GTA_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, 
         filename = paste(plt,"/HLA", allele,"_binomial_pvalue_heatmap.png", sep =''))


data = RES_HLA$residuals_expr[,c("resAbbott_pvalue", "resRoche_pvalue", "resDiasorin_pvalue",
                                 "resGTB_pvalue", "resGTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$residuals_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin", "GTB", "GTC")
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$residuals_expr[,c("resAbbott_Estimate","resRoche_Estimate", "resDiasorin_Estimate",
                                    "resGTB_Estimate", "resGTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$residuals_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$resAbbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"/HLA", allele,"_residuals_pvalue_heatmap.png", sep =''))
dev.off()

# ---- boxplot with ggplot2 -----
# HLA-A*02:01
d = M_hla.a.complete[,c("HLA-A*02:01", "Abbot_semiquantitative", 
                        "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-A*02:01", "Abbot", "Roche", "Diasorin")
d$`HLA-A*02:01` = factor(d$`HLA-A*02:01`)
dd = gather(d, key=antibody, value=value, -`HLA-A*02:01`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-A*02:01`, y = value, 
                     color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-A*02:01") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.A.02.01.png", sep = ''), width = 4, height = 3)

d = M_hla.a.complete[,c("HLA-A*02:01", "resid_abbot", 
                        "resid_roche", "resid_diasorin",
                        "resid_GTB", "resid_GTC")]
colnames(d) = c("HLA-A*02:01", "Res Abbot", "Res Roche", 
                "Res Diasorin", "Res GTB",
                "Res GTC")
d$`HLA-A*02:01` = factor(d$`HLA-A*02:01`)
dd = gather(d, key=antibody, value=value, -`HLA-A*02:01`)
dd$antibody = factor(dd$antibody, levels = c("Res Abbot", "Res Roche", 
                                             "Res Diasorin", "Res GTB",
                                             "Res GTC"))
ggplot(data = dd, aes(x =`HLA-A*02:01`, y = value, 
                      color = antibody, fill = antibody)) + 
  #geom_boxplot(alpha = 0.5, notch=FALSE) + 
  geom_violin(alpha = 0.5) + 
  xlab("allele counts") + 
  ggtitle("HLA-A*02:01") +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values= viridis(5, option = "B", begin = 0, end = 0.8)) + 
  scale_fill_manual(values= viridis(5, option = "B")) + 
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "right")
ggsave(paste(plt,"/boxplot_residuals_HLA.A.02.01.png", sep = ''), width = 8, height = 6.5)

#HLA-A*26:01
d = M_hla.a.complete[,c("HLA-A*26:01", "Abbot_semiquantitative", 
                        "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-A*26:01", "Abbot", "Roche", "Diasorin")
d$`HLA-A*26:01` = factor(d$`HLA-A*26:01`)
dd = gather(d, key=antibody, value=value, -`HLA-A*26:01`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-A*26:01`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-A*26:01") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.A.26.01.png", sep=''), width = 4, height = 3)

d = M_hla.a.complete[,c("HLA-A*26:01", "resid_abbot", 
                        "resid_roche", "resid_diasorin",
                        "resid_GTB", "resid_GTC")]
colnames(d) = c("HLA-A*26:01", "Res Abbot", "Res Roche", 
                "Res Diasorin", "Res GTB",
                "Res GTC")
d$`HLA-A*26:01` = factor(d$`HLA-A*26:01`)
dd = gather(d, key=antibody, value=value, -`HLA-A*26:01`)
dd$antibody = factor(dd$antibody, levels = c("Res Abbot", "Res Roche", 
                                             "Res Diasorin", "Res GTB",
                                             "Res GTC"))

ggplot(data = dd, aes(x =`HLA-A*26:01`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_violin(alpha = 0.5) + 
  xlab("allele counts") + 
  ggtitle("HLA-A*26:01") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values= viridis(5, option = "B", begin = 0, end = 0.8)) + 
  scale_fill_manual(values= viridis(5, option = "B")) + 
  #scale_color_manual(values=c("#999999", "#E69F00", "#FF3300", "magenta", "yellow")) +
  #scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300", "magenta", "yellow")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "right") 
ggsave(paste(plt,"/boxplot_residual_HLA.A.26.01.png", sep = ''), width = 8, height = 6.5)


# ----- mosaic plot ------
d = M_hla.a.complete[,c("HLA-A*26:01", "HLA-A*02:01","HLA-A*24:02","Groundtruth_GTA", 
                        "Groundtruth_direct_contacts_GTB",
                        "Groundtruth_indirect_contacts_GTC", "swabs")]

pdf(paste(plt,"mosaicplot_GTA_HLA.A.26.01.pdf", sep =''), width = 6, height = 6)
mosaicplot(Groundtruth_GTA ~ `HLA-A*26:01`, data = d, 
           main = "HLA-A*26:01", shade = TRUE, xlab = "GTA")
dev.off()


pdf(paste(plt, "mosaicplot_GTA_HLA.A.02.01.pdf", sep =''), width = 6, height = 6)
mosaicplot(Groundtruth_GTA ~ `HLA-A*02:01`, data = d, 
           main = "HLA-A*02:01", shade = TRUE, xlab = "GTA")
dev.off()



# -------- HLA-B --------
# create hla matrix for HLA-B
M_hla.b = Create.HLAMatrix(file = file_hla_typing, allele = "B", 
                           rownames = "Barcode", tool = tool_typing)
# Add metadata
M_hla.b.complete = addMetadata(M_hla.b, metadata = metadata_Vo_HLA, clean_output = TRUE)
# save the results (ignore 6 warning, columns are then renamed in the function)
RES_HLA = results_all(M_hla.complete = M_hla.b.complete, allele = "B", countFilter = 10)

# results dataframe all HLA-B
hlab_all = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = T)]
RESULTS_B_all = Reduce(function(x, y) merge(x, y, by = "hla"), hlab_all)
RESULTS_B_all_f = RESULTS_B_all[, grep(pattern = "Std.Error", colnames(RESULTS_B_all), invert = T)]
write.table(x = RESULTS_B_all_f, file = paste(tabl,"results_hla_b_allmodels.tsv", sep=''), 
            sep = "\t", quote = F, row.names = F)
# results dataframe filtered HLA-B
hlab_expr = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = F)]
RESULTS_B_expr = Reduce(function(x, y) merge(x, y, by = "hla"), hlab_expr)
######### at the moment ignore residuals
RESULTS_B_expr = RESULTS_B_expr[,grep(pattern = "res", x = colnames(RESULTS_B_expr), invert = TRUE)]
atleast_1sign = rowSums(RESULTS_B_expr[,grep(c("pvalue"), 
                                             x = colnames(RESULTS_B_expr))] < 0.1) > 1
signTable_B = RESULTS_B_expr[atleast_1sign,]
Est = signTable_B[,grep(c("Estimate"), x = colnames(signTable_B))]
# protection
Est_p = Est < 0
protection_B = signTable_B[apply(Est_p, 1, function(x) all(x)),"hla"]
# susceptibility
Est_s = Est > 0
susceptibility_B = signTable_B[apply(Est_s, 1, function(x) all(x)),"hla"]
protection_B; susceptibility_B

# ---- HLA-B plots -------
# -------- pheatmap ----------
data = RES_HLA$continue_expr[,c("Abbott_pvalue", "Roche_pvalue", "Diasorin_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$continue_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin")
suppressMessages(library('pheatmap'))
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$continue_expr[,c("Abbott_Estimate","Roche_Estimate", "Diasorin_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$continue_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$Abbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"/HLAB_antibody_pvalue_heatmap.png", sep =''))


data = RES_HLA$binomial_expr[,c("swab_pvalue", "GTA_pvalue", "GTB_pvalue", "GTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$binomial_expr$hla
colnames(data) = c("swab", "GTA", "GTB", "GTC")
ann_data = RES_HLA$binomial_expr[,c("swab_Estimate","GTA_Estimate", "GTB_Estimate", "GTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$binomial_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$swab_Estimate))
colnames(ann_data_2) = c("swab_Coeff","GTA_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  swab_Coeff = mycol,
  GTA_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, 
         filename = paste(plt,"/HLAB_binomial_pvalue_heatmap.png", sep =''))
dev.off()
# --- boxplot ggplot2 --------
d = M_hla.b.complete[,c("HLA-B*51:01", "Abbot_semiquantitative", 
                        "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-B*51:01", "Abbot", "Roche", "Diasorin")
d$`HLA-B*51:01` = factor(d$`HLA-B*51:01`)
dd = gather(d, key=antibody, value=value, -`HLA-B*51:01`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-B*51:01`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-B*51:01") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.B.51.01.png",sep=''), width = 4, height = 3)


d = M_hla.b.complete[,c("HLA-B*51:01", "resid_abbot", 
                        "resid_roche", "resid_diasorin",
                        "resid_GTB", "resid_GTC")]
colnames(d) = c("HLA-B*51:01", "Res Abbot", "Res Roche", 
                "Res Diasorin", "Res GTB",
                "Res GTC")
d$`HLA-B*51:01` = factor(d$`HLA-B*51:01`)
dd = gather(d, key=antibody, value=value, -`HLA-B*51:01`)
dd$antibody = factor(dd$antibody, levels = c("Res Abbot", "Res Roche", 
                                             "Res Diasorin", "Res GTB",
                                             "Res GTC"))
ggplot(data = dd, aes(x =`HLA-B*51:01`, y = value, 
                      color = antibody, fill = antibody)) + 
  #geom_boxplot(alpha = 0.5, notch=FALSE) + 
  geom_violin(alpha = 0.5) + 
  xlab("allele counts") + 
  ggtitle("HLA-B*51:01") +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values= viridis(5, option = "B", begin = 0, end = 0.8)) + 
  scale_fill_manual(values= viridis(5, option = "B")) + 
  #scale_color_manual(values=c("#999999", "#E69F00", "#FF3300", "magenta", "yellow")) +
  #scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300", "magenta", "yellow")) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "right")
ggsave(paste(plt,"/boxplot_residuals_HLA.B.51.01.png",sep=''), width = 8, height = 6.5)

#RES_HLA$residuals_expr[RES_HLA$residuals_expr$hla == "HLA-B*35:02", 
#                      grep(c("pvalue"), x = colnames(RES_HLA$residuals_expr))]
d = M_hla.b.complete[,c("HLA-B*35:02", "Abbot_semiquantitative", 
                        "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-B*35:02", "Abbot", "Roche", "Diasorin")
d$`HLA-B*35:02` = factor(d$`HLA-B*35:02`)
dd = gather(d, key=antibody, value=value, -`HLA-B*35:02`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-B*35:02`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-B*35:02") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.B.35.02.png", sep=''), width = 4, height = 3)


d = M_hla.b.complete[,c("HLA-B*35:02", "resid_abbot", 
                        "resid_roche", "resid_diasorin",
                        "resid_GTB", "resid_GTC")]
colnames(d) = c("HLA-B*35:02", "Res Abbot", "Res Roche", 
                "Res Diasorin", "Res GTB",
                "Res GTC")
d$`HLA-B*35:02` = factor(d$`HLA-B*35:02`)
dd = gather(d, key=antibody, value=value, -`HLA-B*35:02`)
dd$antibody = factor(dd$antibody, levels = c("Res Abbot", "Res Roche", 
                                             "Res Diasorin", "Res GTB",
                                             "Res GTC"))
ggplot(data = dd, aes(x =`HLA-B*35:02`, y = value, 
                      color = antibody, fill = antibody)) + 
  #geom_boxplot(alpha = 0.5, notch=FALSE) + 
  geom_violin(alpha = 0.5) + 
  xlab("allele counts") + 
  ggtitle("HLA-B*35:02") +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values= viridis(5, option = "B", begin = 0, end = 0.8)) + 
  scale_fill_manual(values= viridis(5, option = "B")) + 
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "right")
ggsave(paste(plt,"/boxplot_residuals_HLA.B.35.02.png",sep=''), width = 8, height = 6.5)



# -----  mosaic plots -----
d = M_hla.b.complete[,c("HLA-B*35:02", "HLA-B*51:01","HLA-B*18:01","Groundtruth_GTA", 
                        "Groundtruth_direct_contacts_GTB",
                        "Groundtruth_indirect_contacts_GTC", "swabs")]

pdf(paste(plt,"mosaicplot_GTA_HLA.B.35.02.pdf", sep=''), width = 6, height = 6)
mosaicplot(Groundtruth_GTA ~ `HLA-B*35:02`, data = d, 
           main = "HLA-B*35:02", shade = TRUE, xlab = "GTA")
dev.off()

pdf(paste(plt,"mosaicplot_GTA_HLA.B.51.01.pdf", sep=''), width = 6, height = 6)
mosaicplot(Groundtruth_GTA ~ `HLA-B*51:01`, data = d, 
           main = "HLA-B*51:01", shade = TRUE, xlab = "GTA")
dev.off()

pdf(paste(plt,"mosaicplot_GTA_HLA.B.18.01.pdf", sep=''), width = 6, height = 6)
mosaicplot(Groundtruth_GTA ~ `HLA-B*18:01`, data = d, 
           main = "HLA-B*18:01", shade = TRUE, xlab = "GTA")
dev.off()

# -------- HLA-C -------- 
# create Matrix 
M_hla.c = Create.HLAMatrix(file = file_hla_typing, allele = "C", 
                           rownames = "Barcode", tool = tool_typing)
# add metadat
M_hla.c.complete = addMetadata(M_hla.c, metadata = metadata_Vo_HLA, clean_output = TRUE)
# save the results (ignore 6 warning, columns are then renamed in the function)
RES_HLA = results_all(M_hla.complete = M_hla.c.complete, allele = "C", countFilter = 10)
# results dataframe all HLA-C
hlac_all = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = T)]
RESULTS_C_all = Reduce(function(x, y) merge(x, y, by = "hla"), hlac_all)
RESULTS_C_all_f = RESULTS_C_all[, grep(pattern = "Std.Error", colnames(RESULTS_C_all), invert = T)]
write.table(x = RESULTS_C_all_f, file = paste(tabl,"results_hla_c_allmodels.tsv", sep=''),
            sep = "\t", quote = F, row.names = F)
# results dataframe filtered HLA-C
hlac_expr = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = F)]
RESULTS_C_expr = Reduce(function(x, y) merge(x, y, by = "hla"), hlac_expr)
######### at the moment ignore residuals
RESULTS_C_expr = RESULTS_C_expr[,grep(pattern = "res", x = colnames(RESULTS_C_expr), invert = TRUE)]
atleast_1sign = rowSums(RESULTS_C_expr[,grep(c("pvalue"), 
                                             x = colnames(RESULTS_C_expr))] < 0.1) > 1
signTable_C = RESULTS_C_expr[atleast_1sign,]
Est = signTable_C[,grep(c("Estimate"), x = colnames(signTable_C))]
# protection
Est_p = Est < 0
protection_C = signTable_C[apply(Est_p, 1, function(x) all(x)),"hla"]
# susceptibility
Est_s = Est > 0
susceptibility_C = signTable_C[apply(Est_s, 1, function(x) all(x)),"hla"]
protection_C; susceptibility_C

# ---- HLA-C plots ------
# --- pheatmap -------
data = RES_HLA$continue_expr[,c("Abbott_pvalue", "Roche_pvalue", "Diasorin_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$continue_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin")
suppressMessages(library('pheatmap'))
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$continue_expr[,c("Abbott_Estimate","Roche_Estimate", "Diasorin_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$continue_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$Abbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"/HLAC_antibody_pvalue_heatmap.png", sep =''))




data = RES_HLA$binomial_expr[,c("swab_pvalue", "GTA_pvalue", "GTB_pvalue", "GTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$binomial_expr$hla
colnames(data) = c("swab", "GTA", "GTB", "GTC")
ann_data = RES_HLA$binomial_expr[,c("swab_Estimate","GTA_Estimate", "GTB_Estimate", "GTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$binomial_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$swab_Estimate))
colnames(ann_data_2) = c("swab_Coeff","GTA_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  swab_Coeff = mycol,
  GTA_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, 
         filename = paste(plt,"HLAC_binomial_pvalue_heatmap.png", sep =''))

# ----- boxplot ggplo2 --------
d = M_hla.c.complete[,c("HLA-C*05:01", "Abbot_semiquantitative", 
                        "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-C*05:01", "Abbot", "Roche", "Diasorin")
d$`HLA-C*05:01` = factor(d$`HLA-C*05:01`)
dd = gather(d, key=antibody, value=value, -`HLA-C*05:01`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-C*05:01`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-C*05:01") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.C.05.01.png", sep =''), width = 4, height = 3)

dev.off()

# ----- mosaic ----
d = M_hla.c.complete[,c("HLA-C*05:01","Groundtruth_GTA", 
                        "Groundtruth_direct_contacts_GTB",
                        "Groundtruth_indirect_contacts_GTC", "swabs")]

pdf(paste(plt,"/mosaicplot_GTA_HLA.C.05.01.pdf", sep =''), width = 6, height = 6)
mosaicplot(Groundtruth_GTA ~ `HLA-C*05:01`, data = d, 
           main = "HLA-C*05:01", shade = TRUE, xlab = "GTA")
dev.off()


# -------- ***** Class II ******* --------
# -------- HLA-DQA1 --------
M_hla.dqa1 = Create.HLAMatrix(file = file_hla_typing, allele = "DQA1", 
                              rownames = "Barcode", tool = tool_typing)
# Add metadata
M_hla.dqa1.complete = addMetadata(M_hla.dqa1, metadata = metadata_Vo_HLA, clean_output = TRUE)
# save the results (ignore 6 warning, columns are then renamed in the function)
RES_HLA = results_all(M_hla.complete = M_hla.dqa1.complete, allele = "DQA1", countFilter = 10)

# save all the result lists for all HLA-A in a dataframe
hlaDQA1_all = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = T)]
RESULTS_DQA1_all = Reduce(function(x, y) merge(x, y, by = "hla"), hlaDQA1_all)
RESULTS_DQA1_all_f = RESULTS_DQA1_all[, grep(pattern = "Std.Error", colnames(RESULTS_DQA1_all), invert = T)]
write.table(x = RESULTS_DQA1_all_f, file = paste(tabl,"results_hla_DQA1_allmodels.tsv", sep=''), 
            sep = "\t", quote = F, row.names = F)
# save all the result lists for filtered HLA-A in a dataframe 
hlaDQA1_expr = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = F)]
RESULTS_DQA1_expr = Reduce(function(x, y) merge(x, y, by = "hla"), hlaDQA1_expr)
##### at the moment ignore residuals
RESULTS_DQA1_expr = RESULTS_DQA1_expr[,grep(pattern = "res", x = colnames(RESULTS_DQA1_expr), invert = TRUE)]
# filter HLA with 2 variable with pvalue lowere than .1
atleast_1sign = rowSums(RESULTS_DQA1_expr[,grep(c("pvalue"), 
                                                x = colnames(RESULTS_DQA1_expr))] < 0.1) > 1
signTable_DQA1 = RESULTS_DQA1_expr[atleast_1sign,]
Est = signTable_DQA1[,grep(c("Estimate"), x = colnames(signTable_DQA1))]
# protection
Est_p = Est < 0
protection_DQA1 = signTable_DQA1[apply(Est_p, 1, function(x) all(x)),"hla"]
# susceptibility
Est_s = Est > 0
susceptibility_DQA1 = signTable_DQA1[apply(Est_s, 1, function(x) all(x)),"hla"]
susceptibility_DQA1; protection_DQA1 

# ---- heatmap with pheatmap -----
allele = "DQA1"
data = RES_HLA$continue_expr[,c("Abbott_pvalue", "Roche_pvalue", "Diasorin_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$continue_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin")
suppressMessages(library('pheatmap'))
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$continue_expr[,c("Abbott_Estimate","Roche_Estimate", "Diasorin_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$continue_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$Abbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"HLA", allele, "_antibody_pvalue_heatmap.png", sep =''))


data = RES_HLA$binomial_expr[,c("swab_pvalue", "GTA_pvalue", "GTB_pvalue", "GTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$binomial_expr$hla
colnames(data) = c("swab", "GTA", "GTB", "GTC")
ann_data = RES_HLA$binomial_expr[,c("swab_Estimate","GTA_Estimate", "GTB_Estimate", "GTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$binomial_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$swab_Estimate))
colnames(ann_data_2) = c("swab_Coeff","GTA_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  swab_Coeff = mycol,
  GTA_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, 
         filename = paste(plt,"/HLA", allele,"_binomial_pvalue_heatmap.png", sep =''))


data = RES_HLA$residuals_expr[,c("resAbbott_pvalue", "resRoche_pvalue", "resDiasorin_pvalue",
                                 "resGTB_pvalue", "resGTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$residuals_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin", "GTB", "GTC")
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$residuals_expr[,c("resAbbott_Estimate","resRoche_Estimate", "resDiasorin_Estimate",
                                     "resGTB_Estimate", "resGTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$residuals_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$resAbbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"/HLA", allele,"_residuals_pvalue_heatmap.png", sep =''))
dev.off()


# -------- HLA-DQB1 --------
M_hla.dqb1 = Create.HLAMatrix(file = file_hla_typing, allele = "DQB1", 
                              rownames = "Barcode", tool = tool_typing)
# Add metadata
M_hla.dqb1.complete = addMetadata(M_hla.dqb1, metadata = metadata_Vo_HLA, clean_output = TRUE)
# save the results (ignore 6 warning, columns are then renamed in the function)
RES_HLA = results_all(M_hla.complete = M_hla.dqb1.complete, allele = "DQB1", countFilter = 10)

# save all the result lists for all HLA-DQB1 in a dataframe
hlaDQB1_all = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = T)]
RESULTS_DQB1_all = Reduce(function(x, y) merge(x, y, by = "hla"), hlaDQB1_all)
RESULTS_DQB1_all_f = RESULTS_DQB1_all[, grep(pattern = "Std.Error", colnames(RESULTS_DQB1_all), invert = T)]
write.table(x = RESULTS_DQB1_all_f, file = paste(tabl,"results_hla_DQB1_allmodels.tsv", sep=''), 
            sep = "\t", quote = F, row.names = F)
# save all the result lists for filtered HLA-DQB1 in a dataframe 
hlaDQB1_expr = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = F)]
RESULTS_DQB1_expr = Reduce(function(x, y) merge(x, y, by = "hla"), hlaDQB1_expr)
##### at the moment ignore residuals
RESULTS_DQB1_expr = RESULTS_DQB1_expr[,grep(pattern = "res", x = colnames(RESULTS_DQB1_expr), invert = TRUE)]
# filter HLA with 2 variable with pvalue lowere than .1
atleast_1sign = rowSums(RESULTS_DQB1_expr[,grep(c("pvalue"), 
                                                x = colnames(RESULTS_DQB1_expr))] < 0.1) > 1
signTable_DQB1 = RESULTS_DQB1_expr[atleast_1sign,]
Est = signTable_DQB1[,grep(c("Estimate"), x = colnames(signTable_DQB1))]
# protection
Est_p = Est < 0
protection_DQB1 = signTable_DQB1[apply(Est_p, 1, function(x) all(x)),"hla"]
# susceptibility
Est_s = Est > 0
susceptibility_DQB1 = signTable_DQB1[apply(Est_s, 1, function(x) all(x)),"hla"]
susceptibility_DQB1; protection_DQB1 

# -------- HLA-DQB1 plots --------
# ---- heatmap with pheatmap -----
allele = "DQB1"
data = RES_HLA$continue_expr[,c("Abbott_pvalue", "Roche_pvalue", "Diasorin_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$continue_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin")
suppressMessages(library('pheatmap'))
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$continue_expr[,c("Abbott_Estimate","Roche_Estimate", "Diasorin_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$continue_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$Abbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"HLA", allele, "_antibody_pvalue_heatmap.png", sep =''))


data = RES_HLA$binomial_expr[,c("swab_pvalue", "GTA_pvalue", "GTB_pvalue", "GTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$binomial_expr$hla
colnames(data) = c("swab", "GTA", "GTB", "GTC")
ann_data = RES_HLA$binomial_expr[,c("swab_Estimate","GTA_Estimate", "GTB_Estimate", "GTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$binomial_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$swab_Estimate))
colnames(ann_data_2) = c("swab_Coeff","GTA_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  swab_Coeff = mycol,
  GTA_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, 
         filename = paste(plt,"/HLA", allele,"_binomial_pvalue_heatmap.png", sep =''))


data = RES_HLA$residuals_expr[,c("resAbbott_pvalue", "resRoche_pvalue", "resDiasorin_pvalue",
                                 "resGTB_pvalue", "resGTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$residuals_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin", "GTB", "GTC")
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$residuals_expr[,c("resAbbott_Estimate","resRoche_Estimate", "resDiasorin_Estimate",
                                     "resGTB_Estimate", "resGTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$residuals_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$resAbbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"/HLA", allele,"_residuals_pvalue_heatmap.png", sep =''))
dev.off()

# ---- boxplot with ggplot2 -----
# HLA-DQB1*03:03
d = M_hla.dqb1.complete[,c("HLA-DQB1*03:03", "Abbot_semiquantitative", 
                        "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-DQB1*03:03", "Abbot", "Roche", "Diasorin")
d$`HLA-DQB1*03:03` = factor(d$`HLA-DQB1*03:03`)
dd = gather(d, key=antibody, value=value, -`HLA-DQB1*03:03`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-DQB1*03:03`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-DQB1*03:03") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.DQB1.03.03.png", sep = ''), width = 4, height = 3)


d = M_hla.dqb1.complete[,c("HLA-DQB1*06:03", "Abbot_semiquantitative", 
                           "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-DQB1*06:03", 
                "Abbot", "Roche", "Diasorin")

d$`HLA-DQB1*06:03` = factor(d$`HLA-DQB1*06:03`)
dd = gather(d, key=antibody, value=value, -`HLA-DQB1*06:03`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-DQB1*06:03`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-DQB1*06:03") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.DQB1.06.03.png", sep = ''), width = 4, height = 3)

#"HLA-DQB1*05:01"
d = M_hla.dqb1.complete[,c("HLA-DQB1*05:01", "Abbot_semiquantitative", 
                           "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-DQB1*05:01", "Abbot", "Roche", "Diasorin")
d$`HLA-DQB1*05:01` = factor(d$`HLA-DQB1*05:01`)
dd = gather(d, key=antibody, value=value, -`HLA-DQB1*05:01`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-DQB1*05:01`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-DQB1*05:01") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.DQB1.05.01.png", sep = ''), width = 4, height = 3)



# -------- HLA-DRB1 --------
M_hla.drb1 = Create.HLAMatrix(file = file_hla_typing, allele = "DRB1", 
                              rownames = "Barcode", tool = tool_typing)
# Add metadata
M_hla.drb1.complete = addMetadata(M_hla.drb1, metadata = metadata_Vo_HLA, clean_output = TRUE)
# save the results (ignore 6 warning, columns are then renamed in the function)
RES_HLA = results_all(M_hla.complete = M_hla.drb1.complete, allele = "DRB1", countFilter = 10)

# save all the result lists for all HLA-DRB1 in a dataframe
hlaDRB1_all = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = T)]
RESULTS_DRB1_all = Reduce(function(x, y) merge(x, y, by = "hla"), hlaDRB1_all)
RESULTS_DRB1_all_f = RESULTS_DRB1_all[, grep(pattern = "Std.Error", colnames(RESULTS_DRB1_all), invert = T)]
write.table(x = RESULTS_DRB1_all_f, file = paste(tabl,"results_hla_DRB1_allmodels.tsv", sep=''), 
            sep = "\t", quote = F, row.names = F)
# save all the result lists for filtered HLA-DRB1 in a dataframe 
hlaDRB1_expr = RES_HLA[grep(pattern = "expr", names(RES_HLA), invert = F)]
RESULTS_DRB1_expr = Reduce(function(x, y) merge(x, y, by = "hla"), hlaDRB1_expr)
##### at the moment ignore residuals
RESULTS_DRB1_expr = RESULTS_DRB1_expr[,grep(pattern = "res", x = colnames(RESULTS_DRB1_expr), invert = TRUE)]
# filter HLA with 2 variable with pvalue lowere than .1
atleast_1sign = rowSums(RESULTS_DRB1_expr[,grep(c("pvalue"), 
                                                x = colnames(RESULTS_DRB1_expr))] < 0.1) > 1
signTable_DRB1 = RESULTS_DRB1_expr[atleast_1sign,]
Est = signTable_DRB1[,grep(c("Estimate"), x = colnames(signTable_DRB1))]
# protection
Est_p = Est < 0
protection_DRB1 = signTable_DRB1[apply(Est_p, 1, function(x) all(x)),"hla"]
# susceptibility
Est_s = Est > 0
susceptibility_DRB1 = signTable_DRB1[apply(Est_s, 1, function(x) all(x)),"hla"]
susceptibility_DRB1; protection_DRB1 

# ---- heatmap with pheatmap -----
allele = "DRB1"
data = RES_HLA$continue_expr[,c("Abbott_pvalue", "Roche_pvalue", "Diasorin_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$continue_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin")
suppressMessages(library('pheatmap'))
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$continue_expr[,c("Abbott_Estimate","Roche_Estimate", "Diasorin_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$continue_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$Abbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"HLA", allele, "_antibody_pvalue_heatmap.png", sep =''))


data = RES_HLA$binomial_expr[,c("swab_pvalue", "GTA_pvalue", "GTB_pvalue", "GTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$binomial_expr$hla
colnames(data) = c("swab", "GTA", "GTB", "GTC")
ann_data = RES_HLA$binomial_expr[,c("swab_Estimate","GTA_Estimate", "GTB_Estimate", "GTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$binomial_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$swab_Estimate))
colnames(ann_data_2) = c("swab_Coeff","GTA_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  swab_Coeff = mycol,
  GTA_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, 
         filename = paste(plt,"/HLA", allele,"_binomial_pvalue_heatmap.png", sep =''))


data = RES_HLA$residuals_expr[,c("resAbbott_pvalue", "resRoche_pvalue", "resDiasorin_pvalue",
                                 "resGTB_pvalue", "resGTC_pvalue")]
data = -log10(data)
row.names(data) = RES_HLA$residuals_expr$hla
colnames(data) = c("Abbott", "Roche", "Diasorin", "GTB", "GTC")
minH = 0; maxH=2
myb = seq(minH, maxH, by = 0.01)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
myc <- crp(length(myb))
ann_data = RES_HLA$residuals_expr[,c("resAbbott_Estimate","resRoche_Estimate", "resDiasorin_Estimate",
                                     "resGTB_Estimate", "resGTC_Estimate")]
ann_data_2 = ann_data<0 
ann_data_2[ann_data_2==TRUE] = "protection"; ann_data_2[ann_data_2==FALSE] = "susceptibility"
ann_data_2 = as.data.frame(ann_data_2)
row.names(ann_data_2) = RES_HLA$residuals_expr$hla
mycol = c("#FFCC33","#000000"); names(mycol) = levels(as.factor(ann_data_2$resAbbott_Estimate))
colnames(ann_data_2) = c("A_Coeff","R_Coeff", "D_Coeff", "GTB_Coeff", "GTC_Coeff")
ann_colors = list(
  A_Coeff = mycol,
  R_Coeff = mycol,
  D_Coeff = mycol,
  GTB_Coeff = mycol,
  GTC_Coeff = mycol)
pheatmap(data, annotation_row = ann_data_2, annotation_colors = ann_colors,
         annotation_legend = FALSE, legend = TRUE,
         breaks = myb, color = myc, cellwidth = 30, cellheight = 12, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0,
         treeheight_col = 15, filename = paste(plt,"/HLA", allele,"_residuals_pvalue_heatmap.png", sep =''))
dev.off()

# ---- boxplot with ggplot2 -----
# HLA-DRB1*11:04
d = M_hla.drb1.complete[,c("HLA-DRB1*11:04", "Abbot_semiquantitative", 
                           "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-DRB1*11:04", "Abbot", "Roche", "Diasorin")
d$`HLA-DRB1*11:04` = factor(d$`HLA-DRB1*11:04`)
dd = gather(d, key=antibody, value=value, -`HLA-DRB1*11:04`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-DRB1*11:04`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-DRB1*11:04") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.DRB1.11.04.png", sep = ''), width = 4, height = 3)



# HLA-DRB1*11:01
d = M_hla.drb1.complete[,c("HLA-DRB1*11:01", "Abbot_semiquantitative", 
                           "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA-DRB1*11:01", "Abbot", "Roche", "Diasorin")
d$`HLA-DRB1*11:01` = factor(d$`HLA-DRB1*11:01`)
dd = gather(d, key=antibody, value=value, -`HLA-DRB1*11:01`)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x =`HLA-DRB1*11:01`, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("allele counts") + 
  ggtitle("HLA-DRB1*11:01") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"boxplot_antibody_mean_HLA.DRB1.11.01.png", sep = ''), width = 4, height = 3)




# -------- Protection/Susceptibility table ---------
# Summary
protection_A; susceptibility_A
protection_B; susceptibility_B
protection_C; susceptibility_C

M_hla.a.complete = addMetadata(M_hla.a, metadata = metadata_Vo_HLA, clean_output = TRUE)
M_hla.b.complete = addMetadata(M_hla.b, metadata = metadata_Vo_HLA, clean_output = TRUE)

# option 1
PS_matrix = data.frame(Barcode = row.names(M_hla.a.complete), 
                       HLAA_status = NA,
                       HLAB_status = NA)


allelestatusB = function(x) {
  # x is a patient
  p = sum(M_hla.b.complete[x,protection_B])
  s = sum(M_hla.b.complete[x,susceptibility_B])
  status = p-s
  status
  return(status)
}
allelestatusA = function(x) {
  # x is a patient
  p = sum(M_hla.a.complete[x,protection_A])
  s = sum(M_hla.a.complete[x,susceptibility_A])
  status = p-s
  status
  return(status)
}

PS_matrix$HLAA_status = sapply(PS_matrix$Barcode, function(x) allelestatusA(x))
PS_matrix$HLAB_status = sapply(PS_matrix$Barcode, function(x) allelestatusB(x))
PS_matrix$HLA_sum = (PS_matrix$HLAA_status + PS_matrix$HLAB_status)
PS_matrix$HLA_mean = (PS_matrix$HLAA_status + PS_matrix$HLAB_status)/2

write.xlsx(x = PS_matrix, file = paste(tabl,"/Protection_Susceptibility_Score_Table.xlsx", sep=''),
           row.names = FALSE, asTable = TRUE)

# -------- *** plots Protection/Susceptibility table ***  ---------
M_hla.a.complete$Barcode = row.names(M_hla.a.complete)
MPS = merge(M_hla.a.complete, PS_matrix, by.y = "Barcode")
d = MPS[,c("HLA_sum", "Abbot_semiquantitative", 
                        "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA_status", "Abbot", "Roche", "Diasorin")
d$HLA_status = factor(d$HLA_status)
dd = gather(d, key=antibody, value=value, -HLA_status)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x=HLA_status, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("(HLA-A SPscore + HLA-B SPscore)/2") + 
  ggtitle("Susceptibility Protection Score HLA-A + HLA-B") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"/boxplot_antibody_mean_HLA_AB_status.png", sep =''), width = 5, height = 3)


d = MPS[,c("HLAA_status", "Abbot_semiquantitative", 
           "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA_status", "Abbot", "Roche", "Diasorin")
d$HLA_status = factor(d$HLA_status)
dd = gather(d, key=antibody, value=value, -HLA_status)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x=HLA_status, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("HLA-A SPscore") + 
  ggtitle("Susceptibility Protection Score HLA-A") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"/boxplot_antibody_mean_HLA_A_status.png", sep =''), width = 8, height = 6)

d = MPS[,c("HLAB_status", "Abbot_semiquantitative", 
           "Roche_Total_ICO", "Diasorin_IgG_semiquantitative")]
colnames(d) = c("HLA_status", "Abbot", "Roche", "Diasorin")
d$HLA_status = factor(d$HLA_status)
dd = gather(d, key=antibody, value=value, -HLA_status)
dd$antibody = factor(dd$antibody, levels = c("Abbot", "Roche", "Diasorin"))
ggplot(data = dd, aes(x=HLA_status, y = value, 
                      color = antibody, fill = antibody)) + 
  geom_boxplot(alpha = 0.5, notch=FALSE) + 
  xlab("HLA-B SPscore") + 
  ggtitle("Susceptibility Protection Score HLA-B") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#FF3300")) +
  geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(antibody ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"/boxplot_antibody_mean_HLA_B_status.png", sep =''), width = 8, height = 6)

d = MPS[,c("HLA_sum", "swabs", 
           "Groundtruth_GTA", "Groundtruth_direct_contacts_GTB")]
colnames(d) = c("HLA_status", "swabs", "GTA", "GTB")
d$HLA_status = factor(d$HLA_status)
dd = gather(d, key=VAR, value=value, -HLA_status)
dd$VAR = factor(dd$VAR, levels = c("swabs", "GTA", "GTB"))
ggplot(data = dd, aes(x=HLA_status, y = value, 
                      color = VAR, fill = VAR)) + 
  geom_violin(alpha = 0.5) + 
  xlab("(HLA-A SPscore + HLA-B SPscore)/2") + 
  ggtitle("Susceptibility Protection Score") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("purple","blue", "orange")) +
  scale_fill_manual(values=c("purple", "blue", "orange")) +
  #geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(VAR ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"/violinplot_GT_mean_HLA_AB_status.png", sep =''), width = 8, height = 6)


mosaicplot(GTA ~ HLA_status, data = d, 
           main = "HLA_status", shade = TRUE, xlab = "GTA")

pdf(paste(plt,"/mosaicplot_GTA_HLAsum.pdf", sep =''), width = 6, height = 6)
mosaicplot(GTA ~ HLA_status, data = d, 
           main = "HLA_status", shade = TRUE, xlab = "GTA")
dev.off()


d = MPS[,c("HLAA_status", "swabs", 
           "Groundtruth_GTA", "Groundtruth_direct_contacts_GTB")]
colnames(d) = c("HLA_status", "swabs", "GTA", "GTB")
d$HLA_status = factor(d$HLA_status)
dd = gather(d, key=VAR, value=value, -HLA_status)
dd$VAR = factor(dd$VAR, levels = c("swabs", "GTA", "GTB"))
ggplot(data = dd, aes(x=HLA_status, y = value, 
                      color = VAR, fill = VAR)) + 
  geom_violin(alpha = 0.5) + 
  xlab("HLA-A SPscore") + 
  ggtitle("Susceptibility Protection Score HLA-A") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("purple","blue", "orange")) +
  scale_fill_manual(values=c("purple", "blue", "orange")) +
  #geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(VAR ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"/violinplot_GT_mean_HLA_A_status.png", sep =''), width = 8, height = 6)
mosaicplot(GTA ~ HLA_status, data = d, 
           main = "HLA_status", shade = TRUE, xlab = "GTA")

d = MPS[,c("HLAB_status", "swabs", 
           "Groundtruth_GTA", "Groundtruth_direct_contacts_GTB")]
colnames(d) = c("HLA_status", "swabs", "GTA", "GTB")
d$HLA_status = factor(d$HLA_status)
dd = gather(d, key=VAR, value=value, -HLA_status)
dd$VAR = factor(dd$VAR, levels = c("swabs", "GTA", "GTB"))
ggplot(data = dd, aes(x=HLA_status, y = value, 
                      color = VAR, fill = VAR)) + 
  geom_violin(alpha = 0.5) + 
  xlab("HLA-B SPscore") + 
  ggtitle("Susceptibility Protection Score HLA-B") +
  stat_summary(fun=mean, geom="point", shape=23, size=3) +
  scale_color_manual(values=c("purple","blue", "orange")) +
  scale_fill_manual(values=c("purple", "blue", "orange")) +
  #geom_jitter(shape=20, position=position_jitter(0.2)) +
  facet_grid(VAR ~ ., scales = "free_y") + 
  theme(legend.title = element_blank(), legend.position = "top")
ggsave(paste(plt,"/violinplot_GT_mean_HLA_B_status.png", sep =''), width = 8, height = 6)
mosaicplot(GTA ~ HLA_status, data = d, 
           main = "HLA_status", shade = TRUE, xlab = "GTA")



# ----- Fisher Test and Chisquare Test --------
#' function to perform chi2 e Fisher test on contingency matrix
#' - M_hla.complete is the matrix produced by CreateMatrix with metadata
#' - var is the variable to use for the phenotype (must be a column of M_hla.complete)
#' - countFilter is the number of min HLA counts use to filter data before FDR evalutation
hlaTests = function(M_hla.complete, var, countFilter = 10) {
  hla = colnames(M_hla.complete)[grepl(colnames(M_hla.complete),pattern = "^HLA")]
  hla_statistics_df = data.frame()
  for (hla_i in hla) {
  	# create contingency matrix
    homo = table(M_hla.complete[M_hla.complete[,hla_i]==2,var])
    etero = table(M_hla.complete[M_hla.complete[,hla_i]==1,var])
    null = table(M_hla.complete[M_hla.complete[,hla_i]==0,var])
    CM = as.matrix(rbind(homo, etero, null))
    # perform tests
    chisq <- chisq.test(CM)
    fisher <- fisher.test(CM)
    hla_statistics_df_entry = data.frame(hla = hla_i,
                                         pvalue.chisq = chisq$p.value, 
                                         pvalue.fisher = fisher$p.value)
    hla_statistics_df = rbind(hla_statistics_df, hla_statistics_df_entry)
  }
  filter = colSums(M_hla.complete[,grepl(x = colnames(M_hla.complete), 
                                         pattern = "^HLA")]) > countFilter
  HLAexp = colnames(M_hla.complete)[grepl("^HLA", colnames(M_hla.complete))][filter]
  hla_statistics_df_exp = hla_statistics_df[hla_statistics_df$hla %in% HLAexp,]
  hla_statistics_df_exp$FDR.chisq = p.adjust(hla_statistics_df_exp$pvalue.chisq, method = "BH")
  hla_statistics_df_exp$FDR.fisher = p.adjust(hla_statistics_df_exp$pvalue.fisher, method = "BH")
  return(hla_statistics_df_exp)
}

hlaTests(M_hla.a.complete, var = "Groundtruth_GTA")
hlaTests(M_hla.b.complete, var = "Groundtruth_GTA")
hlaTests(M_hla.c.complete, var = "Groundtruth_GTA")


# -------- save image ----------
save.image(file = "hla.RData")

sessionInfo()

