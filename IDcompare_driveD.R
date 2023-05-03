#Identifiction camparism: TMT16 and hyperplexed

setwd("D:/Hyperplex_manuscript/ID_with_no_SILAC/ID_compare/")
figure_save = "D:/Hyperplex_manuscript/ID_with_no_SILAC/"

#function
read_inputfile = function(PDrunName, HorL){
  file_assign = read_tsv(paste0(PDrunName, "_", HorL, "_InputFiles.txt")) %>%
    separate(`File Name`, into = c("drive","F1","F2","F3","F4", "rawfile","raw"), sep = "[\\\\\\.]" )%>%
    dplyr::mutate(RawFile = rawfile)%>%
    separate(rawfile, into = c("Experiment", "Sample", "Chip", "field", "MSmethod"), sep = "_" )  %>%
    dplyr::mutate(FileID_SILAC = paste0(HorL, `File ID`))
  return(file_assign)
}

read_channel_TMT16 = function(channel_assign, field, file){
  TMTchannels = c("X126","X127N","X127C","X128N","X128C","X129N","X129C","X130N","X130C",
                  "X131N","X131C", "X132N","X132C","X133N","X133C","X134N")
  SILAC_TMTchannels = c("SILAC_126",	"SILAC_127N",	"SILAC_127C",	"SILAC_128N",
                        "SILAC_128C",	"SILAC_129N",	"SILAC_129C",	"SILAC_130N",
                        "SILAC_130C",	"SILAC_131N",	"SILAC_131C",	"SILAC_132N",
                        "SILAC_132C",	"SILAC_133N",	"SILAC_133C",	"SILAC_134N")
  if (str_detect(file, "L")){
    channel_assign = channel_assign %>%
      dplyr::select(Location, all_of(TMTchannels))
  }
  else if (str_detect(file, "H")){
    channel_assign = channel_assign %>%
      dplyr::select(Location, all_of(SILAC_TMTchannels)) 
    colnames(channel_assign) = str_replace_all(colnames(channel_assign), "SILAC_", "X")
  }
  assigned_Channels = channel_assign[which(channel_assign$Location == paste0("F", field)),]
  return(assigned_Channels)
}

read_Cells_in_channel = function(file_name){
  fileIDs = unique(file_name$FileID_SILAC)
  cells_in_channel = data.frame()
  for (file in fileIDs){
    row = which(file_name$FileID_SILAC == file)
    field = file_name[[row, "field"]]
    sample = file_name[[row, "Sample"]]
    fileID = file_name[[row, "File ID"]]
    method = file_name[[row,"Experiment"]]
    if (sample == "HeLa"){
      channel_assign = channel_assign_16_1cell
    } else if (sample == "HeLaDigest"){
      channel_assign = channel_assign_16_1cell
      channel_assign = data.frame (lapply(channel_assign, function(x){
        gsub("HeLa", "HeLaDigest", x)
      }))
    } else if (sample == "HeLaK562"&file_name[[row, "Experiment"]] == "TMT16"){
      channel_assign = channel_assign_16
    } else if (sample == "HeLaK562"&file_name[[row, "Experiment"]] == "Hyperplex32"){
      channel_assign = channel_assign_32
    } else if (sample == "digest"){
      channel_assign = channel_assign_32
      channel_assign = data.frame (lapply(channel_assign, function(x){
        gsub("HeLa", "HeLaDigest", x)
      }))
      channel_assign = data.frame (lapply(channel_assign, function(x){
        gsub("K562", "K562Digest", x)
      }))
      channel_assign = data.frame (lapply(channel_assign, function(x){
        gsub("Control.+Digest", "Control.PBS", x)
      }))
    }
      cell_in_channel = read_channel_TMT16(channel_assign, field, file) %>%
        dplyr::mutate(FileID_SILAC = file, FileID = fileID, Rawfile = file_name[[row, "RawFile"]], 
                      Field = as.numeric(str_replace(Location,"F","")), Method = method)
      cells_in_channel = rbind(cells_in_channel, cell_in_channel)
  }
      return(cells_in_channel)
  }
    
read_protein = function(PDrunName, HorL){
  df_protein = read_tsv(paste0(PDrunName, "_", HorL, "_Proteins.txt")) 
  df_protein = df_protein %>%
    filter(Contaminant == "FALSE")  %>%
    filter(Master == "IsMasterProtein") %>%
    filter(`Protein FDR Confidence: Combined` == "High") 
  colnames(df_protein) = str_replace_all(colnames(df_protein), 
                                         c(", Sample"="",    #delete ", Sample"
                                           "Found in Sample: \\[S[0-9]+\\] " = paste0("ID_", HorL), 
                                           "Found in Sample: F" = paste0("ID_", HorL, "F"), 
                                           "Abundance: " = paste0("Abundance_", HorL), 
                                           ": " = "_X"))
  return(df_protein)
}

find_quan_pecentage = function(channel_list, protein_list){
  i = 1
  pct = data.frame()
  for (channel in channel_list){
    IDcol = paste0("ID_", channel)
    abundancecol = paste0("Abundance_", channel)
    protein_list_by_channel = protein_list %>%
      select(Accession, `# Unique Peptides`, all_of(IDcol), all_of(abundancecol))
    colnames(protein_list_by_channel)[which(colnames(protein_list_by_channel) == IDcol)] = "confidence"
    colnames(protein_list_by_channel)[which(colnames(protein_list_by_channel) == abundancecol)] = "abundance"
    
    protein_list_by_channel_1 = protein_list_by_channel %>%
      filter(`# Unique Peptides` > 0, confidence == "High")
    protein_list_by_channel_2 = protein_list_by_channel %>%
      filter(`# Unique Peptides` > 1, confidence == "High")
    protein_list_quan_1 = protein_list_by_channel_1 %>%
      filter(abundance > 0)
    protein_list_quan_2 = protein_list_by_channel_2 %>%
      filter(abundance > 0)
    quan_pct_uni_1 = nrow(protein_list_quan_1)/nrow(protein_list_by_channel_1)*100
    quan_pct_uni_2 = nrow(protein_list_quan_2)/nrow(protein_list_by_channel_1)*100
    
    pct[i, "Identifier"] = channel
    pct[i, "quan_pct_1"] = quan_pct_uni_1
    pct[i, "quan_pct_2"] = quan_pct_uni_2
    pct[i, "HighConf_uniqPep1plus"] = nrow(protein_list_by_channel_1)
    i = i +1
  }
  return(pct)
}


 #####   

library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)
  
  channel_assign_16 = read_csv("TMT16 channels 2 cell.csv")
  colnames(channel_assign_16)[-1] = paste0("X", colnames(channel_assign_16))[-1]
  channel_assign_16_1cell = read_csv("TMT16 channels-1 cell.csv")
  colnames(channel_assign_16_1cell)[-1] = paste0("X", colnames(channel_assign_16_1cell))[-1]
  channel_assign_32 = read_csv("TMT32 design _ 081621.csv")
  colnames(channel_assign_32)[2:17] = paste0("X", colnames(channel_assign_32))[2:17]
  
  PDrunName = "16and32"   #fill here
  
  file_name_L = read_inputfile(PDrunName,HorL = "L") 
  file_name_H = read_inputfile(PDrunName,HorL = "H") %>%
    filter(Experiment== "Hyperplex32")
  file_name = rbind(file_name_L, file_name_H) %>%
    filter(Sample == "HeLaK562")  
  # In this comparison, only the 2 cell runs are compared. Because in SILAC-TMT, "digest" means 2 cell line digestions;
  # But in TMT only, the "HeLadigest" is only HeLa digestion
  
  rm(file_name_L, file_name_H)
  
  file_name_16 = file_name %>% filter(Experiment == "TMT16")
  file_name_32 = file_name %>% filter(Experiment == "Hyperplex32")

  cells_in_channel_16 = read_Cells_in_channel(file_name_16)
  cells_in_channel_32 = read_Cells_in_channel(file_name_32)
  cells_in_channel = rbind(cells_in_channel_16, cells_in_channel_32)
  
  cells_in_channel_longer = cells_in_channel %>%
    pivot_longer(cols = starts_with("x1"), names_to = "Channel", values_to = "Cell")%>%
    dplyr::mutate(Identifier = paste0(FileID_SILAC, "_", Channel))
  cells_in_channel_longer$method = str_replace_all(cells_in_channel_longer$Method, c("Hyperplex32" = "HyperSCP", "TMT16" = "TMTpro"))
  
  rm(channel_assign_16,channel_assign_16_1cell, channel_assign_32,file_name_16,  file_name_32,
     cells_in_channel_16, cells_in_channel_32)

  proteinL = read_protein(PDrunName, HorL = "L")
  proteinH = read_protein(PDrunName, HorL = "H")
  ProteinL_ID = proteinL %>%
    dplyr::select(Accession, str_subset(colnames(proteinL), "ID_"))
  ProteinH_ID = proteinH %>%
    dplyr::select(Accession, str_subset(colnames(proteinH), "ID_"))
  Protein_ID = dplyr::full_join(ProteinL_ID, ProteinH_ID, by = "Accession")
  rm(ProteinL_ID, ProteinH_ID)
  
  proteinID = data.frame()
  i = 1
  for(identifier in colnames(Protein_ID[,-1])){
    proteinID[i, "Identifier"] = str_replace(identifier, "ID_", "")
    proteinID[i, "Protein_ID"] = sum(Protein_ID[,which(colnames(Protein_ID) == identifier)] == "High", na.rm = TRUE)
    i = i + 1
  }
  
  proteinID = left_join(proteinID, cells_in_channel_longer[,c("Identifier", "Cell", "method")], by = "Identifier")
  proteinID = proteinID[which(!is.na(proteinID$Cell)),]
  
  proteinID_sc = proteinID %>%
    filter(!str_detect(Cell, "Control|Reference|Null"))
  
  ggplot(data = proteinID_sc, aes(x = Cell, y = Protein_ID, fill = method))+
    #geom_violin(na.rm = TRUE)+
    geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
    xlab("Cell Line")+
    ylab("Identified Protein Groups")+
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5,size = rel(1.05), color = "black"),
          axis.text.y = element_text(size = 16),
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.background = element_rect(size=0.5
                                           #, linetype="solid", colour ="black",fill="white"
          )) +
    coord_cartesian(ylim=c(200,650)) +
    theme(legend.position = c(0.78, 0.15))
  #    + geom_jitter(color="black", size=0.4, alpha=0.9)
  
  ggsave(filename = paste0(figure_save, "IDcomparison.png"), width = 6, height = 5, dpi = 500)
  
  
  
  #quantifiable ID percentage based on PD output
  
  all_channels = proteinID_sc$Identifier
  H_channels = all_channels[which(str_detect(all_channels,"H"))]
  L_channels = all_channels[which(str_detect(all_channels,"L"))]
 
  H_quan_pct = find_quan_pecentage(channel_list = H_channels, protein_list = proteinH)
  L_quan_pct = find_quan_pecentage(channel_list = L_channels, protein_list = proteinL)
  
  quan_pct = rbind(H_quan_pct, L_quan_pct)
  
  quan_pct = full_join(quan_pct, proteinID_sc, by = "Identifier") 
  quan_pct[which(str_detect(quan_pct$Identifier, "HF")), "Label"] = "HyperSCP-heavy"
  quan_pct[which(str_detect(quan_pct$Identifier, "LF") & str_detect(quan_pct$method, "HyperSCP")), "Label"] = "HyperSCP-light"
  quan_pct[which(quan_pct$method == "TMTpro"), "Label"] = "TMTpro only"
  
  
  
  QuanPct = ggplot(data = quan_pct, aes(x = Label, y = quan_pct_2, fill = Label))+
    #geom_violin(na.rm = TRUE)+
    geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
    xlab("")+
    ylab("Quantifiable Protein Group (%)")+
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 0.5, vjust = 0.7, size = rel(1.05), color = "black"),
          axis.text.y = element_text(size = 16)) +
    coord_cartesian(ylim=c(50, 100)) +
    theme(legend.position = c(0.78, 0.15))
  
  
  QuanPct = QuanPct + theme(legend.position="none" )
  QuanPct
  
  ggsave(QuanPct, filename = paste0(figure_save, "QuanPct_UNI pep more than 1.png"), width = 5, height = 5, dpi = 500)