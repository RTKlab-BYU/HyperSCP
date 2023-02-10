# functions #

library(readr)
library(ggplot2)
library(gridGraphics)
library(stringr)
library(matrixStats)
library("sva")
library(tidyverse)
library("genefilter")
library(sqldf)
library(ggfortify)
library(ggpubr)
library(ggridges)


# read the input files and file number, field and chip info
read_file_name = function(PDrunName, HorL){
  file_assign = read_tsv(paste0(PDrunName, "_", HorL, "_InputFiles.txt")) %>%
    separate(`File Name`, into = c("drive","year","month","project", "rawfile","raw"), sep = "[\\\\\\.]" )%>%
    dplyr::mutate(RawFile = rawfile)%>%
    separate(rawfile, into = c("chip", "Sample", "field", "Note"), sep = "_" )  %>%
    separate(`RT Range [min]`, into = c("start","end"), sep = " - " ) %>%
    separate(end, into = c("RunTime", "pluspoint"), sep = "\\.") %>%
    dplyr::mutate(Gradient = str_replace_all(RunTime, c("145" = "60min", "185" = "90min", "210" = "120min", "144" = "60min")), 
           FileID_SILAC = paste0(HorL, `File ID`))
return(file_assign)
}

# match channels with columns
read_channel = function(channel_assign, field, file){
  TMTchannels = c("126","127N","127C","128N","128C","129N","129C","130N","130C",
                  "131N","131C", "132N","132C","133N","133C","134N", "134C", "135N")
  SILAC_TMTchannels = c("SILAC_126",	"SILAC_127N",	"SILAC_127C",	"SILAC_128N",
                        "SILAC_128C",	"SILAC_129N",	"SILAC_129C",	"SILAC_130N",
                        "SILAC_130C",	"SILAC_131N",	"SILAC_131C",	"SILAC_132N",
                        "SILAC_132C",	"SILAC_133N",	"SILAC_133C",	"SILAC_134N",  
                        "SILAC_134C", "SILAC_135N")
  if (str_detect(file, "L")){
    channel_assign = channel_assign %>%
      dplyr::select(Location, all_of(TMTchannels))
  }
  else if (str_detect(file, "H")){
    channel_assign = channel_assign %>%
      dplyr::select(Location, all_of(SILAC_TMTchannels)) 
    colnames(channel_assign) = str_replace_all(colnames(channel_assign), "SILAC_", "")
  }
  assigned_Channels = channel_assign[which(channel_assign$Location == field),]
  return(assigned_Channels)
}


#replace cells with cell names:
actual_cell = function(data, Cell1, Cell2){
  data = data.frame (lapply(data, function(x){
    gsub("Cell_1", Cell1, x)
  }))
  data = data.frame (lapply(data, function(x){
    gsub("Cell_2", Cell2, x)
  }))
  return(data)
}

#read protein files, for both ID and abundance
read_protein = function(PDrunName, HorL){
  df_protein = read_tsv(paste0(PDrunName, "_", HorL, "_Proteins.txt")) 
  df_protein = df_protein %>%
    filter(Contaminant == "FALSE")  %>%
    filter(Master == "IsMasterProtein") %>%
    filter(`Protein FDR Confidence: Combined` == "High") 
  colnames(df_protein) = str_replace_all(colnames(df_protein), 
                                         c(", Sample"="",    #delete ", Sample"
                                           "Found in Sample: \\[S[0-9]+\\] " = paste0("ID_", HorL), 
                                           "Found in Sample: F" = paste0("ID_", HorL, "F"), #IF PD3.0 output used 
                                           "Abundance: " = paste0("Abundance_", HorL), 
                                           ": " = "_X"))
  return(df_protein)
}

#read protein ratio by PD output PD3.0
read_protein_ratio = function(PDrunName, HorL){
  df_protein = read_tsv(paste0(PDrunName, "_", HorL, "_Proteins.txt")) 
  df_protein = df_protein %>%
    filter(Contaminant == "FALSE")  %>%
    filter(Master == "IsMasterProtein") %>%
    filter(`Protein FDR Confidence: Combined` == "High") %>%
    dplyr::select(Accession, str_subset(colnames(df_protein),"Abundance Ratio"))
  colnames(df_protein) = str_replace_all(colnames(df_protein), 
                                         c("Abundance Ratio: \\(" = HorL, 
                                           "\\) \\/ \\(.+" = "", 
                                           ", " = "_X"))
  return(df_protein)
}

#read peptide file, for manual normalization, and peptide ID
read_peptide = function(PDrunName, HorL){
  df_peptide = read_tsv(paste0(PDrunName, "_", HorL, "_PeptideGroups.txt")) 
  df_peptide = df_peptide %>%
    filter(Contaminant == "FALSE")  %>%
    filter(Confidence == "High") 
  colnames(df_peptide) = str_replace_all(colnames(df_peptide), 
                                         c(", Sample"="",    #delete ", Sample"
                                           "Found in Sample: \\[S[0-9]+\\] " = paste0("ID_", HorL), 
                                           "Found in Sample: F" = paste0("ID_", HorL, "F"), #IF PD3.0 output used 
                                           "Abundance: " = paste0("Abundance_", HorL), 
                                           ": " = "_X", 
                                         "Annotated Sequence" = "Annotated_Sequence"))
  return(df_peptide)
}

# plot ID 
plot_ID_protein_peptide = function(df, output_protein, output_peptide, figure_save){
  box_protein = ggplot(data = df, aes(x = Cell, y = Protein_ID, fill = Carrier))+
    #geom_violin(na.rm = TRUE)+
    geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
    xlab("Cell Line")+
    ylab("Identified Protein Groups")+
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5,size = rel(1.05), color = "black"),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14))
  #  + geom_jitter(color="black", size=0.4, alpha=0.9)
  ggsave(plot = box_protein, filename = paste0(figure_save,
                                               output_protein), 
         width = 6, height = 5, dpi = 500)
  
  box_peptide = ggplot(data = df, aes(x = Cell, y = Peptide_ID, fill = Carrier))+
    #geom_violin(na.rm = TRUE)+
    geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
    xlab("Cell Line")+
    ylab("Identified Peptides")+
    theme_bw(base_size = 18) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5,size = rel(1.05), color = "black"),
          axis.text.y = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)) 
 #  +geom_jitter(color="black", size=0.4, alpha=0.9)
  ggsave(plot = box_peptide, filename = paste0(figure_save,
                           output_peptide), 
         width = 6, height = 5, dpi = 500)
}

#read psm file
read_PSM = function(PDrunName, HorL){
  df_PSM = read_tsv(paste0(PDrunName, "_", HorL, "_PSMs.txt")) %>%
    filter(Contaminant == "FALSE", Confidence == "High") %>%
    filter(is.na(`Quan Info`)) %>%
    dplyr::mutate(FileID_SILAC = paste0(HorL, `File ID`))
    colnames(df_PSM) = str_replace_all(colnames(df_PSM), "Abundance: ", "X")
  return(df_PSM)
}

#normalization based on column (by median) and row (by mean):
cr_norm = function(df,n){
  if(n > 0){
  df_m = as.matrix(df[, c(-1:-n)])
  for(k in 1:ncol(df_m)){
    df_m[,k] = df_m[,k]/median(df_m[,k], na.rm = T)
  }
  for(k in 1:nrow(df_m)){
    df_m[k,] = df_m[k,]/mean(df_m[k,], na.rm = T)
  }
  df = as.data.frame(cbind(df[, c(1:n)], df_m))
  } else if (n == 0){
    df_m = as.matrix(df)
    for(k in 1:ncol(df_m)){
      df_m[,k] = df_m[,k]/median(df_m[,k], na.rm = T)
    }
    for(k in 1:nrow(df_m)){
      df_m[k,] = df_m[k,]/mean(df_m[k,], na.rm = T)
    }
    df = as.data.frame(df_m)
  }
  return(df)
}

#for log-transfered samples
cr_norm_logged = function(df,n){
  if(n > 0){
    df_m = as.matrix(df[, c(-1:-n)])
      for(k in 1:ncol(df_m)){
        df_m[,k] = df_m[,k]-median(df_m[,k], na.rm = T)
      }
      for(k in 1:nrow(df_m)){
        df_m[k,] = df_m[k,]-mean(df_m[k,], na.rm = T)
      }
      df = as.data.frame(cbind(df[, c(1:n)], df_m))
  } else if (n == 0){
    df_m = as.matrix(df)
    for(k in 1:ncol(df_m)){
      df_m[,k] = df_m[,k]-median(df_m[,k], na.rm = T)
    }
    for(k in 1:nrow(df_m)){
      df_m[k,] = df_m[k,]-mean(df_m[k,], na.rm = T)
    }
    df = as.data.frame(df_m)
  }
  return(df)
}



# Remove a row all of NA
RemoveRowsAllNA = function(df){
  df = df[-which(apply(df,1,function(x) all(is.na(x)))),]
  return(df)
}

# filter column and row with a certen percentage of NAs
#n: how many cols are not numbers from the left
filterNApercent = function(df, RowPct, ColPct, n){
  matrix = as.matrix(df[, c(-1:-n)])
  keep_row = c()
  cols = ncol(matrix)
  for(row in 1:nrow(matrix)){
    RowNAPct = sum(is.na(matrix[row,])) / cols
    if(RowNAPct <= RowPct){
      keep_row = c(keep_row,row)
    }
    }
  matrix = matrix[keep_row,]
  df = df[keep_row,] 
  keep_col = c()
  rows = nrow(matrix)
  for(col in 1:ncol(matrix)){
    ColNAPct = sum(is.na(matrix[,col])) / rows
    if(ColNAPct <= ColPct){
      keep_col = c(keep_col,col)
    }
  }
  df = df[,c(1:n,keep_col+n)]
  return(df)
}

# K - nearest neighbors imputation (this part is from ScoPE2)
hknn = function(dat, k=3){
  # Create a copy of the data, NA values to be filled in later
  dat.imp = dat
  # Calculate similarity metrics for all column pairs (default is Euclidean distance)
  dist.mat = as.matrix( dist(t(dat)) )
  #dist.mat =  1-as.matrix(cor((dat), use="pairwise.complete.obs"))
  #dist.mat = as.matrix(as.dist( dist.cosine(t(dat)) ))
  # Column names of the similarity matrix, same as data matrix
  cnames = colnames(dist.mat)
  # For each column in the data... 
  for(X in cnames){
    # Find the distances of all other columns to that column 
    distances = dist.mat[, X]
    # Reorder the distances, smallest to largest (this will reorder the column names as well)
    distances.ordered = distances[order(distances, decreasing = F)]
    # Reorder the data matrix columns, smallest distance to largest from the column of interest
    # Obviously, first column will be the column of interest, column X
    dat.reordered = dat[ , names(distances.ordered ) ]
    # Take the values in the column of interest
    vec = dat[, X]
    # Which entries are missing and need to be imputed...
    na.index = which( is.na(vec) )
    # For each of the missing entries (rows) in column X...
    for(i in na.index){
      # Find the most similar columns that have a non-NA value in this row
      closest.columns = names( which( !is.na(dat.reordered[i, ])  ) )
      #print(length(closest.columns))
      # If there are more than k such columns, take the first k most similar
      if( length(closest.columns)>k ){
        # Replace NA in column X with the mean the same row in k of the most similar columns
        vec[i] = mean( dat[ i, closest.columns[1:k] ] )
      }
      # If there are less that or equal to k columns, take all the columns
      if( length(closest.columns)<=k ){
        # Replace NA in column X with the mean the same row in all of the most similar columns
        vec[i] = mean( dat[ i, closest.columns ] )
      }
    }
    # Populate a the matrix with the new, imputed values
    dat.imp[,X] = vec
  }
  return(dat.imp)
}

#imputation by different cells:
imp_by_cell = function(df_protein1, k = 3){
  imp_input = df_protein1
  Protein = imp_input[, 1]
  imp_input = as.matrix(imp_input[, -1])
  rownames(imp_input) = Protein 
  sc_imputated = hknn(imp_input, k = k)
  sum(is.na(sc_imputated))
  sc_imputated[is.na(sc_imputated)] = 0
  Protein = rownames(sc_imputated)
  sc_imp_df = data.frame(Protein,sc_imputated)
  return(sc_imp_df)
}

#means for correlation of cells
mean_of_cell = function(df, cell){
  df = as.data.frame(df)
  df_filter = df[,c(1, which(str_detect(colnames(df), cell)))]
  Average = as.data.frame(apply(df_filter[,-1], 1, mean))
  colnames(Average) = cell
  return(Average)
}
#filter one cells for correlation within one cell
Within_Cell = function(df, cell){
  df = as.data.frame(df)
  df_filter = df[,c(1, which(str_detect(colnames(df), cell)))]
  return(df_filter)
}



