#boost:126
#not used: 127C
#reference: 134N
#for carrier-free samples: no 126 and 127C
setwd("E:/Hyperplex/PD2.5_EclipseData/R_processing") # source code location
source("HyperSCP_Function.R", echo=TRUE)
figure_save = "E:/Hyperplex/PD2.5_ExplorisData/Exploris_ByGradient/Figure/" # where to save figures

setwd("E:/Hyperplex/PD2.5_ExplorisData/Exploris_ByGradient") # where the PD output saves

ref_channel = "134N" #  which channel is the reference channel
TMTchannels = c("126","127N","127C","128N","128C","129N","129C","130N","130C",
                "131N","131C", "132N","132C","133N","133C","134N", "134C", "135N")
sc_channels = paste0("X", TMTchannels[-which(TMTchannels %in% c("126","127C",ref_channel))])
ref_channel = paste0("X", ref_channel)
TMTchannels = paste0("X", TMTchannels)

### read and assign channels, report annotation H and L as 2 different files.

channel_assign = read_csv("E:/Hyperplex/PD2.5_ExplorisData/Table S1. 2_cells in channles.csv")
PDrunName = "ExMloading_60min"   #fill here

file_name_L = read_file_name(PDrunName,HorL = "L") 
file_name_H = read_file_name(PDrunName,HorL = "H") 
file_name = rbind(file_name_L, file_name_H)
rm(file_name_L, file_name_H)


fileIDs = unique(file_name$FileID_SILAC)

cells_in_channel = data.frame()
for (file in fileIDs){
  row = which(file_name$FileID_SILAC == file)
  field = file_name[[row, "field"]]
  chip = file_name[[row, "chip"]]
  fileID = file_name[[row, "File ID"]]
  cell_in_channel = read_channel(channel_assign, field, file) %>%
    dplyr::mutate(FileID_SILAC = file, FileID = fileID, Rawfile = file_name[[row, "RawFile"]], 
                  Gradient = file_name[[row, "Gradient"]], 
                  Field = as.numeric(str_replace(Location,"F",""))) %>%
    dplyr::mutate(Carrier = ifelse(Field <= 21, "10 ng", "None")) %>%
    dplyr::select(FileID, FileID_SILAC,  Carrier, Rawfile, Gradient, starts_with("1"))
  if (str_detect(chip,"Chip1|Chip2")){
    cell_in_channel = actual_cell(cell_in_channel, "HeLa", "K562")
  }
  else if(chip == "Chip3"){
    cell_in_channel = actual_cell(cell_in_channel, "HeLa", "A549")
  }
  else if(chip == "Chip4"){
    cell_in_channel = actual_cell(cell_in_channel, "HFL1", "A549")
  }
  cells_in_channel = rbind(cells_in_channel, cell_in_channel)
}
cells_in_channel_longer = cells_in_channel %>%
  pivot_longer(cols = starts_with("x1"), names_to = "Channel", values_to = "Cell")%>%
  dplyr::mutate(Identifier = paste0(FileID_SILAC, "_", Channel), Raw = Rawfile) %>%
  separate(Raw, into=c("Chip","left","Field"), sep = "_") %>%
  select(FileID_SILAC,Carrier, Rawfile, Gradient,Channel,Cell,Identifier, Chip)
rm(cell_in_channel, chip, field, file, fileID, row)
# write_csv(cells_in_channel, paste0("annotation_PDrunName", HorL, ".csv"))

#cells_in_channel_chip4 = cells_in_channel_longer %>% 
#  filter(Chip == "Chip4", Gradient == "60min") %>%
#  mutate(Id_cell = paste0(Identifier, ".", Cell))
#cells_in_channel_chip4 = cells_in_channel_chip4[-which(str_detect(cells_in_channel_chip4$Cell, "NotUsed|Boost|Reference")),]
#Chip4_sample = cells_in_channel_chip4$Id_cell

### read PSM and normalization  ###xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
PSM_L = read_PSM(PDrunName, HorL = "L")
PSM_H = read_PSM(PDrunName, HorL = "H")
PSM = rbind(PSM_L, PSM_H) 
PSM = dplyr::full_join(PSM,cells_in_channel[,c("Carrier","FileID_SILAC", "Gradient")], by = "FileID_SILAC")
rm(PSM_L, PSM_H)


#filter and select the carrier
PSM_noCarrier = PSM %>%
  dplyr::filter(Carrier == "None")
PSM_Carrier = PSM %>%
  dplyr::filter(Carrier == "10 ng")  #reduce possible contaminated channels
PSM_Carrier$CarrierRatio = apply(PSM_Carrier[, sc_channels]/PSM_Carrier$X126, 1, median, na.rm = TRUE)
PSM_Carrier = PSM_Carrier %>%
  dplyr::filter(CarrierRatio <= 0.05)  #reduce possible contaminated channels

# start data including carrier/or not
PSM_all = as.data.frame(PSM_Carrier)   ##change here: "PSM_Carrier"  or "PSM_noCarrier"

# Reporter ion intensity distributions #
PSM_sc_longer1 = PSM_all %>%
  dplyr::mutate(`Annotated Sequence` = paste0(`Annotated Sequence`, "_", Charge)) %>%
  dplyr::select(`Spectrum File`, Carrier, FileID_SILAC,`Annotated Sequence`, `Master Protein Accessions`, all_of(sc_channels))
colnames(PSM_sc_longer1) = c("RawFile", "Carrier", "FileID_SILAC", "Sequence", "Protein", sc_channels)
PSM_sc_longer1 = PSM_sc_longer1 %>%
  pivot_longer(cols = starts_with("X1"), names_to = "Channel", values_to = "Quan")%>%
  dplyr::mutate(Identifier = paste0(FileID_SILAC, "_", Channel)) %>%
  dplyr::mutate(Quan = log(Quan, base = 2))
PSM_sc_longer1[which(is.na(PSM_sc_longer1$Quan)), "Quan"] = 0
PSM_sc_longer1 = dplyr::left_join(PSM_sc_longer1, cells_in_channel_longer[, c("Cell", "Identifier")], by = "Identifier")
Q =  PSM_sc_longer1 %>%
  dplyr::group_by(Cell) %>%
  summarise(Q99 = quantile(Quan, 0.99, na.rm = TRUE)) 
PSM_sc_99p = data.frame()
for(cell in Q$Cell){
  PSM_sc_99p_temp = PSM_sc_longer1 %>%
    dplyr::filter(Cell == cell, Quan <= dplyr::pull(Q[which(Q$Cell == cell), "Q99"]))
  PSM_sc_99p = rbind(PSM_sc_99p, PSM_sc_99p_temp)
  rm(PSM_sc_99p_temp)
}

ggplot(PSM_sc_99p, aes(x = Cell, y = Quan, fill = Cell, color = Cell)) +
  geom_violin(scale = "area") +   theme_bw(base_size = 18) +   coord_cartesian(ylim=c(0,20))+
  ylab("Log2(Reporter Ion Itensity)") +
  xlab("Cell Type") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        legend.position = "None")+
  theme(axis.text.x = element_text(angle = 0, vjust = 1, size = rel(1.4), color = "black"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1,size = rel(1.4), color = "black"))+
  theme(axis.title.x = element_text(angle = 0, vjust = 0.5,size = rel(1.5), color = "black"), 
        axis.title.y = element_text(angle = 90, hjust = 0.5,size = rel(1.5), color = "black"))+
  scale_fill_hue(c=55, l=80)
  
ggsave(filename = paste0(figure_save,
                         "RI wcarrier 60min_area.png"), width = 7, height = 6, dpi = 500) 
# Reporter ion intensity distributions ends #

#normalization using reference channel:
PSM_all[,TMTchannels] = PSM_all[,TMTchannels]/PSM_all[,ref_channel]

# filter the single-cell channels for later data processing:
PSM_sc = PSM_all %>%
  dplyr::mutate(`Annotated Sequence` = paste0(`Annotated Sequence`, "_", Charge)) %>%
  dplyr::select(`Spectrum File`, Carrier, FileID_SILAC,`Annotated Sequence`, `Master Protein Accessions`, all_of(sc_channels))
colnames(PSM_sc) = c("RawFile", "Carrier", "FileID_SILAC", "Sequence", "Protein", sc_channels)
PSM_sc_longer = PSM_sc %>% 
  pivot_longer(cols = starts_with("X1"), names_to = "Channel", values_to = "Quan")%>%
  dplyr::mutate(Identifier = paste0(FileID_SILAC, "_", Channel))
PSM_sc_longer = dplyr::left_join(PSM_sc_longer, cells_in_channel_longer[, c("Cell", "Identifier")], by = "Identifier")
PSM_sc_longer = PSM_sc_longer %>%
  dplyr::filter(!is.na(Protein)) %>%
  dplyr::mutate(find_duplicate = paste0(Sequence, "_", Identifier), Id_cell = paste0(Identifier, ".", Cell)) 

#find out the duplicate sequence from the same channel and then keep the max quan item
#duplicate_name = PSM_sc_longer %>%   group_by(find_duplicate) %>%   summarise(freq = n()) %>%
#  filter(freq>1) %>%   select (find_duplicate)

# This is the same but faster
# Take the mean of the replicates
duplicate_name = sqldf("select find_duplicate, count(find_duplicate) 
                       count from PSM_sc_longer group by find_duplicate having count >1")$find_duplicate
duplicate_data = PSM_sc_longer[PSM_sc_longer$find_duplicate %in% duplicate_name, ]
duplicate_data = duplicate_data %>%
  dplyr::group_by(find_duplicate) %>%
  mutate(Quan = mean(Quan))
duplicate_data = duplicate_data[!duplicated(duplicate_data[,"find_duplicate"]),]
PSM_sc_longer = rbind(PSM_sc_longer[!PSM_sc_longer$find_duplicate %in% duplicate_name, ], duplicate_data)  

rm(duplicate_name,duplicate_data)

########### This is to find the missing value rates between H & L ###################

pep_longer = PSM_sc_longer  %>%
  mutate(pep_find = paste0(Sequence, "__", FileID_SILAC))
pep_list = as.data.frame(unique(pep_longer$pep_find)) %>%
  separate(`unique(pep_longer$pep_find)`, sep = "__", into = c("Sequence", "FileID_SILAC")) %>%
  mutate(count = 1)
pep_list = pep_list %>%
  pivot_wider(names_from = "FileID_SILAC", values_from = count)

#with carrier or not
files = unique(cells_in_channel[which(cells_in_channel$Carrier == "10 ng"),"FileID"])
# if not:
files = unique(cells_in_channel[which(cells_in_channel$Carrier == "None"),"FileID"])
#count
find_missing_rate = data.frame()
i = 1
for (file in files){
  Lfile = paste0("L", file)
  Hfile = paste0("H", file)
  pep_file = pep_list[, c("Sequence", Lfile, Hfile)]
  pep_file[, "NAno"] = is.na(pep_file[,Lfile]) + is.na(pep_file[,Hfile])
  pep_file = pep_file %>% filter(NAno != 2)
  find_missing_rate[i, "File"] = file
  find_missing_rate[i, "Light"] = sum(is.na(pep_file[,Lfile]))/ nrow(pep_file)*100
  find_missing_rate[i, "Heavy"] = sum(is.na(pep_file[,Hfile]))/ nrow(pep_file)*100
  i = i + 1
}

find_missing_rate = find_missing_rate %>%
  pivot_longer(cols = c("Light", "Heavy"), names_to = "label", values_to = "value")

# if combine 60min and 90min separately
MR_60mim = find_missing_rate 
MR_60mim = MR_60mim %>%
  mutate(Gradient = "60 min")
MR_90mim = find_missing_rate
MR_60mim = MR_90mim %>%
  mutate(Gradient = "90 min")
MR = rbind(MR_60mim, MR_90mim)

ggplot(data = MR, aes(x = label, y = value, fill = Gradient))+
  geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
  xlab("SILAC Label")+
  ylab("Missing Rate (%)")+
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5,size = rel(1.1), color = "black"),
        axis.text.y = element_text(angle = 0, hjust = 0.5,size = rel(1.1), color = "black"),
        axis.title.x = element_text(angle = 0, vjust = 0.5,size = rel(1.1), color = "black"), 
        axis.title.y = element_text(angle = 90, hjust = 0.5,size = rel(1.1), color = "black"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 17),
        legend.position = c(0.75, 0.8),
        legend.background = element_rect(fill = "transparent"))
ggsave(filename = paste0(figure_save, "Missing Value 60and90 min wCarrier.png"), 
       width = 5, height = 5, dpi = 500)
# end

ggplot(data = find_missing_rate, aes(x = label, y = value))+
         geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
         xlab("SILAC Label")+
         ylab("Missing Rate (%)")+
         theme_bw(base_size = 18) +
         theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
               axis.text.x = element_text(angle = 0, hjust = 0.5,size = rel(1.05), color = "black"),#legend.position = "top",
               axis.text.y = element_text(size = 16),
               legend.title = element_text(size = 16),
               legend.text = element_text(size = 14))
       #  + geom_jitter(color="black", size=0.4, alpha=0.9)
ggsave(filename = paste0(figure_save, "Missing Value 90 min Carrier.png"), 
              width = 5, height = 5, dpi = 500)




#Select data in Chip4:
#PSM_sc_longer = PSM_sc_longer[which(PSM_sc_longer$Id_cell %in% Chip4_sample),]

#PSM_sc_chip4_output = PSM_sc_longer %>%
#  pivot_wider(names_from = "Id_cell", values_from = "Quan", id_cols = c("Sequence", "Protein"))
#write_tsv(as.data.frame(PSM_sc_chip4_output), "PSM_noCarrier_sc_chip4_all_ref_norm.txt")
#rm(PSM_sc_chip4_output)
################################################################################
# calculate CV -PSM (actually protein level, based on SCoPE2) #
PSM_sc_cv = PSM_sc_longer 
PSM_sc_cv$Quan[PSM_sc_cv$Quan == Inf] = NA
PSM_sc_cv$Quan[PSM_sc_cv$Quan == 0] = NA

PSM_sc_cv = PSM_sc_cv %>% 
  dplyr::group_by(Identifier) %>% 
  dplyr::mutate(cell_median = median(Quan, na.rm= TRUE), norm_cell = Quan/cell_median) %>%    # norm_cell = Reporter ion intensity divided by the median of each cell
  dplyr::group_by(Sequence, RawFile) %>% #FileID_SILAC   RawFile
  dplyr::mutate(norm_raw_seq = Quan / mean(norm_cell, na.rm=TRUE)) %>%  # for each FileID_SILAC, and each sequence, Quan data was normalized by the mean of norm_cell
  dplyr::group_by(Protein, Identifier) %>%
  dplyr::mutate(cv = sd(norm_raw_seq, na.rm=T) / mean(norm_raw_seq, na.rm=T))  %>%
  dplyr::group_by(Identifier) %>% 
  dplyr::mutate(cv_median = median(cv, na.rm=T))

# PSM_sc_cv1 is the kept peptides for cv calculation (has more than 1 normalization value /sequence/rawfile, otherwise no cv)
PSM_sc_cv1 = PSM_sc_cv%>% 
  dplyr::group_by(Protein, Identifier) %>%
  dplyr::mutate(cvn = sum(!is.na(norm_raw_seq)))%>%
  dplyr::filter(cvn > 1)

hist(unique(PSM_sc_cv1$cv_median[PSM_sc_cv1$Cell != "Control"]), col=rgb(0,1,0,1/4), prob=T, breaks=50, main = "X single cells ", xlab="CV")
hist(unique(PSM_sc_cv1$cv_median[PSM_sc_cv1$Cell == "Control"]), col=rgb(1,0,0,1/4), prob=T, add=T, breaks=50)


# Filter out variable wells and controls
cvBar = 0.3
sc_kept = unique(PSM_sc_cv1$Id_cell[PSM_sc_cv1$Cell != "Control" & PSM_sc_cv1$cv_median < cvBar])
control_kept = unique(PSM_sc_cv1$Id_cell[PSM_sc_cv1$Cell == "Control"& PSM_sc_cv1$cv_median > cvBar]) 
sc_total = unique(PSM_sc_cv1$Id_cell[PSM_sc_cv1$Cell != "Control"])
control_total = unique(PSM_sc_cv1$Id_cell[PSM_sc_cv1$Cell == "Control"])

PSM_sc_cv1$control = "sc"
PSM_sc_cv1$control[PSM_sc_cv1$Cell == "Control"] = "Control"

# Plot CV
Colors = c(  "#fdae61", "#2c7bb6", "#a6dba0", "#d7191c")
plot_cv = ggplot(data = PSM_sc_cv1, aes(x = cv_median)) + 
  geom_density(aes(fill = control, alpha = 0.6), adjust = 3) + 
  theme_pubr( ) +
  scale_fill_manual(values = Colors[c(1,2)]) +
  xlab("Quantification variability") + ylab("Density") + rremove("y.ticks") + rremove("y.text") +
  font("xylab", size=25) +
  font("x.text", size=20) +
  #  font("y.text", size=20) +
  coord_cartesian(xlim=c(0.0,0.7))+
  annotate("text", x=0.07, y= 16, label=paste0(length(sc_kept)," Cells"), size=8, color=Colors[c(2)])+
  annotate("text", x=0.55, y= 14, label=paste0(length(control_kept)," Controls"), size=8, color=Colors[c(1)])+
  annotate("text", x=0.51, y= 16, label=paste0(length(sc_total) -length(sc_kept)," Cells"), size=8, color=Colors[c(2)])+
  annotate("text", x=0.08, y= 14, label=paste0(length(control_total) - length(control_kept)," Controls"), size=8, color=Colors[c(1)])+
  rremove("legend") + 
  geom_vline(xintercept=cvBar, lty=2, size=1.5, color="gray50")
plot_cv
ggsave(plot_cv, filename = paste0(figure_save,
                                  "CVfilter_Carrier_90min_chip4.png"), 
       width = 7, height = 5, dpi = 500)

#build longer data frame
#PSM_sc_unmelt = dcast(PSM_sc_longer, Sequence ~ Id_cell, value.var = "Quan", fill=NA)
PSM_sc_wider = PSM_sc_longer %>%
  dplyr::filter(!is.na(Protein)) %>%
  dplyr::select(Sequence, Protein, Quan, Id_cell) %>%
  pivot_wider(names_from = "Id_cell", values_from = "Quan")
PSM_sc_wider[PSM_sc_wider == 0] = NA
PSM_sc_wider[PSM_sc_wider == Inf] = NA
PSM_sc_wider[PSM_sc_wider == -Inf] = NA

# filter based on the cv-kept single cells
PSM_sc_sconly = PSM_sc_wider[,c(1,2, which(colnames(PSM_sc_wider) %in% sc_kept))]
PSM_sc_scncon = PSM_sc_wider[,c(1,2,which(colnames(PSM_sc_wider) %in% c(sc_kept,control_kept)))]

# histgram of sc data 
par(mfrow=c(3,3))
# normalized by ref
scPSM_by_ref = PSM_sc_sconly
scPSM_by_ref_m = as.matrix(scPSM_by_ref[,-c(1,2)])
hist(as.matrix(scPSM_by_ref_m), breaks="FD", xlim=c(-2,2))
rm(scPSM_by_ref_m)
# normalized by column (by median) and row (by mean)
scPSM_by_CR = cr_norm(scPSM_by_ref,2)
scPSM_by_CR_m = as.matrix(scPSM_by_CR[,-c(1,2)])
hist(scPSM_by_CR_m, breaks="FD", xlim=c(-2,2))
rm(scPSM_by_CR_m)
# filter by missing value 
scPSM_filter_NA = filterNApercent(scPSM_by_CR, 0.95, 0.99, n=2) #n: how many cols are not numbers from the left
scPSM_filter_NA_m = as.matrix(scPSM_filter_NA[,-c(1,2)])
hist(scPSM_filter_NA_m, breaks="FD", xlim=c(-2,2))
rm(scPSM_filter_NA_m)

# collapse to protein level by mean 
scPSM_longer = scPSM_filter_NA %>% 
  dplyr::filter(!str_detect(Protein, "\\;")) %>%    # delete non-unique peptides?
  pivot_longer(cols = starts_with("LF")|starts_with("HF"), names_to = "Id_cell", values_to = "Quan")
scPSM_longer2 = scPSM_longer %>%
  filter(!is.na(Quan))
scPSM_longer2 =  scPSM_longer2 %>%
  mutate(ProteinIdentifier = paste0(Protein, "-", Id_cell))%>%
  dplyr::group_by(ProteinIdentifier) %>% 
  dplyr::summarize(mean = mean(Quan, na.rm=T), freq = n()) 
scPSM_longer2 = scPSM_longer2 %>%
 # dplyr::filter( freq > 1) %>%                  #  if keep >=2 unique peptides
  separate(ProteinIdentifier, into = c("Protein", "Id_cell"), sep = "-" ) %>%
  select("Protein", "Id_cell", "mean") 
scProtein_protein = scPSM_longer2 %>%
  pivot_wider(names_from = "Id_cell", values_from = "mean")
rm(scPSM_longer, scPSM_longer2)

# log2 transform
scPSM_log2 = scProtein_protein
scPSM_log2[, -1] = log(scPSM_log2[, -1], base = 2)
scPSM_log2_m = as.matrix(scPSM_log2[, -1])
hist(scPSM_log2_m, breaks="FD", xlim=c(-2,2))
rm(scPSM_log2_m)

#write_tsv(scPSM_log2, "PSM_Carrier_sc_chip4_norm_log2.txt")

# filter proteins
# all data
scProtein = filterNApercent(scPSM_log2, 0.9, 0.99, n=1) #pre-filter, not necessary

#by group - cells, keep the protein at least 60% valid in at least a group
group_HeLa = filterNApercent(Within_Cell(scProtein, "HeLa"), 0.4, 0.99, n = 1)
group_K562 = filterNApercent(Within_Cell(scProtein, "K562"), 0.4, 0.99, n = 1)
group_A549 = filterNApercent(Within_Cell(scProtein, "A549"), 0.4, 0.99, n = 1)
group_HFL1 = filterNApercent(Within_Cell(scProtein, "HFL1"), 0.4, 0.99, n = 1)

scProtein0 = full_join(group_HeLa,group_K562, by = "Protein")
scProtein0 = full_join(scProtein0,group_HFL1, by = "Protein")
scProtein0 = full_join(scProtein0,group_A549, by = "Protein")
# scProtein0 = full_join(group_HFL1,group_A549, by = "Protein")

#imputation of proteins found in each cell types, KNN, by different cell types
imp_HeLa = imp_by_cell(group_HeLa, k=5)
imp_K562 = imp_by_cell(group_K562, k=5)
imp_A549 = imp_by_cell(group_A549, k=5)
imp_HFL1 = imp_by_cell(group_HFL1, k=5)

scProtein_imp = full_join(imp_HeLa,imp_K562, by = "Protein")
scProtein_imp = full_join(scProtein_imp,imp_HFL1, by = "Protein")
scProtein_imp = full_join(scProtein_imp,imp_A549, by = "Protein")
# scProtein_imp = full_join(imp_A549,imp_HFL1, by = "Protein") #chip4 only

#write_tsv(scProtein, "Protein_Carrier_sc_chip4_norm_log2.txt")

#Quantifiable protein ID
par(mfrow=c(1,1))
proteinID_PSM = data.frame()
i = 1
for (colname in colnames(scProtein0)[-1]){
  proteinID_PSM[i,"sample"] = colname
  proteinID_PSM[i,"ProteinID"] = sum(!is.na(scProtein0[, colname]))
  i = i + 1
}
proteinID_PSM = proteinID_PSM %>%
  separate(sample, into = c("file", "Cell"), sep = "\\.")
ggplot(data = proteinID_PSM, aes(x = Cell, y = ProteinID, fill = Cell)) +
  #geom_violin(na.rm = TRUE)+
  geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
  xlab("Cell Line")+
  ylab("Quantifiable Protein Groups")+
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = rel(0.9), color = "black"),
        axis.text.y = element_text(size = 16)) + rremove("legend") 
#  + geom_jitter(color="black", size=0.4, alpha=0.9)
ggsave(filename = paste0(figure_save, "QuanProteinID_noCarrier_chip4_60vv.png"), 
       width = 6, height = 5, dpi = 500)



#filter protein list found in different cell types
scProtein = filterNApercent(scProtein_imp, 0.01, 0.99, n = 1) 

scProtein_m = as.matrix(scProtein[, -1])
rownames(scProtein_m) = scProtein[,1]
hist(scProtein_m, breaks="FD", xlim=c(-2,2))
#found in one cell line only
HeLa_special = group_HeLa[which(!group_HeLa[,"Protein"] %in% scProtein[,"Protein"]),"Protein"]
K562_special = group_K562[which(!group_K562[,"Protein"] %in% scProtein[,"Protein"]),"Protein"] 
A549_special = group_A549[which(!group_A549[,"Protein"] %in% scProtein[,"Protein"]),"Protein"] 
HFL1_special = group_HFL1[which(!group_HFL1[,"Protein"] %in% scProtein[,"Protein"]),"Protein"] 

rm(group_HeLa, group_K562, group_A549, group_HFL1,
   imp_HeLa, imp_K562, imp_A549, imp_HFL1)

# Batch correction with ComBat
Id_cell = as.data.frame(colnames(scProtein_m)) %>%
  separate(`colnames(scProtein_m)`, into = c("File","Channel","Cell"), sep = "[_\\.]" )
Combat_mod = as.data.frame(Id_cell$Cell)
colnames(Combat_mod) = "Cell"
Combat_mod = model.matrix(~as.factor(Cell), data = Combat_mod)

ComBat_norm = ComBat(scProtein_m, batch = Id_cell[,"File"], mod = Combat_mod)

hist(scProtein_m, breaks="FD", xlim=c(-2,2))
hist(ComBat_norm, breaks="FD", xlim=c(-2,2))

## PCA ##############################
#use the imputated data for PCA and volcano (that means, proteins found in 1/2/3 cell types only not included)

# If H or L doesn't need to be shown in figure
protein_List = as.data.frame(t(ComBat_norm))
Id_cell = as.data.frame(rownames(protein_List)) %>%
  separate(`rownames(protein_List)`, into = c("File","Channel","Cell"), sep = "[_\\.]" ) 

# no frame
protein_List = cbind(Cell = Id_cell$Cell, protein_List)
protein_List.pca = prcomp(protein_List[, -1])
summary(protein_List.pca)

autoplot(protein_List.pca, data = protein_List, colour = 'Cell', frame = T, frame.type = 'norm') + #shape = FALSE, label.size = 3
  theme_bw(base_size = 18, ) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size=20),
        legend.title = element_text(size=24))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = rel(1.7), color = "black"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1,size = rel(1.7), color = "black"))+
  theme(axis.title.x = element_text(angle = 0, vjust = 0.5,size = rel(1.5), color = "black"), 
        axis.title.y = element_text(angle = 90, hjust = 0.5,size = rel(1.5), color = "black"))
ggsave(filename = paste0(figure_save,
                         "PCA carrier 90min_.png"), width = 9.4, height = 6, dpi = 500)         

### if H or L need to be shown in figure
protein_List = as.data.frame(t(ComBat_norm))
Id_cell = as.data.frame(rownames(protein_List)) %>%
  separate(`rownames(protein_List)`, into = c("File","Channel","Cell"), sep = "[_\\.]" ) %>%
  separate(File, into = c("SILAC","File"),sep = "F" ) %>%
  dplyr::mutate(Label = str_replace_all(SILAC, c("L" = "", "H" = "SILAC_"))) %>%
  dplyr::mutate(Cell = paste0(Label, Cell))

# with frame 
protein_List = cbind(Cell = Id_cell$Cell, protein_List)
protein_List.pca = prcomp(protein_List[, -1])
summary(protein_List.pca)

cPalette <- c("#a50026", "#1b7837", "#6baed6", "#542788", "#f46d43", "#41ae76", "#08519c", "#dd3497")
autoplot(protein_List.pca, data = protein_List, colour = 'Cell', scale = 0,frame = T, frame.type = 'norm')+
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size=20),
        legend.title = element_text(size=24))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = rel(2), color = "black"))+
  theme(axis.text.y = element_text(angle = 0, hjust = 1,size = rel(2), color = "black"))+
  theme(axis.title.x = element_text(angle = 0, vjust = 0.5,size = rel(1.7), color = "black"), 
        axis.title.y = element_text(angle = 90, hjust = 0.5,size = rel(1.7), color = "black"))
 # + scale_fill_manual(values=cPalette) +
#  scale_colour_manual(values=cPalette)

ggsave(filename = paste0(figure_save,
                         "PCA nocarrier Chip4_HandL.png"), width = 10, height = 6, dpi = 500)


######################################################################
# correlation

#correlation can use the data with unique proteins found only in 1/2/3 cell lines scProtein0
#however, the combat normalization has to be used to eliminate the batch effect, which requires no NAs

# Replicates

Accession = as.data.frame(rownames(ComBat_norm))
names(Accession) = "Protein"
dfComBat_norm = cbind(Accession, as.data.frame(ComBat_norm))
# Replicates
cells = c("HeLa", "K562","HFL1","A549")
scProtein_mean = as.data.frame(dfComBat_norm$Protein); colnames(scProtein_mean) = "Protein"
for (Cell in cells){
  Within = Within_Cell(dfComBat_norm, cell = Cell)
  Average_H = mean_of_cell(Within, cell = "HF[1-9]+"); colnames(Average_H) = paste0("SILAC ", Cell)
  Average_L = mean_of_cell(Within, cell = "LF[1-9]+"); colnames(Average_L) = Cell
  scProtein_mean = cbind(scProtein_mean, Average_H, Average_L)
}

rm(Cell,Within,Average_H,Average_L)

library(GGally)
ggpairs(scProtein_mean[,-1], upper = list(continuous = wrap(ggally_cor, method = "spearman",size= 4, stars = F)), 
        lower = list(continuous = wrap("points",size = 0.2))) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = paste0(figure_save,
                         "Cor_nocarrier_chip4.png"), width = 7, height = 7, dpi = 500)


library(corrplot)
library(RColorBrewer)
M <-cor(scProtein_mean[,-1])
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=8, name="RdYlBu"))

# chip 4 : volcano

ComBat_norm_chip4 = as.data.frame(ComBat_norm)

HFL1_col = which(str_detect(colnames(ComBat_norm_chip4), "HFL1"))
A549_col = which(str_detect(colnames(ComBat_norm_chip4), "A549"))

Chip4_dif = ComBat_norm_chip4 %>%
  mutate(average.HFL1 = apply(., 1, function(x) mean(x[HFL1_col])),
         average.A549 = apply(., 1, function(x) mean(x[A549_col])),
         difference = average.A549 - average.HFL1,
         pvalue = apply(., 1, function(x) t.test(x[HFL1_col], x[A549_col], var.equal = TRUE)$p.value),
         padjust = p.adjust(pvalue, method = "BH", n = length(pvalue)),
         nlogpadjust = -log(padjust, base = 10),
         nlogp = -log(pvalue, base = 10))
Chip4_dif$diffexpressed = "F"
Chip4_dif$diffexpressed[Chip4_dif$difference > 0.6 & Chip4_dif$pvalue < 0.05] = "up"
Chip4_dif$diffexpressed[Chip4_dif$difference < -0.6 & Chip4_dif$pvalue < 0.05] = "down"
Chip4_dif$BHdiffexpressed = "F"
Chip4_dif$BHdiffexpressed[Chip4_dif$difference > 0.6 & Chip4_dif$padjust < 0.05] = "up"
Chip4_dif$BHdiffexpressed[Chip4_dif$difference < -0.6 & Chip4_dif$padjust < 0.05] = "down"
Chip4_dif$label = NA
Chip4_dif$label[Chip4_dif$diffexpressed != "F"| Chip4_dif$BHdiffexpressed != "F"] = 
  rownames(Chip4_dif)[Chip4_dif$diffexpressed != "F"| Chip4_dif$BHdiffexpressed != "F"]
Chip4_dif$Protein = rownames(Chip4_dif)

#functions:
proteinlist_L = read_protein(PDrunName, HorL = "L") %>%
  select(Accession, Description, `Biological Process`, `Biological Process`, `Molecular Function`, WikiPathways, `Gene ID`)
proteinlist_H = read_protein(PDrunName, HorL = "H") %>%
  select(Accession, Description, `Biological Process`, `Biological Process`, `Molecular Function`, WikiPathways, `Gene ID`)
proteinlist = full_join(proteinlist_L, proteinlist_H, by = "Accession")
rm(proteinlist_L, proteinlist_H)

Chip4_dif_output = left_join(Chip4_dif,proteinlist, by = c("Protein" = "Accession") )

write_tsv(Chip4_dif_output, "Protein_Carrier_sc_Chip4_different.txt")

only_A549 = proteinlist[which(proteinlist$Accession %in% A549_special),] 
only_HFL1 = proteinlist[which(proteinlist$Accession %in% HFL1_special),]
write_tsv(only_A549, "Chip4_only_A549.txt")
write_tsv(only_HFL1, "Chip4_only_HFL1.txt")


library(ggrepel)
ggplot(data = Chip4_dif, aes(x = difference, y = nlogpadjust, col = diffexpressed, label = label))+
  geom_point()+
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept= -log10(0.05), col="red")+
  scale_color_manual(values=c("blue", "black", "red"))+
  geom_text_repel() + 
  rremove("legend") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = rel(1.1), color = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1,size = rel(1.1), color = "black"),
        axis.title.x = element_text(angle = 0, vjust = 0.5,size = rel(1.1), color = "black"), 
        axis.title.y = element_text(angle = 90, hjust = 0.5,size = rel(1.1), color = "black"))+
        xlab("Log2(fold change A549/HFL1)") + ylab("-log(p-value)") 
ggsave(filename = paste0(figure_save, "Chip4Volcano_carrier_BHadjust.png"), width = 7, height = 5, dpi = 500)
#ggsave(filename = paste0(figure_save,"Chip4Volcano_noLabel.png"), width = 5, height = 4, dpi = 500)



### Protein and peptide ID from Protein file ###########################################################

# Protein ID identified: 

get_ID = function(PDrunName, cell.names){
  #generate a list of files:
  file_name_L = read_file_name(PDrunName,HorL = "L") 
  file_name_H = read_file_name(PDrunName,HorL = "H") 
  file_name = rbind(file_name_L, file_name_H)
  rm(file_name_L, file_name_H)
  fileIDs = unique(file_name$FileID_SILAC)
  cells_in_channel = data.frame()
  for (file in fileIDs){
    row = which(file_name$FileID_SILAC == file)
    field = file_name[[row, "field"]]
    chip = file_name[[row, "chip"]]
    fileID = file_name[[row, "File ID"]]
    cell_in_channel = read_channel(channel_assign, field, file) %>%
      dplyr::mutate(FileID_SILAC = file, FileID = fileID, Rawfile = file_name[[row, "RawFile"]], 
                    Gradient = file_name[[row, "Gradient"]], 
                    Field = as.numeric(str_replace(Location,"F",""))) %>%
      dplyr::mutate(Carrier = ifelse(Field <= 21, "10 ng", "None")) %>%
      dplyr::select(FileID, FileID_SILAC,  Carrier, Rawfile, Gradient, starts_with("1"))
    if (str_detect(chip,"Chip1|Chip2")){
      cell_in_channel = actual_cell(cell_in_channel, "HeLa", "K562")
    }
    else if(chip == "Chip3"){
      cell_in_channel = actual_cell(cell_in_channel, "HeLa", "A549")
    }
    else if(chip == "Chip4"){
      cell_in_channel = actual_cell(cell_in_channel, "HFL1", "A549")
    }
    cells_in_channel = rbind(cells_in_channel, cell_in_channel)
  }
  cells_in_channel_longer = cells_in_channel %>%
    pivot_longer(cols = starts_with("x1"), names_to = "Channel", values_to = "Cell")%>%
    dplyr::mutate(Identifier = paste0(FileID_SILAC, "_", Channel), Raw = Rawfile) %>%
    separate(Raw, into=c("Chip","left","Field"), sep = "_") %>%
    select(FileID_SILAC,Carrier, Rawfile, Gradient,Channel,Cell,Identifier, Chip)
  
  #read protein ID identified:
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
  #read peptide ID identified:
  peptideL = read_peptide(PDrunName, HorL = "L")
  peptideH = read_peptide(PDrunName, HorL = "H")
  peptideL_ID = peptideL %>%
    dplyr::select(Annotated_Sequence, str_subset(colnames(peptideL), "ID_"))
  peptideH_ID = peptideH %>%
    dplyr::select(Annotated_Sequence, str_subset(colnames(peptideH), "ID_"))
  peptide_ID = dplyr::full_join(peptideL_ID, peptideH_ID, by = "Annotated_Sequence")
  rm(peptideL_ID, peptideH_ID)
  peptideID = data.frame()
  i = 1
  for(identifier in colnames(peptide_ID[,-1])){
    peptideID[i, "Identifier"] = str_replace(identifier, "ID_", "")
    peptideID[i, "Peptide_ID"] = sum(peptide_ID[,which(colnames(peptide_ID) == identifier)] == "High", na.rm = TRUE)
    i = i + 1
  }
  
  Identifications = dplyr::full_join(cells_in_channel_longer, proteinID, by = "Identifier")
  Identifications = dplyr::full_join(Identifications, peptideID, by = "Identifier")
  rm(peptide_ID, Protein_ID, peptideID, proteinID, identifier, i)
  
  # select cells 
  Identifications = Identifications %>%
    dplyr::filter(str_detect(Cell, cell.names))%>% 
    mutate(Carrier_gradient = paste0(Carrier, ", ", Gradient))
  
  return(Identifications)
}

#combine 60 and 90 min
Identifications_90 = get_ID(PDrunName = "ExMloading_90min", cell.names = "HeLa|K562|A549")#|HFL1
Identifications_60 = get_ID(PDrunName = "ExMloading_60min", cell.names = "HeLa|K562|A549")#|HFL1
Identifications = rbind(Identifications_90, Identifications_60)

# plot protein and peptide identifications for different gradient
ggplot(data = Identifications, aes(x = Cell, y = Protein_ID, fill = Carrier_gradient))+
  geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
  xlab("Cell Type")+
  ylab("Identified Protein Groups")+
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = rel(1), color = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1,size = rel(1), color = "black"),
        axis.title.x = element_text(angle = 0, vjust = 0.5,size = rel(1.1), color = "black"), 
        axis.title.y = element_text(angle = 90, hjust = 0.5,size = rel(1.1), color = "black"),
        legend.position = c(0.75, 0.4),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_blank(),
        legend.text = element_text(size = 17)) 
#  + geom_jitter(color="black", size=0.4, alpha=0.9)
ggsave(filename = paste0(figure_save, "Protein ID 60 and 90 min.png"), 
       width = 5, height = 5, dpi = 500)

ggplot(data = Identifications, aes(x = Cell, y = Peptide_ID, fill = Carrier_gradient))+
  #geom_violin(na.rm = TRUE)+
  geom_boxplot(width = 0.3, alpha = 0.5, na.rm = TRUE)+
  xlab("Cell Type")+
  ylab("Identified Peptide Groups")+
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = "top",
        axis.text.x = element_text(angle = 0, hjust = 0.5,size = rel(1.05), color = "black"),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))  + 
  guides(fill=guide_legend(title=""))
#  +geom_jitter(color="black", size=0.4, alpha=0.9)
ggsave(filename = paste0(figure_save, "Pepetide ID2.png"), 
       width = 6, height = 5, dpi = 500)