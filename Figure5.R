
library(readxl)
library(stringr)
library(tidyverse)
library(tidypaleo)
library(cowplot)
library(ggeasy)

marine_cores <- X2023_08_03_shotgun_marine_cores_families_percentage_timeslices
Hauptunterteilung_neu <- read_xlsx("C:/Data/Hauptunterteilung_modified.xlsx", sheet="Sheet1", col_names=TRUE, trim_ws=TRUE)

Hauptunterteilung_neu <- Hauptunterteilung_neu[ ,-which(names(Hauptunterteilung_neu) == "Hauptunterteilung" | names(Hauptunterteilung_neu) == "Ulrike2")]

Hauptunterteilung_neu$Ulrike[Hauptunterteilung_neu$Ulrike == "Southern Hemisphere (< -30 degrees)"] <- "Tropical"
Hauptunterteilung_neu$Ulrike[Hauptunterteilung_neu$Ulrike == "Coastal"] <- "Aquatic"

colnames(Hauptunterteilung_neu)[8] <- "Hauptunterteilung"

Hauptunterteilung_tab <- Hauptunterteilung_neu %>% filter(Hauptunterteilung != "Unspecific")


inputpath_pollen   <- "/Volumes/projects/LegacyPollen_Dataset_timeslices/"

LegacyPollen2_percentages_asia          <- read.csv(paste0(inputpath_pollen,"LegacyPollen2_percentages_asia.csv"), header=TRUE, sep="\t")
LegacyPollen2_percentages_europe        <- read.csv(paste0(inputpath_pollen,"LegacyPollen2_percentages_europe.csv"), header=TRUE, sep="\t")
LegacyPollen2_percentages_north_america <- read.csv(paste0(inputpath_pollen,"LegacyPollen2_percentages_north_america.csv"), header=TRUE, sep="\t")
LegacyPollen2_percentages_south_america <- read.csv(paste0(inputpath_pollen,"LegacyPollen2_percentages_south_america.csv"), header=TRUE, sep="\t")
LegacyPollen2_percentages_africa        <- read.csv(paste0(inputpath_pollen,"LegacyPollen2_percentages_africa.csv"), header=TRUE, sep="\t")
LegacyPollen2_percentages_indopacific   <- read.csv(paste0(inputpath_pollen,"LegacyPollen2_percentages_indopacific.csv"), header=TRUE, sep="\t")



# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

# timeslices: 
timeslices <- seq(0,50000, 2000)



### Preparation Pollen datasets:

{

###### ASIA #####

{
  
  # define help.i
  help.i <- unique(LegacyPollen2_percentages_asia$meanAge)
  
  # create empty output
  output_asia <- data.frame(timeslice=NA, count=NA)
  Meanoutput_asia <- dim(LegacyPollen2_percentages_asia[,-c(1:11)])[2]
  Matrixwithtaxa_asia <- data.frame(matrix(NA, nrow=length(help.i), ncol=Meanoutput_asia +1))
  
  colnames(Matrixwithtaxa_asia) <- c("meanAge", names(LegacyPollen2_percentages_asia[,-c(1:11)]))
  Matrixwithtaxa_asia$meanAge <- help.i
  
  for (i in 1:length(help.i)) {
    
    print(paste0(i,"/",length(help.i)))
    
    sub_asia <- subset(LegacyPollen2_percentages_asia, meanAge == help.i[i])
    
    output_asia[i,"timeslice"] <- help.i[i]
    output_asia[i, "count"] <- nrow(sub_asia)
    Matrixwithtaxa_asia[i, -1] <- unlist(colMeans(sub_asia[,-c(1:11)]))
    
  } 
  
  # -------------------------------------------------------------------------------------------------
  
  asia_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(Matrixwithtaxa_asia, meanAge >= timeslices[i] & meanAge < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-1])))
      
      timeslice_df <-cbind(timeslice = timeslices[i]/1000, y1)
      
      asia_timeslices <- rbind(asia_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

###### EUROPE ##########

{
  
  # define help.i
  help.i <- unique(LegacyPollen2_percentages_europe$meanAge)
  
  # create empty output 
  output_europe <- data.frame(timeslice=NA, count=NA)
  Meanoutput_europe <- dim(LegacyPollen2_percentages_europe[,-c(1:11)])[2]
  Matrixwithtaxa_europe <- data.frame(matrix(NA, nrow=length(help.i), ncol=Meanoutput_europe +1))
  
  colnames(Matrixwithtaxa_europe) <- c("meanAge", names(LegacyPollen2_percentages_europe[,-c(1:11)]))
  Matrixwithtaxa_europe$meanAge <- help.i
  
  for (i in 1:length(help.i)) {
    
    print(paste0(i,"/",length(help.i)))
    
    sub_europe <- subset(LegacyPollen2_percentages_europe, meanAge == help.i[i])
    
    output_europe[i,"timeslice"] <- help.i[i]
    output_europe[i, "count"] <- nrow(sub_europe)
    Matrixwithtaxa_europe[i, -1] <- unlist(colMeans(sub_europe[,-c(1:11)]))
    
  } 
  
  # -------------------------------------------------------------------------------------------------
  
  europe_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(Matrixwithtaxa_europe, meanAge >= timeslices[i] & meanAge < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-1])))
      
      timeslice_df <-cbind(timeslice = timeslices[i]/1000, y1)
      
      europe_timeslices <- rbind(europe_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

###### NORTH AMERICA ########

{
  
  # define help.i
  help.i <- unique(LegacyPollen2_percentages_north_america$meanAge)
  
  # create empty output 
  output_north_america <- data.frame(timeslice=NA, count=NA)
  Meanoutput_north_america <- dim(LegacyPollen2_percentages_north_america[,-c(1:11)])[2]
  Matrixwithtaxa_north_america <- data.frame(matrix(NA, nrow=length(help.i), ncol=Meanoutput_north_america +1))
  
  colnames(Matrixwithtaxa_north_america) <- c("meanAge", names(LegacyPollen2_percentages_north_america[,-c(1:11)]))
  Matrixwithtaxa_north_america$meanAge <- help.i
  
  for (i in 1:length(help.i)) {
    
    print(paste0(i,"/",length(help.i)))
    
    sub_north_america <- subset(LegacyPollen2_percentages_north_america, meanAge == help.i[i])
    
    output_north_america[i,"timeslice"] <- help.i[i]
    output_north_america[i, "count"] <- nrow(sub_north_america)
    Matrixwithtaxa_north_america[i, -1] <- unlist(colMeans(sub_north_america[,-c(1:11)]))
    
  } 
  
  # -------------------------------------------------------------------------------------------------
  
  namerica_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(Matrixwithtaxa_north_america, meanAge >= timeslices[i] & meanAge < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-1])))
      
      timeslice_df <-cbind(timeslice = timeslices[i]/1000, y1)
      
      namerica_timeslices <- rbind(namerica_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

###### SOUTH AMERICA #######

{
  
  # define help.i
  help.i <- unique(LegacyPollen2_percentages_south_america$meanAge)
  
  # create empty output 
  output_south_america <- data.frame(timeslice=NA, count=NA)
  Meanoutput_south_america <- dim(LegacyPollen2_percentages_south_america[,-c(1:11)])[2]
  Matrixwithtaxa_south_america <- data.frame(matrix(NA, nrow=length(help.i), ncol=Meanoutput_south_america +1))
  
  colnames(Matrixwithtaxa_south_america) <- c("meanAge", names(LegacyPollen2_percentages_south_america[,-c(1:11)]))
  Matrixwithtaxa_south_america$meanAge <- help.i
  
  for (i in 1:length(help.i)) {
    
    print(paste0(i,"/",length(help.i)))
    
    sub_south_america <- subset(LegacyPollen2_percentages_south_america, meanAge == help.i[i])
    
    output_south_america[i,"timeslice"] <- help.i[i]
    output_south_america[i, "count"] <- nrow(sub_south_america)
    Matrixwithtaxa_south_america[i, -1] <- unlist(colMeans(sub_south_america[,-c(1:11)]))
    
  } 
  
  # -------------------------------------------------------------------------------------------------
  
  samerica_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(Matrixwithtaxa_south_america, meanAge >= timeslices[i] & meanAge < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-1])))
      
      timeslice_df <-cbind(timeslice = timeslices[i]/1000, y1)
      
      samerica_timeslices <- rbind(samerica_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
   
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

###### AFRICA ######

{
  
  # define help.i
  help.i <- unique(LegacyPollen2_percentages_africa$meanAge)
  
  # create empty output 
  output_africa <- data.frame(timeslice=NA, count=NA)
  Meanoutput_africa <- dim(LegacyPollen2_percentages_africa[,-c(1:11)])[2]
  Matrixwithtaxa_africa <- data.frame(matrix(NA, nrow=length(help.i), ncol=Meanoutput_africa +1))
  
  colnames(Matrixwithtaxa_africa) <- c("meanAge", names(LegacyPollen2_percentages_africa[,-c(1:11)]))
  Matrixwithtaxa_africa$meanAge <- help.i
  
  for (i in 1:length(help.i)) {
    
    print(paste0(i,"/",length(help.i)))
    
    sub_africa <- subset(LegacyPollen2_percentages_africa, meanAge == help.i[i])
    
    output_africa[i,"timeslice"] <- help.i[i]
    output_africa[i, "count"] <- nrow(sub_africa)
    Matrixwithtaxa_africa[i, -1] <- unlist(colMeans(sub_africa[,-c(1:11)]))
    
  } 
  
  # -------------------------------------------------------------------------------------------------
  
  africa_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(Matrixwithtaxa_africa, meanAge >= timeslices[i] & meanAge < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-1])))
      
      timeslice_df <-cbind(timeslice = timeslices[i]/1000, y1)
      
      africa_timeslices <- rbind(africa_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

###### INDOPACIFIC #########

{
  
  # define help.i
  help.i <- unique(LegacyPollen2_percentages_indopacific$meanAge)
  
  # create empty output 
  output_indopacific <- data.frame(timeslice=NA, count=NA)
  Meanoutput_indopacific <- dim(LegacyPollen2_percentages_indopacific[,-c(1:11)])[2]
  Matrixwithtaxa_indopacific <- data.frame(matrix(NA, nrow=length(help.i), ncol=Meanoutput_indopacific +1))
  
  colnames(Matrixwithtaxa_indopacific) <- c("meanAge", names(LegacyPollen2_percentages_indopacific[,-c(1:11)]))
  Matrixwithtaxa_indopacific$meanAge <- help.i
  
  for (i in 1:length(help.i)) {
    
    print(paste0(i,"/",length(help.i)))
    
    sub_indopacific <- subset(LegacyPollen2_percentages_indopacific, meanAge == help.i[i])
    
    output_indopacific[i,"timeslice"] <- help.i[i]
    output_indopacific[i, "count"] <- nrow(sub_indopacific)
    Matrixwithtaxa_indopacific[i, -1] <- unlist(colMeans(sub_indopacific[,-c(1:11)]))
    
  } 
  
  # -------------------------------------------------------------------------------------------------
  
  indopacific_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(Matrixwithtaxa_indopacific, meanAge >= timeslices[i] & meanAge < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-1])))
      
      timeslice_df <-cbind(timeslice = timeslices[i]/1000, y1)
      
      indopacific_timeslices <- rbind(indopacific_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### combine all continental datasets:

pollen_all_continents <- dplyr::bind_rows(asia_timeslices,
                                          europe_timeslices,
                                          namerica_timeslices)

pollen_all_continents[is.na(pollen_all_continents)] <- 0 # replace NA with zero


# select Betulaceae and Salicaceae:
pollen_Alnus_Betula_Populus_Salix <- pollen_all_continents[ ,c(1,
                                                 which(names(pollen_all_continents) == "Alnus"),
                                                 which(names(pollen_all_continents) == "Betula"),
                                                 which(names(pollen_all_continents) == "Populus"),
                                                 which(names(pollen_all_continents) == "Salix"))]



test1 <- pollen_Alnus_Betula_Populus_Salix[pollen_Alnus_Betula_Populus_Salix$timeslice == 0,]


pollen_Betulaceae_Salicaceae_asia <- europe_timeslices %>%
  group_by(timeslice) %>%
  summarise(Alnus_all = sum(Alnus),
            Betula_all = sum(Betula),
            Populus_all = sum(Populus),
            Salix_all = sum(Salix)) %>% 
  mutate(Betulaceae = Alnus_all + Betula_all,
         Salicaceae = Populus_all + Salix_all) %>% 
  mutate(age = timeslice+1)

pollen_Betulaceae_Salicaceae_europe <- asia_timeslices %>%
  group_by(timeslice) %>%
  summarise(Alnus_all = sum(Alnus),
            Betula_all = sum(Betula),
            Populus_all = sum(Populus),
            Salix_all = sum(Salix)) %>% 
  mutate(Betulaceae = Alnus_all + Betula_all,
         Salicaceae = Populus_all + Salix_all) %>% 
  mutate(age = timeslice+1)

pollen_Betulaceae_Salicaceae_namerica <- namerica_timeslices %>%
  group_by(timeslice) %>%
  summarise(Alnus_all = sum(Alnus),
            Betula_all = sum(Betula),
            Populus_all = sum(Populus),
            Salix_all = sum(Salix)) %>% 
  mutate(Betulaceae = Alnus_all + Betula_all,
         Salicaceae = Populus_all + Salix_all) %>% 
  mutate(age = timeslice+1)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

}

# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

### Preparation marine cores:


####

{

ts <- marine_cores_combined[,1]
core_names <- data.frame(timeslice=ts, startage=NA, core_name=NA)

for(i in 1:length(ts)){
  
  if(str_detect(ts[i],"To")){core_names[i,"core_name"] <- "off-Tobago"
  }else if(str_detect(ts[i],"Au")){core_names[i,"core_name"] <- "off-Australia"
  }else if(str_detect(ts[i],"Ant")){core_names[i,"core_name"] <- "Bransfield Strait"
  }else if(str_detect(ts[i],"Ka")){core_names[i,"core_name"] <- "off-Kamchatka"
  }else if(str_detect(ts[i],"BS")){core_names[i,"core_name"] <- "Beringian Sea"
  }else if(str_detect(ts[i],"FS")){core_names[i,"core_name"] <- "Fram Strait"
  }

  core_names[i,"startage"] <- as.numeric(str_extract(ts[i], "[0-9]+"))

}

# -------------------------------------------------------------------------------------------------

}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Plotting:

marine_cores <- marine_cores_combined
level_order <- c("Bransfield Strait", "off-Australia", "Fram Strait", "Beringian Sea", "off-Kamchatka", "off-Tobago")  

colramp <- c("#c8586c","#a3ad62","#bd925a","#C6D9F1","#17365D","#86776F") 

# -------------------------------------------------------------------------------------------------

# Betulaceae Pollen:

{
  
  Betulas_pollen <- ggplot(pollen_Betulaceae_Salicaceae_global, aes(x=age, y=Betulaceae)) +
    geom_bar(stat="identity", width=1.5, fill="grey70") + 
    scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
    easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=6),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  
}

# -------------------------------------------------------------------------------------------------

# Betulaceae:

{
  
  all_cores_Betulaceae <- marine_cores[ ,c(1,which(names(marine_cores) == "Betulaceae"))]
  
  Betulaceae_df <- left_join(core_names, all_cores_Betulaceae)
  
  Betulaceae_df <- Betulaceae_df %>% 
    mutate(age = case_when(core_name == "Bransfield Strait" ~ startage+0.3,
                           core_name == "Fram Strait" ~ startage+0.6,
                           core_name == "off-Kamchatka" ~ startage+0.9,
                           core_name == "off-Australia" ~ startage+1.2,
                           core_name == "Beringian Sea" ~ startage+1.5,
                           core_name == "off-Tobago" ~ startage+1.8)) 
  
  Betulaceae_df <- Betulaceae_df[order(Betulaceae_df$age), ]
  
  
  Betulas <- ggplot(Betulaceae_df, aes(x=age, y=Betulaceae,fill=core_name)) +
    geom_bar(stat="identity", width=0.3) + 
    scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
    scale_y_continuous(limits=c(0,20), breaks=round(seq(min(0), max(20), by =10),1), labels=round(seq(min(0), max(20), by =10),1), expand=c(0,0)) +
    easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=6),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  
  
  Betulas1 <- Betulas + scale_fill_manual(name = "", values = colramp, breaks = level_order) 
  
}

# -------------------------------------------------------------------------------------------------

# Salicaceae Pollen:

{
  
  Salicas_pollen_1 <- ggplot() +
    #geom_bar(stat="identity", width=1.5, fill="grey70") +
    geom_line(pollen_Betulaceae_Salicaceae_asia,aes(x=age, y=Salix_all), color="grey")+
    geom_line(pollen_Betulaceae_Salicaceae_europe,aes(x=age, y=Salix_all), color="darkseagreen")+
    geom_line(pollen_Betulaceae_Salicaceae_namerica,aes(x=age, y=Salix_all), color="forestgreen")+
    geom_line(pollen_Betulaceae_Salicaceae_global, aes(x=age, y=Salix_all), color="black")+
    scale_x_continuous(limits=c(0,30), breaks=round(seq(min(0), max(30), by =2),1), labels=round(seq(min(0), max(30), by =2),1), expand=c(0.025,0)) +
    scale_y_continuous(limits=c(0,8), breaks=round(seq(min(0), max(8), by =2),1), labels=round(seq(min(0), max(8), by =2),1), expand=c(0,0)) +
    #easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=6),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 

}

# -------------------------------------------------------------------------------------------------

# Salicaceae:

{
  
  all_cores_Salicaceae <- marine_cores[ ,c(1,which(names(marine_cores) == "Salicaceae"))]
  
  Salicaceae_df <- left_join(core_names, all_cores_Salicaceae)
  
  Salicaceae_df <- Salicaceae_df %>% 
    mutate(age = case_when(core_name == "Bransfield Strait" ~ startage+0.3,
                           core_name == "Fram Strait" ~ startage+0.6,
                           core_name == "off-Kamchatka" ~ startage+0.9,
                           core_name == "off-Australia" ~ startage+1.2,
                           core_name == "Beringian Sea" ~ startage+1.5,
                           core_name == "off-Tobago" ~ startage+1.8)) 
  
  Salicaceae_df <- Salicaceae_df[order(Salicaceae_df$age), ]
  
  
  Salicas <- ggplot(Salicaceae_df, aes(x=age, y=Salicaceae,fill=core_name)) +
    geom_bar(stat="identity", width=0.3) + 
    scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
    easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=6),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  
  
  Salicas1 <- Salicas + scale_fill_manual(name = "", values = colramp, breaks = level_order) 
  
}

# -------------------------------------------------------------------------------------------------

# Poaceae:

{
  
  all_cores_Poaceae <- marine_cores[ ,c(1,which(names(marine_cores) == "Poaceae"))]
  
  Poaceae_df <- left_join(core_names, all_cores_Poaceae)
  
  Poaceae_df <- Poaceae_df %>% 
    mutate(age = case_when(core_name == "Bransfield Strait" ~ startage+0.3,
                           core_name == "Fram Strait" ~ startage+0.6,
                           core_name == "off-Kamchatka" ~ startage+0.9,
                           core_name == "off-Australia" ~ startage+1.2,
                           core_name == "Beringian Sea" ~ startage+1.5,
                           core_name == "off-Tobago" ~ startage+1.8)) 
  
  Poaceae_df <- Poaceae_df[order(Poaceae_df$age), ]
  
  ### highest value (76.2%) is manually set to 35%:
  Poaceae_df$Poaceae[Poaceae_df$Poaceae > 35] <- 35
  
  
  Poas <- ggplot(Poaceae_df, aes(x=age, y=Poaceae,fill=core_name)) +
    geom_bar(stat="identity", width=0.3) + 
    scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
    scale_y_continuous(limits=c(0,40), breaks=round(seq(min(0), max(40), by =20),1), labels=round(seq(min(0), max(40), by =20),1), expand=c(0,0)) +
    easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=6),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  
  
  Poas1 <- Poas + scale_fill_manual(name = "", values = colramp, breaks = level_order) 
  
}

# -------------------------------------------------------------------------------------------------

# Mosses/Ferns:

{
  
  moss_fern_taxa <- subset(Hauptunterteilung_tab, Hauptunterteilung == "Moss/Fern")
  
  moss_fern_marine_cores <- marine_cores[ ,c(which(names(marine_cores) %in% moss_fern_taxa$family))]
  moss_fern_marine_cores <- moss_fern_marine_cores %>% 
    mutate(MossFern = rowSums(moss_fern_marine_cores))
  
  all_cores_MossFern <- data.frame(timeslice=marine_cores[ ,1], MossFern=moss_fern_marine_cores[ ,which(names(moss_fern_marine_cores) == "MossFern")])
  
  MossFern_df <- left_join(core_names, all_cores_MossFern)
  
  MossFern_df <- MossFern_df %>% 
    mutate(age = case_when(core_name == "Bransfield Strait" ~ startage+0.3,
                           core_name == "Fram Strait" ~ startage+0.6,
                           core_name == "off-Kamchatka" ~ startage+0.9,
                           core_name == "off-Australia" ~ startage+1.2,
                           core_name == "Beringian Sea" ~ startage+1.5,
                           core_name == "off-Tobago" ~ startage+1.8)) 
  
  MossFern_df <- MossFern_df[order(MossFern_df$age), ]
  
  
  MossFerns <- ggplot(MossFern_df, aes(x=age, y=MossFern,fill=core_name)) +
    geom_bar(stat="identity", width=0.3) + 
    scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
    easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=6),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  
  
  MossFerns1 <- MossFerns + scale_fill_manual(name = "", values = colramp, breaks = level_order) 
  
}

# -------------------------------------------------------------------------------------------------

# Aquatics:

{
  
  aquatic_taxa <- subset(Hauptunterteilung_tab, Hauptunterteilung == "Aquatic")
  
  aquatic_marine_cores <- marine_cores[ ,c(which(names(marine_cores) %in% aquatic_taxa$family))]
  aquatic_marine_cores <- aquatic_marine_cores %>% 
    mutate(Aquatic = rowSums(aquatic_marine_cores))
  
  all_cores_Aquatic <- data.frame(timeslice=marine_cores[ ,1], Aquatic=aquatic_marine_cores[ ,which(names(aquatic_marine_cores) == "Aquatic")])
  
  Aquatic_df <- left_join(core_names, all_cores_Aquatic)
  
  Aquatic_df <- Aquatic_df %>% 
    mutate(age = case_when(core_name == "Bransfield Strait" ~ startage+0.3,
                           core_name == "Fram Strait" ~ startage+0.6,
                           core_name == "off-Kamchatka" ~ startage+0.9,
                           core_name == "off-Australia" ~ startage+1.2,
                           core_name == "Beringian Sea" ~ startage+1.5,
                           core_name == "off-Tobago" ~ startage+1.8)) 
  
  Aquatic_df <- Aquatic_df[order(Aquatic_df$age), ]
  
  ### highest value (50.0%) is manually set to 35%:
  Aquatic_df$Aquatic[Aquatic_df$Aquatic > 35] <- 35
  
  
  # Aquatics <- ggplot(Aquatic_df, aes(x=age, y=Aquatic,fill=core_name)) +
  #   geom_bar(stat="identity", width=0.3) + 
  #   scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
  #   easy_remove_x_axis() + 
  #   theme(plot.title=element_text(size=10),
  #         plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
  #         axis.title=element_text(size=6),
  #         axis.text.y=element_text(size=10)) + 
  #   theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
  #   theme(panel.background = element_blank())+ 
  #   theme(panel.grid.major=element_blank(),
  #         panel.grid.minor=element_blank()) +
  #   theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  
  
  
  Aquatics <- ggplot(Aquatic_df, aes(x=age, y=Aquatic,fill=core_name)) +
    geom_bar(stat="identity", width=0.3) + 
    scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
    scale_y_continuous(limits=c(0,35), breaks=round(seq(min(0), max(40), by =20),1), labels=round(seq(min(0), max(40), by =20),1), expand=c(0,0)) +
    #easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=6),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.x = element_line(color="grey50", linewidth = 0.2), 
          axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  
  
  Aquatics1 <- Aquatics + scale_fill_manual(name = "", values = colramp, breaks = level_order) 
  
  
  
  
  # 
  # 
  Aquatics_temp <- ggplot(Aquatic_df, aes(x=age, y=Aquatic,fill=core_name)) +
    geom_bar(stat="identity", width=0.3) + 
    scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
    #  scale_y_continuous(limits=c(0,32), breaks=round(seq(min(0), max(30), by =10),1), labels=round(seq(min(0), max(30), by =10),1), expand=c(0,0)) +
    easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=6),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "bottom") 



  Aquatics_temp1 <- Aquatics_temp + scale_fill_manual(name = "", values = colramp, breaks = level_order)


  
}

# -------------------------------------------------------------------------------------------------

##### plot all #######
plot_grid(Betulas_pollen,
          Betulas1,
          Salicas_pollen,
          Salicas1,
          Poas1,
          MossFerns1,
          Aquatics1,
          Aquatics_temp1,
          ncol=1,nrow=8)


ggsave("C:/Data/2023-05-31_stratigraphic_plot_marine_cores_taxa.png", width = 18, height = 18, units = "cm", dpi = 300)


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------




tobago     <- marine_cores[str_detect(marine_cores$timeslice, "To"), ]
australia  <- marine_cores[str_detect(marine_cores$timeslice, "Au"), ]
antarctica <- marine_cores[str_detect(marine_cores$timeslice, "Ant"), ]
kamchatka  <- marine_cores[str_detect(marine_cores$timeslice, "Ka"), ]
beringia   <- marine_cores[str_detect(marine_cores$timeslice, "BS"), ]
fram       <- marine_cores[str_detect(marine_cores$timeslice, "FS"), ]
  



tobago_means     <- colMeans(tobago[,-1])
australia_means  <- colMeans(australia[,-1])
antarctica_means <- colMeans(antarctica[,-1])
kamchatka_means  <- colMeans(kamchatka[,-1])
beringia_means   <- colMeans(beringia[,-1])
fram_means       <- colMeans(fram[,-1])

marine_cores_ts_means <- rbind(tobago_means,
                               australia_means,
                               antarctica_means,
                               kamchatka_means,
                               beringia_means,
                               fram_means)


marine_cores_overall_mean <- colMeans(marine_cores_ts_means)


library(treemapify)

marine_cores_overall_mean_df <- as.data.frame(marine_cores_overall_mean)
marine_cores_overall_mean_df$families <- row.names(marine_cores_overall_mean_df) 
marine_cores_overall_mean_df$category <- 1
colnames(marine_cores_overall_mean_df)[1] <- "perc_means"




marine_cores_treeplot <-ggplot(marine_cores_overall_mean_df, aes(area=perc_means, label=families, fill=category)) +  ##fill=Hauptunterteilung, , subgroup = Hauptunterteilung
  # labs(title="Tobago Holocene") +
  geom_treemap()+
  geom_treemap_text(colour = "white")+
  theme(legend.position = "none")+
#  scale_fill_manual(name = "", values= "#089099") + 
  #scale_fill_manual(name = "", values= c("#46aea0","#089099","#00718b","#045275"), breaks = level_order, labels=labs) + 
  # scale_fill_viridis_d(option = "viridis", begin=0,end=0.6) +
  theme(plot.title=element_text(size=10))

marine_cores_treeplot


ggsave("C:/Data/2023-05-31_treemaps_marine_cores_overall_means.jpg", width = 30, height = 30, units = "cm", dpi = 300)







