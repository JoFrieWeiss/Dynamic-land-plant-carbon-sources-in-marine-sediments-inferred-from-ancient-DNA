
library(readxl)
library(data.table)
library(stringr)
library(tidyverse)
library(tidypaleo)
library(cowplot)
library(ggeasy)

inputpath <- "/Volumes/projects/user/shotgun_data/05.2_shotgun_data_nt05/" 

# -------------------------------------------------------------------------------------------------

Hauptunterteilung_neu <- read_xlsx("C:/Data/Hauptunterteilung_modified.xlsx", sheet="Sheet1", col_names=TRUE, trim_ws=TRUE)
Hauptunterteilung_neu <- Hauptunterteilung_neu[ ,-which(names(Hauptunterteilung_neu) == "Hauptunterteilung" | names(Hauptunterteilung_neu) == "Ulrike2")]
Hauptunterteilung_neu$Ulrike[Hauptunterteilung_neu$Ulrike == "Southern Hemisphere (< -30 degrees)"] <- "Tropical"
Hauptunterteilung_neu$Ulrike[Hauptunterteilung_neu$Ulrike == "Coastal"] <- "Aquatic"
colnames(Hauptunterteilung_neu)[8] <- "Hauptunterteilung"
Hauptunterteilung_tab <- Hauptunterteilung_neu %>% filter(Hauptunterteilung != "Unspecific")

# aquatic taxa that should be excluded:
coastal <- as.vector(unlist(Hauptunterteilung_tab[Hauptunterteilung_tab$Hauptunterteilung == "Aquatic","family"]))

# classes that should be excluded:
target <- c("Chlorokybophyceae", "Klebsormidiophyceae", "Mesostigmatophyceae", "Charophyceae", "Coleochaetophyceae", "Zygnemophyceae")

# timeslices: 
timeslices <- seq(0,50000, 2000)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### marine multicore:

{
  
  shotgun_multi <- fread(paste0(inputpath,"marine_NextSeq_BHV/marineNextSeqBHV_nt0.5_Kraken2.txt")) # the readin data is a file of combined krona reports 
  names(shotgun_multi) <- c("samples", "Percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")
  
  # load lineage using viktor's python script 
  lineage_multi <- fread(paste0(inputpath,"marine_NextSeq_BHV/marineNextSeqBHV_nt0.5_Kraken2_lineageDB.csv"))  
  
  # add lineage to the shotgun data
  shotgun.lineage_multi <- plyr::join(shotgun_multi, lineage_multi, by = "taxID")
  
  # load metadata containing sample name, age and depth
  metadata.readin_multi <- openxlsx::read.xlsx(xlsxFile = paste0(inputpath,"marine_NextSeq_BHV/JoW008L_nt05_shotgun_sample_name_change.xlsx"))
  unique(metadata.readin_multi$core) # check how many different cores in the sequencing run
  
  # -------------------------------------------------------------------------------------------------
  
  ### Tobago Bay:
  
  {
    
    metadata_tobago <- filter(metadata.readin_multi, core == "M78-1-235") # filter to core M78-1-235 #Tobago Bay
    
    # remove whitespace at the end of the string:
    metadata_tobago$samples <- str_trim(metadata_tobago$samples)
    
    # add metadata to the shotgun data
    shotgun.join_tobago <- plyr::join(shotgun.lineage_multi, metadata_tobago, by = "samples")
    
    shotgun.join_tobago_samples <- shotgun.join_tobago %>% filter(read_fraction == "merged" & type == "sample")
    
    # -------------------------------------------------------------------------------------------------
    
    shotgun.join_tobago_eukaryota <- shotgun.join_tobago_samples %>% filter(superkingdom == "Eukaryota")
    
    shotgun.join_tobago_eukaryota_total <- shotgun.join_tobago_eukaryota %>% 
      group_by(age) %>% 
      summarise(eukaryota_total = sum(CladeCount))
    
    shotgun.join_tobago_streptophyta <- shotgun.join_tobago_eukaryota %>% filter(phylum == "Streptophyta")
    
    shotgun.join_tobago_streptophyta_total <- shotgun.join_tobago_streptophyta %>% 
      group_by(age) %>% 
      summarise(streptophyta_total = sum(CladeCount))
    
    shotgun.join_tobago_other <- shotgun.join_tobago_eukaryota %>% filter(!phylum == "Streptophyta")
    
    shotgun.join_tobago_other_total <- shotgun.join_tobago_other %>% 
      group_by(age) %>% 
      summarise(other_total = sum(CladeCount)) 
    
    shotgun.join_tobago_embryophyta  <- shotgun.join_tobago_streptophyta %>% filter(!class %in% target)
    
    shotgun.join_tobago_embryophyta_terra <- shotgun.join_tobago_embryophyta %>% filter(!family %in% coastal)
    shotgun.join_tobago_embryophyta_aqua  <- shotgun.join_tobago_embryophyta %>% filter(family %in% coastal)
    
    shotgun.join_tobago_embryophyta_terra_total <- shotgun.join_tobago_embryophyta_terra %>% 
      group_by(age) %>% 
      summarise(embryophyta_terrestric_total = sum(CladeCount))
    
    shotgun.join_tobago_embryophyta_aqua_total <- shotgun.join_tobago_embryophyta_aqua %>% 
      group_by(age) %>% 
      summarise(embryophyta_aquatic_total = sum(CladeCount))
    
    # -------------------------------------------------------------------------------------------------
    
    tobago_eukaryota_embryophyta_df1 <- left_join(shotgun.join_tobago_eukaryota_total, shotgun.join_tobago_streptophyta_total)
    tobago_eukaryota_embryophyta_df2 <- left_join(tobago_eukaryota_embryophyta_df1, shotgun.join_tobago_other_total)
    tobago_eukaryota_embryophyta_df3 <- left_join(shotgun.join_tobago_embryophyta_terra_total, shotgun.join_tobago_embryophyta_aqua_total)
    
    tobago_eukaryota_embryophyta_df  <- left_join(tobago_eukaryota_embryophyta_df2, tobago_eukaryota_embryophyta_df3)
    
    tobago_eukaryota_embryophyta_df[is.na(tobago_eukaryota_embryophyta_df)] <- 0
    
    tobago_eukaryota_embryophyta_df <- tobago_eukaryota_embryophyta_df %>% 
      group_by(age) %>% 
      mutate(embryophyta_total            = (embryophyta_terrestric_total + embryophyta_aquatic_total),
             ratio_eukaryota_streptophyta = (streptophyta_total/eukaryota_total)*100,
             ratio_eukaryota_other        = (other_total/eukaryota_total)*100,
             ratio_eukaryota_terrestric   = (embryophyta_terrestric_total/eukaryota_total)*100,
             ratio_eukaryota_aquatic      = (embryophyta_aquatic_total/eukaryota_total)*100,
             ratio_eukaryota_embryophyta  = (embryophyta_total/eukaryota_total)*100) 
    
    tobago_eukaryota_embryophyta_df$core_name <- "off-Tobago"
    
    # -------------------------------------------------------------------------------------------------
    
    tobago_eukaryota_embryophyta_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(tobago_eukaryota_embryophyta_df, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-c(1,13)])))
        
        timeslice_df <-unique(cbind(timeslice=(timeslices[i]/1000),y1,core_name=x1$core_name))
        
        tobago_eukaryota_embryophyta_timeslices <- rbind(tobago_eukaryota_embryophyta_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
     
    # write.csv(tobago_eukaryota_embryophyta_timeslices, file="C:/Data/Marine_cores/2023-06-01_eukaryota_embryophyta_tobago.csv", row.names=FALSE)
    
    # -------------------------------------------------------------------------------------------------
    
    # for PCA - terrestrial:
    
    {
      
      shotgun.join_tobago_embryophyta_terra_families <- shotgun.join_tobago_embryophyta_terra %>% filter(Rank == "F")
      
      tobago_embryophyta_terra_select <- select(shotgun.join_tobago_embryophyta_terra_families, c("age", "Name", "CladeCount"))
      
      tobago_embryophyta_terra_taxatab <- pivot_wider(tobago_embryophyta_terra_select, names_from = Name, values_from = CladeCount)
      tobago_embryophyta_terra_taxatab[is.na(tobago_embryophyta_terra_taxatab)] <- 0 # replace NA with zero
      
      tobago_embryophyta_terra_counts <- as.data.frame(tobago_embryophyta_terra_taxatab[order(tobago_embryophyta_terra_taxatab$age), ])
      
      tobago_embryophyta_terra_perc <- data.frame(age = tobago_embryophyta_terra_counts[ ,1], 
                                                  (tobago_embryophyta_terra_counts[ ,-1] / rowSums(tobago_embryophyta_terra_counts[ ,-1]) * 100))
      
      # -------------------------------------------------------------------------------------------------
      
      tobago_embryophyta_terra_perc_timeslices <- NULL
      
      
      for(i in 1:(length(timeslices)-1)){
        
        x1 <- subset(tobago_embryophyta_terra_perc, age >= timeslices[i] & age < timeslices[i+1])
        
        if(dim(x1)[1] > 0){
          
          y1 <- as.data.frame(t(colMeans(x1[,-1])))
          
          timeslice_df <-unique(cbind(timeslice = paste0("To",timeslices[i]/1000),y1))
          
          tobago_embryophyta_terra_perc_timeslices <- rbind(tobago_embryophyta_terra_perc_timeslices, timeslice_df)
          
        }
      }
      
      # -------------------------------------------------------------------------------------------------
      
    }
    
    # -------------------------------------------------------------------------------------------------
    
    # for PCA - aquatic:
    
    {
      
      shotgun.join_tobago_embryophyta_aqua_families <- shotgun.join_tobago_embryophyta_aqua %>% filter(Rank == "F")
      
      tobago_embryophyta_aqua_select <- select(shotgun.join_tobago_embryophyta_aqua_families, c("age", "Name", "CladeCount"))
      
      tobago_embryophyta_aqua_taxatab <- pivot_wider(tobago_embryophyta_aqua_select, names_from = Name, values_from = CladeCount)
      tobago_embryophyta_aqua_taxatab[is.na(tobago_embryophyta_aqua_taxatab)] <- 0 # replace NA with zero
      
      tobago_embryophyta_aqua_counts <- as.data.frame(tobago_embryophyta_aqua_taxatab[order(tobago_embryophyta_aqua_taxatab$age), ])
      
      tobago_embryophyta_aqua_perc <- data.frame(age = tobago_embryophyta_aqua_counts[ ,1], 
                                                 (tobago_embryophyta_aqua_counts[ ,-1] / rowSums(tobago_embryophyta_aqua_counts[ ,-1]) * 100))
      
      # -------------------------------------------------------------------------------------------------
      
      tobago_embryophyta_aqua_perc_timeslices <- NULL
      
      
      for(i in 1:(length(timeslices)-1)){
        
        x1 <- subset(tobago_embryophyta_aqua_perc, age >= timeslices[i] & age < timeslices[i+1])
        
        if(dim(x1)[1] > 0){
          
          y1 <- as.data.frame(t(colMeans(x1[,-1])))
          
          timeslice_df <-unique(cbind(timeslice = paste0("To",timeslices[i]/1000),y1))
          
          tobago_embryophyta_aqua_perc_timeslices <- rbind(tobago_embryophyta_aqua_perc_timeslices, timeslice_df)
          
        }
      }
      
      # -------------------------------------------------------------------------------------------------
      
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------
  
  ### Australia:
  
  {
    
    metadata_australia <- filter(metadata.readin_multi, core == "MD03-2614G") # filter to site MD03-2614G
    
    # remove whitespace at the end of the string:
    metadata_australia$samples <- str_trim(metadata_australia$samples)
    
    # add metadata to the shotgun data
    shotgun.join_australia <- plyr::join(shotgun.lineage_multi, metadata_australia, by = "samples")
    
    shotgun.join_australia_samples <- shotgun.join_australia %>% filter(read_fraction == "merged" & type == "sample")
    
    # -------------------------------------------------------------------------------------------------
    
    shotgun.join_australia_eukaryota <- shotgun.join_australia_samples %>% filter(superkingdom == "Eukaryota")
    
    shotgun.join_australia_eukaryota_total <- shotgun.join_australia_eukaryota %>% 
      group_by(age) %>% 
      summarise(eukaryota_total = sum(CladeCount))
    
    shotgun.join_australia_streptophyta <- shotgun.join_australia_eukaryota %>% filter(phylum == "Streptophyta")
    
    shotgun.join_australia_streptophyta_total <- shotgun.join_australia_streptophyta %>% 
      group_by(age) %>% 
      summarise(streptophyta_total = sum(CladeCount))
    
    shotgun.join_australia_other <- shotgun.join_australia_eukaryota %>% filter(!phylum == "Streptophyta")
    
    shotgun.join_australia_other_total <- shotgun.join_australia_other %>% 
      group_by(age) %>% 
      summarise(other_total = sum(CladeCount))
    
    shotgun.join_australia_embryophyta  <- shotgun.join_australia_streptophyta %>% filter(!class %in% target)
    
    shotgun.join_australia_embryophyta_terra <- shotgun.join_australia_embryophyta %>% filter(!family %in% coastal)
    shotgun.join_australia_embryophyta_aqua  <- shotgun.join_australia_embryophyta %>% filter(family %in% coastal)
    
    shotgun.join_australia_embryophyta_terra_total <- shotgun.join_australia_embryophyta_terra %>% 
      group_by(age) %>% 
      summarise(embryophyta_terrestric_total = sum(CladeCount))
    
    shotgun.join_australia_embryophyta_aqua_total <- shotgun.join_australia_embryophyta_aqua %>% 
      group_by(age) %>% 
      summarise(embryophyta_aquatic_total = sum(CladeCount))
    
    # -------------------------------------------------------------------------------------------------
    
    australia_eukaryota_embryophyta_df1 <- left_join(shotgun.join_australia_eukaryota_total, shotgun.join_australia_streptophyta_total)
    australia_eukaryota_embryophyta_df2 <- left_join(australia_eukaryota_embryophyta_df1, shotgun.join_australia_other_total)
    australia_eukaryota_embryophyta_df3 <- left_join(shotgun.join_australia_embryophyta_terra_total, shotgun.join_australia_embryophyta_aqua_total)
    
    australia_eukaryota_embryophyta_df  <- left_join(australia_eukaryota_embryophyta_df2, australia_eukaryota_embryophyta_df3)
    
    australia_eukaryota_embryophyta_df[is.na(australia_eukaryota_embryophyta_df)] <- 0
    
    australia_eukaryota_embryophyta_df <- australia_eukaryota_embryophyta_df %>% 
      group_by(age) %>% 
      mutate(embryophyta_total            = (embryophyta_terrestric_total + embryophyta_aquatic_total),
             ratio_eukaryota_streptophyta = (streptophyta_total/eukaryota_total)*100,
             ratio_eukaryota_other        = (other_total/eukaryota_total)*100,
             ratio_eukaryota_terrestric   = (embryophyta_terrestric_total/eukaryota_total)*100,
             ratio_eukaryota_aquatic      = (embryophyta_aquatic_total/eukaryota_total)*100,
             ratio_eukaryota_embryophyta  = (embryophyta_total/eukaryota_total)*100) 
    
    australia_eukaryota_embryophyta_df$core_name <- "off-Australia"
    
    # -------------------------------------------------------------------------------------------------
    
    australia_eukaryota_embryophyta_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(australia_eukaryota_embryophyta_df, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-c(1,13)])))
        
        timeslice_df <-unique(cbind(timeslice=(timeslices[i]/1000),y1,core_name=x1$core_name))
        
        australia_eukaryota_embryophyta_timeslices <- rbind(australia_eukaryota_embryophyta_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
    # write.csv(australia_eukaryota_embryophyta_timeslices, file="C:/Data/Marine_cores/2023-06-01_eukaryota_embryophyta_australia.csv", row.names=FALSE)
    
    # -------------------------------------------------------------------------------------------------
    
    # for PCA - terrestrial:
    
    {
      
      shotgun.join_australia_embryophyta_terra_families <- shotgun.join_australia_embryophyta_terra %>% filter(Rank == "F")
      
      australia_embryophyta_terra_select <- select(shotgun.join_australia_embryophyta_terra_families, c("age", "Name", "CladeCount"))
      
      australia_embryophyta_terra_taxatab <- pivot_wider(australia_embryophyta_terra_select, names_from = Name, values_from = CladeCount)
      australia_embryophyta_terra_taxatab[is.na(australia_embryophyta_terra_taxatab)] <- 0 # replace NA with zero
      
      australia_embryophyta_terra_counts <- as.data.frame(australia_embryophyta_terra_taxatab[order(australia_embryophyta_terra_taxatab$age), ])
      
      australia_embryophyta_terra_perc <- data.frame(age = australia_embryophyta_terra_counts[ ,1], 
                                                     (australia_embryophyta_terra_counts[ ,-1] / rowSums(australia_embryophyta_terra_counts[ ,-1]) * 100))
      
      # -------------------------------------------------------------------------------------------------
      
      australia_embryophyta_terra_perc_timeslices <- NULL
      
      
      for(i in 1:(length(timeslices)-1)){
        
        x1 <- subset(australia_embryophyta_terra_perc, age >= timeslices[i] & age < timeslices[i+1])
        
        if(dim(x1)[1] > 0){
          
          y1 <- as.data.frame(t(colMeans(x1[,-1])))
          
          timeslice_df <-unique(cbind(timeslice = paste0("Au",timeslices[i]/1000),y1))
          
          australia_embryophyta_terra_perc_timeslices <- rbind(australia_embryophyta_terra_perc_timeslices, timeslice_df)
          
        }
      }
      
      # -------------------------------------------------------------------------------------------------
      
    }
    
    # -------------------------------------------------------------------------------------------------
    
    # for PCA - aquatic:
    
    {
      
      shotgun.join_australia_embryophyta_aqua_families <- shotgun.join_australia_embryophyta_aqua %>% filter(Rank == "F")
      
      australia_embryophyta_aqua_select <- select(shotgun.join_australia_embryophyta_aqua_families, c("age", "Name", "CladeCount"))
      
      australia_embryophyta_aqua_taxatab <- pivot_wider(australia_embryophyta_aqua_select, names_from = Name, values_from = CladeCount)
      australia_embryophyta_aqua_taxatab[is.na(australia_embryophyta_aqua_taxatab)] <- 0 # replace NA with zero
      
      australia_embryophyta_aqua_counts <- as.data.frame(australia_embryophyta_aqua_taxatab[order(australia_embryophyta_aqua_taxatab$age), ])
      
      australia_embryophyta_aqua_perc <- data.frame(age = australia_embryophyta_aqua_counts[ ,1], 
                                                 (australia_embryophyta_aqua_counts[ ,-1] / rowSums(australia_embryophyta_aqua_counts[ ,-1]) * 100))
      
      # -------------------------------------------------------------------------------------------------
      
      australia_embryophyta_aqua_perc_timeslices <- NULL
      
      
      for(i in 1:(length(timeslices)-1)){
        
        x1 <- subset(australia_embryophyta_aqua_perc, age >= timeslices[i] & age < timeslices[i+1])
        
        if(dim(x1)[1] > 0){
          
          y1 <- as.data.frame(t(colMeans(x1[,-1])))
          
          timeslice_df <-unique(cbind(timeslice = paste0("Au",timeslices[i]/1000),y1))
          
          australia_embryophyta_aqua_perc_timeslices <- rbind(australia_embryophyta_aqua_perc_timeslices, timeslice_df)
          
        }
      }
      
      # -------------------------------------------------------------------------------------------------
      
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------
  
  ### Fram Strait:   !!!! incomplete; only new sequencing
  
  {
    
    metadata_fram <- filter(metadata.readin_multi, core == "MSM05-5-712-2") # filter to core MSM05-5-712-2
    
    # remove whitespace at the end of the string:
    metadata_fram$samples <- str_trim(metadata_fram$samples)
    
    # add metadata to the shotgun data
    shotgun.join_fram <- plyr::join(shotgun.lineage_multi, metadata_fram, by = "samples")
    
    shotgun.join_fram_samples <- shotgun.join_fram %>% filter(read_fraction == "merged" & type == "sample")
    
  }
  
  # -------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Antarctica:

{
  
  shotgun_antarctica <- fread(paste0(inputpath,"PS97_NextSeqBHV/PS97_nt0.5_Kraken2.txt")) # the readin data is a file of combined krona reports 
  names(shotgun_antarctica) <- c("samples", "Percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")
  
  # load lineage using viktor's python script 
  lineage_antarctica <- fread(paste0(inputpath,"PS97_NextSeqBHV/PS97_nt0.5_Kraken2_lineageDB.csv"))  
  
  # add lineage to the shotgun data
  shotgun.lineage_antarctica <- plyr::join(shotgun_antarctica, lineage_antarctica, by = "taxID")
  
  # load metadata containing sample name, age and depth
  metadata.readin_antarctica <- openxlsx::read.xlsx(xlsxFile = paste0(inputpath,"PS97_NextSeqBHV/PS97-72-01_shotgun_sample_name_nt0.5.xlsx"))
  unique(metadata.readin_antarctica$core) # check how many different cores in the sequencing run
  
  # -------------------------------------------------------------------------------------------------
  
  metadata_antarctica <- metadata.readin_antarctica
  
  # remove whitespace at the end of the string:
  metadata_antarctica$samples <- str_trim(metadata_antarctica$samples)
  
  # add metadata to the shotgun data
  shotgun.join_antarctica <- plyr::join(shotgun.lineage_antarctica, metadata_antarctica, by = "samples")
  
  shotgun.join_antarctica_samples <- shotgun.join_antarctica %>% filter(read_fraction == "merged" & type == "sample")
  
  # -------------------------------------------------------------------------------------------------
  
  shotgun.join_antarctica_eukaryota <- shotgun.join_antarctica_samples %>% filter(superkingdom == "Eukaryota")
  
  shotgun.join_antarctica_eukaryota_total <- shotgun.join_antarctica_eukaryota %>% 
    group_by(age) %>% 
    summarise(eukaryota_total = sum(CladeCount))
  
  shotgun.join_antarctica_streptophyta <- shotgun.join_antarctica_eukaryota %>% filter(phylum == "Streptophyta")
  
  shotgun.join_antarctica_streptophyta_total <- shotgun.join_antarctica_streptophyta %>% 
    group_by(age) %>% 
    summarise(streptophyta_total = sum(CladeCount))
  
  shotgun.join_antarctica_other <- shotgun.join_antarctica_eukaryota %>% filter(!phylum == "Streptophyta")
  
  shotgun.join_antarctica_other_total <- shotgun.join_antarctica_other %>% 
    group_by(age) %>% 
    summarise(other_total = sum(CladeCount))
  
  shotgun.join_antarctica_embryophyta <- shotgun.join_antarctica_streptophyta %>% filter(!class %in% target)
  
  shotgun.join_antarctica_embryophyta_terra <- shotgun.join_antarctica_embryophyta %>% filter(!family %in% coastal)
  shotgun.join_antarctica_embryophyta_aqua  <- shotgun.join_antarctica_embryophyta %>% filter(family %in% coastal)
  
  shotgun.join_antarctica_embryophyta_terra_total <- shotgun.join_antarctica_embryophyta_terra %>% 
    group_by(age) %>% 
    summarise(embryophyta_terrestric_total = sum(CladeCount))
  
  shotgun.join_antarctica_embryophyta_aqua_total <- shotgun.join_antarctica_embryophyta_aqua %>% 
    group_by(age) %>% 
    summarise(embryophyta_aquatic_total = sum(CladeCount))
 
  # -------------------------------------------------------------------------------------------------
  
  antarctica_eukaryota_embryophyta_df1 <- left_join(shotgun.join_antarctica_eukaryota_total, shotgun.join_antarctica_streptophyta_total)
  antarctica_eukaryota_embryophyta_df2 <- left_join(antarctica_eukaryota_embryophyta_df1, shotgun.join_antarctica_other_total)
  antarctica_eukaryota_embryophyta_df3 <- left_join(shotgun.join_antarctica_embryophyta_terra_total, shotgun.join_antarctica_embryophyta_aqua_total)
  
  antarctica_eukaryota_embryophyta_df  <- left_join(antarctica_eukaryota_embryophyta_df2, antarctica_eukaryota_embryophyta_df3)
  
  antarctica_eukaryota_embryophyta_df[is.na(antarctica_eukaryota_embryophyta_df)] <- 0
  
  antarctica_eukaryota_embryophyta_df <- antarctica_eukaryota_embryophyta_df %>% 
    group_by(age) %>% 
    mutate(embryophyta_total            = (embryophyta_terrestric_total + embryophyta_aquatic_total),
           ratio_eukaryota_streptophyta = (streptophyta_total/eukaryota_total)*100,
           ratio_eukaryota_other        = (other_total/eukaryota_total)*100,
           ratio_eukaryota_terrestric   = (embryophyta_terrestric_total/eukaryota_total)*100,
           ratio_eukaryota_aquatic      = (embryophyta_aquatic_total/eukaryota_total)*100,
           ratio_eukaryota_embryophyta  = (embryophyta_total/eukaryota_total)*100) 
  
  antarctica_eukaryota_embryophyta_df$core_name <- "Bransfield Strait"
 
  # -------------------------------------------------------------------------------------------------
  
  antarctica_eukaryota_embryophyta_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(antarctica_eukaryota_embryophyta_df, age >= timeslices[i] & age < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-c(1,13)])))
      
      timeslice_df <-unique(cbind(timeslice=(timeslices[i]/1000),y1,core_name=x1$core_name))
      
      antarctica_eukaryota_embryophyta_timeslices <- rbind(antarctica_eukaryota_embryophyta_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # write.csv(antarctica_eukaryota_embryophyta_timeslices, file="C:/Data/Marine_cores/2023-06-01_eukaryota_embryophyta_antarctica.csv", row.names=FALSE)
  
  # -------------------------------------------------------------------------------------------------
  
  # for PCA - terrestrial:
  
  {
    
    shotgun.join_antarctica_embryophyta_terra_families <- shotgun.join_antarctica_embryophyta_terra %>% filter(Rank == "F")
    
    antarctica_embryophyta_terra_select <- select(shotgun.join_antarctica_embryophyta_terra_families, c("age", "Name", "CladeCount"))
    
    antarctica_embryophyta_terra_taxatab <- pivot_wider(antarctica_embryophyta_terra_select, names_from = Name, values_from = CladeCount)
    antarctica_embryophyta_terra_taxatab[is.na(antarctica_embryophyta_terra_taxatab)] <- 0 # replace NA with zero
    
    antarctica_embryophyta_terra_counts <- as.data.frame(antarctica_embryophyta_terra_taxatab[order(antarctica_embryophyta_terra_taxatab$age), ])
    
    antarctica_embryophyta_terra_perc <- data.frame(age = antarctica_embryophyta_terra_counts[ ,1], 
                                                (antarctica_embryophyta_terra_counts[ ,-1] / rowSums(antarctica_embryophyta_terra_counts[ ,-1]) * 100))
    
    # -------------------------------------------------------------------------------------------------
    
    antarctica_embryophyta_terra_perc_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(antarctica_embryophyta_terra_perc, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-1])))
        
        timeslice_df <-unique(cbind(timeslice = paste0("Ant",timeslices[i]/1000),y1))
        
        antarctica_embryophyta_terra_perc_timeslices <- rbind(antarctica_embryophyta_terra_perc_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # for PCA - aquatic:
  
  {
    
    shotgun.join_antarctica_embryophyta_aqua_families <- shotgun.join_antarctica_embryophyta_aqua %>% filter(Rank == "F")
    
    antarctica_embryophyta_aqua_select <- select(shotgun.join_antarctica_embryophyta_aqua_families, c("age", "Name", "CladeCount"))
    
    antarctica_embryophyta_aqua_taxatab <- pivot_wider(antarctica_embryophyta_aqua_select, names_from = Name, values_from = CladeCount)
    antarctica_embryophyta_aqua_taxatab[is.na(antarctica_embryophyta_aqua_taxatab)] <- 0 # replace NA with zero
    
    antarctica_embryophyta_aqua_counts <- as.data.frame(antarctica_embryophyta_aqua_taxatab[order(antarctica_embryophyta_aqua_taxatab$age), ])
    
    antarctica_embryophyta_aqua_perc <- data.frame(age = antarctica_embryophyta_aqua_counts[ ,1], 
                                               (antarctica_embryophyta_aqua_counts[ ,-1] / rowSums(antarctica_embryophyta_aqua_counts[ ,-1]) * 100))
    
    # -------------------------------------------------------------------------------------------------
    
    antarctica_embryophyta_aqua_perc_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(antarctica_embryophyta_aqua_perc, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-1])))
        
        timeslice_df <-unique(cbind(timeslice = paste0("Ant",timeslices[i]/1000),y1))
        
        antarctica_embryophyta_aqua_perc_timeslices <- rbind(antarctica_embryophyta_aqua_perc_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Beringia KL-12:

{
  
  shotgun_kl12 <- fread(paste0(inputpath,"KL12_APMG689/APMG689_nt0.5_Kraken2.txt")) # the readin data is a file of combined krona reports 
  names(shotgun_kl12) <- c("samples", "Percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")
  
  # load lineage using viktor's python script 
  lineage_kl12 <- fread(paste0(inputpath,"KL12_APMG689/APMG689_nt0.5_Kraken2_lineageDB.csv"))  
  
  # add lineage to the shotgun data
  shotgun.lineage_kl12 <- plyr::join(shotgun_kl12, lineage_kl12, by = "taxID")
  
  # load metadata containing sample name, age and depth
  metadata.readin_kl12 <- openxlsx::read.xlsx(xlsxFile = paste0(inputpath,"KL12_APMG689/APMG689_sample_name_change_nt05.xlsx"))
  colnames(metadata.readin_kl12)[3] <- "core"
  unique(metadata.readin_kl12$core) # check how many different cores in the sequencing run
  
  # -------------------------------------------------------------------------------------------------
  
  metadata_kl12 <- metadata.readin_kl12
  
  # remove whitespace at the end of the string:
  metadata_kl12$samples <- str_trim(metadata_kl12$samples)
  
  # add metadata to the shotgun data
  shotgun.join_kl12 <- plyr::join(shotgun.lineage_kl12, metadata_kl12, by = "samples")
  
  shotgun.join_kl12_samples <- shotgun.join_kl12 %>% filter(read_fraction == "merged" & type == "sample")
  
  # -------------------------------------------------------------------------------------------------
  
  shotgun.join_kl12_eukaryota <- shotgun.join_kl12_samples %>% filter(superkingdom == "Eukaryota")
  
  shotgun.join_kl12_eukaryota_total <- shotgun.join_kl12_eukaryota %>% 
    group_by(age) %>% 
    summarise(eukaryota_total = sum(CladeCount))
  
  shotgun.join_kl12_streptophyta <- shotgun.join_kl12_eukaryota %>% filter(phylum == "Streptophyta")
  
  shotgun.join_kl12_streptophyta_total <- shotgun.join_kl12_streptophyta %>% 
    group_by(age) %>% 
    summarise(streptophyta_total = sum(CladeCount))
  
  shotgun.join_kl12_other <- shotgun.join_kl12_eukaryota %>% filter(!phylum == "Streptophyta")
  
  shotgun.join_kl12_other_total <- shotgun.join_kl12_other %>% 
    group_by(age) %>% 
    summarise(other_total = sum(CladeCount))
  
  shotgun.join_kl12_embryophyta  <- shotgun.join_kl12_streptophyta %>% filter(!class %in% target)
  
  shotgun.join_kl12_embryophyta_terra <- shotgun.join_kl12_embryophyta %>% filter(!family %in% coastal)
  shotgun.join_kl12_embryophyta_aqua  <- shotgun.join_kl12_embryophyta %>% filter(family %in% coastal)
  
  shotgun.join_kl12_embryophyta_terra_total <- shotgun.join_kl12_embryophyta_terra %>% 
    group_by(age) %>% 
    summarise(embryophyta_terrestric_total = sum(CladeCount))
  
  shotgun.join_kl12_embryophyta_aqua_total <- shotgun.join_kl12_embryophyta_aqua %>% 
    group_by(age) %>% 
    summarise(embryophyta_aquatic_total = sum(CladeCount))
  
  # -------------------------------------------------------------------------------------------------
  
  kamchatka_eukaryota_embryophyta_df1 <- left_join(shotgun.join_kl12_eukaryota_total, shotgun.join_kl12_streptophyta_total)
  kamchatka_eukaryota_embryophyta_df2 <- left_join(kamchatka_eukaryota_embryophyta_df1, shotgun.join_kl12_other_total)
  kamchatka_eukaryota_embryophyta_df3 <- left_join(shotgun.join_kl12_embryophyta_terra_total, shotgun.join_kl12_embryophyta_aqua_total)
  
  kamchatka_eukaryota_embryophyta_df  <- left_join(kamchatka_eukaryota_embryophyta_df2, kamchatka_eukaryota_embryophyta_df3)
  
  kamchatka_eukaryota_embryophyta_df[is.na(kamchatka_eukaryota_embryophyta_df)] <- 0
  
  kamchatka_eukaryota_embryophyta_df <- kamchatka_eukaryota_embryophyta_df %>% 
    group_by(age) %>% 
    mutate(embryophyta_total            = (embryophyta_terrestric_total + embryophyta_aquatic_total),
           ratio_eukaryota_streptophyta = (streptophyta_total/eukaryota_total)*100,
           ratio_eukaryota_other        = (other_total/eukaryota_total)*100,
           ratio_eukaryota_terrestric   = (embryophyta_terrestric_total/eukaryota_total)*100,
           ratio_eukaryota_aquatic      = (embryophyta_aquatic_total/eukaryota_total)*100,
           ratio_eukaryota_embryophyta  = (embryophyta_total/eukaryota_total)*100) 
  
  kamchatka_eukaryota_embryophyta_df$core_name <- "off-Kamchatka"
  
  # -------------------------------------------------------------------------------------------------
  
  kamchatka_eukaryota_embryophyta_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(kamchatka_eukaryota_embryophyta_df, age >= timeslices[i] & age < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-c(1,13)])))
      
      timeslice_df <-unique(cbind(timeslice=(timeslices[i]/1000),y1,core_name=x1$core_name))
      
      kamchatka_eukaryota_embryophyta_timeslices <- rbind(kamchatka_eukaryota_embryophyta_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # write.csv(kamchatka_eukaryota_embryophyta_timeslices, file="C:/Data/Marine_cores/2023-06-01_eukaryota_embryophyta_kamchatka.csv", row.names=FALSE)
  
  # -------------------------------------------------------------------------------------------------
  
  # for PCA - terrestrial:
  
  {
    
    shotgun.join_kl12_embryophyta_terra_families <- shotgun.join_kl12_embryophyta_terra %>% filter(Rank == "F")
    
    kamchatka_embryophyta_terra_select <- select(shotgun.join_kl12_embryophyta_terra_families, c("age", "Name", "CladeCount"))
    
    kamchatka_embryophyta_terra_taxatab <- pivot_wider(kamchatka_embryophyta_terra_select, names_from = Name, values_from = CladeCount)
    kamchatka_embryophyta_terra_taxatab[is.na(kamchatka_embryophyta_terra_taxatab)] <- 0 # replace NA with zero
    
    kamchatka_embryophyta_terra_counts <- as.data.frame(kamchatka_embryophyta_terra_taxatab[order(kamchatka_embryophyta_terra_taxatab$age), ])
    
    kamchatka_embryophyta_terra_perc <- data.frame(age = kamchatka_embryophyta_terra_counts[ ,1], 
                                                (kamchatka_embryophyta_terra_counts[ ,-1] / rowSums(kamchatka_embryophyta_terra_counts[ ,-1]) * 100))
    
    # -------------------------------------------------------------------------------------------------
    
    kamchatka_embryophyta_terra_perc_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(kamchatka_embryophyta_terra_perc, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-1])))
        
        timeslice_df <-unique(cbind(timeslice = paste0("Ka",timeslices[i]/1000),y1))
        
        kamchatka_embryophyta_terra_perc_timeslices <- rbind(kamchatka_embryophyta_terra_perc_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # for PCA - aquatic:
  
  {
    
    shotgun.join_kl12_embryophyta_aqua_families <- shotgun.join_kl12_embryophyta_aqua %>% filter(Rank == "F")
    
    kamchatka_embryophyta_aqua_select <- select(shotgun.join_kl12_embryophyta_aqua_families, c("age", "Name", "CladeCount"))
    
    kamchatka_embryophyta_aqua_taxatab <- pivot_wider(kamchatka_embryophyta_aqua_select, names_from = Name, values_from = CladeCount)
    kamchatka_embryophyta_aqua_taxatab[is.na(kamchatka_embryophyta_aqua_taxatab)] <- 0 # replace NA with zero
    
    kamchatka_embryophyta_aqua_counts <- as.data.frame(kamchatka_embryophyta_aqua_taxatab[order(kamchatka_embryophyta_aqua_taxatab$age), ])
    
    kamchatka_embryophyta_aqua_perc <- data.frame(age = kamchatka_embryophyta_aqua_counts[ ,1], 
                                               (kamchatka_embryophyta_aqua_counts[ ,-1] / rowSums(kamchatka_embryophyta_aqua_counts[ ,-1]) * 100))
    
    # -------------------------------------------------------------------------------------------------
    
    kamchatka_embryophyta_aqua_perc_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(kamchatka_embryophyta_aqua_perc, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-1])))
        
        timeslice_df <-unique(cbind(timeslice = paste0("Ka",timeslices[i]/1000),y1))
        
        kamchatka_embryophyta_aqua_perc_timeslices <- rbind(kamchatka_embryophyta_aqua_perc_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
   
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Beringia KL-77:

{
  
  # load shotgun data from concatenated kraken reports
  shotgun.one_kl77 <- read.table(paste0(inputpath,"KL77_NextSeqBHV/KL77-1_nt0.5_Kraken2.txt"), header = FALSE, sep = "\t", quote = "") # the readin data is a file of combined krona reports 
  names(shotgun.one_kl77) <- c("samples", "Percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")
  shotgun.two_kl77 <- read.table(paste0(inputpath,"KL77_NextSeqBHV/KL77-2_nt0.5_Kraken2.txt"), header = FALSE, sep = "\t", quote = "") # the readin data is a file of combined krona reports 
  names(shotgun.two_kl77) <- c("samples", "Percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")
  shotgun_kl77 <- rbind(shotgun.one_kl77, shotgun.two_kl77)
  
  # load lineage using viktor's python script 
  lineage.one_kl77 <- read.csv(paste0(inputpath,"KL77_NextSeqBHV/KL77-1_nt0.5_Kraken2_lineageDB.csv"))  
  lineage.two_kl77 <- read.csv(paste0(inputpath,"KL77_NextSeqBHV/KL77-2_nt0.5_Kraken2_lineageDB.csv"))  
  lineage_kl77 <- rbind(lineage.one_kl77, lineage.two_kl77)
  
  # add lineage to the shotgun data
  shotgun.lineage_kl77 <- plyr::join(shotgun_kl77, lineage_kl77, by = "taxID")
  
  # load metadata containing sample name, age and depth
  metadata.one_kl77 <- openxlsx::read.xlsx(xlsxFile = paste0(inputpath,"KL77_NextSeqBHV/KL77-1_sample_name_change_nt0.5.xlsx"))
  metadata.one_kl77$core <- "KL77-1"
  metadata.two_kl77 <- openxlsx::read.xlsx(xlsxFile = paste0(inputpath,"KL77_NextSeqBHV/KL77-1_sample_name_change_nt0.5.xlsx"))
  metadata.two_kl77$core <- "KL77-2"
  metadata.readin_kl77 <- rbind(metadata.one_kl77, metadata.two_kl77)
  unique(metadata.readin_kl77$core) # check how many different cores in the sequencing run
  
  # -------------------------------------------------------------------------------------------------
  
  metadata_kl77 <- metadata.readin_kl77
  
  # remove whitespace at the end of the string:
  metadata_kl77$samples <- str_trim(metadata_kl77$samples)
  
  # add metadata to the shotgun data
  shotgun.join_kl77 <- plyr::join(shotgun.lineage_kl77, metadata_kl77, by = "samples")
  
  shotgun.join_kl77_samples <- shotgun.join_kl77 %>% filter(read_fraction == "merged" & type == "sample")
  
  # -------------------------------------------------------------------------------------------------
  
  shotgun.join_kl77_eukaryota <- shotgun.join_kl77_samples %>% filter(superkingdom == "Eukaryota")
  
  shotgun.join_kl77_eukaryota_total <- shotgun.join_kl77_eukaryota %>% 
    group_by(age) %>% 
    summarise(eukaryota_total = sum(CladeCount))
  
  shotgun.join_kl77_streptophyta <- shotgun.join_kl77_eukaryota %>% filter(phylum == "Streptophyta")
  
  shotgun.join_kl77_streptophyta_total <- shotgun.join_kl77_streptophyta %>% 
    group_by(age) %>% 
    summarise(streptophyta_total = sum(CladeCount))
  
  shotgun.join_kl77_other <- shotgun.join_kl77_eukaryota %>% filter(!phylum == "Streptophyta")
  
  shotgun.join_kl77_other_total <- shotgun.join_kl77_other %>% 
    group_by(age) %>% 
    summarise(other_total = sum(CladeCount))
  
  shotgun.join_kl77_embryophyta  <- shotgun.join_kl77_streptophyta %>% filter(!class %in% target)
  
  shotgun.join_kl77_embryophyta_terra <- shotgun.join_kl77_embryophyta %>% filter(!family %in% coastal)
  shotgun.join_kl77_embryophyta_aqua  <- shotgun.join_kl77_embryophyta %>% filter(family %in% coastal)
  
  shotgun.join_kl77_embryophyta_terra_total <- shotgun.join_kl77_embryophyta_terra %>% 
    group_by(age) %>% 
    summarise(embryophyta_terrestric_total = sum(CladeCount))
  
  shotgun.join_kl77_embryophyta_aqua_total <- shotgun.join_kl77_embryophyta_aqua %>% 
    group_by(age) %>% 
    summarise(embryophyta_aquatic_total = sum(CladeCount))
  
  # -------------------------------------------------------------------------------------------------
  
  beringia_eukaryota_embryophyta_df1 <- left_join(shotgun.join_kl77_eukaryota_total, shotgun.join_kl77_streptophyta_total)
  beringia_eukaryota_embryophyta_df2 <- left_join(beringia_eukaryota_embryophyta_df1, shotgun.join_kl77_other_total)
  beringia_eukaryota_embryophyta_df3 <- left_join(shotgun.join_kl77_embryophyta_terra_total, shotgun.join_kl77_embryophyta_aqua_total)
  
  beringia_eukaryota_embryophyta_df  <- left_join(beringia_eukaryota_embryophyta_df2, beringia_eukaryota_embryophyta_df3)
  
  beringia_eukaryota_embryophyta_df[is.na(beringia_eukaryota_embryophyta_df)] <- 0
  
  beringia_eukaryota_embryophyta_df <- beringia_eukaryota_embryophyta_df %>% 
    group_by(age) %>% 
    mutate(embryophyta_total            = (embryophyta_terrestric_total + embryophyta_aquatic_total),
           ratio_eukaryota_streptophyta = (streptophyta_total/eukaryota_total)*100,
           ratio_eukaryota_other        = (other_total/eukaryota_total)*100,
           ratio_eukaryota_terrestric   = (embryophyta_terrestric_total/eukaryota_total)*100,
           ratio_eukaryota_aquatic      = (embryophyta_aquatic_total/eukaryota_total)*100,
           ratio_eukaryota_embryophyta  = (embryophyta_total/eukaryota_total)*100) 
  
  beringia_eukaryota_embryophyta_df$core_name <- "Beringian Sea"
  
  # -------------------------------------------------------------------------------------------------
  
  beringia_eukaryota_embryophyta_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(beringia_eukaryota_embryophyta_df, age >= timeslices[i] & age < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-c(1,13)])))
      
      timeslice_df <-unique(cbind(timeslice=(timeslices[i]/1000),y1,core_name=x1$core_name))
      
      beringia_eukaryota_embryophyta_timeslices <- rbind(beringia_eukaryota_embryophyta_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # write.csv(beringia_eukaryota_embryophyta_timeslices, file="C:/Data/Marine_cores/2023-06-01_eukaryota_embryophyta_beringia.csv", row.names=FALSE)
  
  # -------------------------------------------------------------------------------------------------
  
  # for PCA - terrestrial:
  
  {
    
    shotgun.join_kl77_embryophyta_terra_families <- shotgun.join_kl77_embryophyta_terra %>% filter(Rank == "F")
    
    beringia_embryophyta_terra_select <- select(shotgun.join_kl77_embryophyta_terra_families, c("age", "Name", "CladeCount"))
    beringia_embryophyta_terra_select <- unique(beringia_embryophyta_terra_select)
    
    beringia_embryophyta_terra_taxatab <- pivot_wider(beringia_embryophyta_terra_select, names_from = Name, values_from = CladeCount)
    beringia_embryophyta_terra_taxatab[is.na(beringia_embryophyta_terra_taxatab)] <- 0 # replace NA with zero
    
    beringia_embryophyta_terra_counts <- as.data.frame(beringia_embryophyta_terra_taxatab[order(beringia_embryophyta_terra_taxatab$age), ])
    
    beringia_embryophyta_terra_perc <- data.frame(age = beringia_embryophyta_terra_counts[ ,1], 
                                                (beringia_embryophyta_terra_counts[ ,-1] / rowSums(beringia_embryophyta_terra_counts[ ,-1]) * 100))
    
    # -------------------------------------------------------------------------------------------------
    
    beringia_embryophyta_terra_perc_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(beringia_embryophyta_terra_perc, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-1])))
        
        timeslice_df <-unique(cbind(timeslice = paste0("BS",timeslices[i]/1000),y1))
        
        beringia_embryophyta_terra_perc_timeslices <- rbind(beringia_embryophyta_terra_perc_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # for PCA - aquatic:
  
  {
    
    shotgun.join_kl77_embryophyta_aqua_families <- shotgun.join_kl77_embryophyta_aqua %>% filter(Rank == "F")
    
    beringia_embryophyta_aqua_select <- select(shotgun.join_kl77_embryophyta_aqua_families, c("age", "Name", "CladeCount"))
    beringia_embryophyta_aqua_select <- unique(beringia_embryophyta_aqua_select)
    
    beringia_embryophyta_aqua_taxatab <- pivot_wider(beringia_embryophyta_aqua_select, names_from = Name, values_from = CladeCount)
    beringia_embryophyta_aqua_taxatab[is.na(beringia_embryophyta_aqua_taxatab)] <- 0 # replace NA with zero
    
    beringia_embryophyta_aqua_counts <- as.data.frame(beringia_embryophyta_aqua_taxatab[order(beringia_embryophyta_aqua_taxatab$age), ])
    
    beringia_embryophyta_aqua_perc <- data.frame(age = beringia_embryophyta_aqua_counts[ ,1], 
                                               (beringia_embryophyta_aqua_counts[ ,-1] / rowSums(beringia_embryophyta_aqua_counts[ ,-1]) * 100))
    
    # -------------------------------------------------------------------------------------------------
    
    beringia_embryophyta_aqua_perc_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(beringia_embryophyta_aqua_perc, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-1])))
        
        timeslice_df <-unique(cbind(timeslice = paste0("BS",timeslices[i]/1000),y1))
        
        beringia_embryophyta_aqua_perc_timeslices <- rbind(beringia_embryophyta_aqua_perc_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Fram Strait (other sequencing):

{
  
  shotgun_fram2 <- fread("//projects/biodiv/user/shotgun_data/05.2_shotgun_data_nt05/Lama/Lama_BHV2/lama_bhv2_nt_0.5.txt") # the readin data is a file of combined krona reports 
  names(shotgun_fram2) <- c("samples", "Percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")
  
  # load lineage using viktor's python script 
  lineage_fram2 <- fread("//smb.isipd.dmawi.de/projects/biodiv/user/shotgun_data/05.2_shotgun_data_nt05/Lama/Lama_BHV2/lineage_lama_bhv2_nt_0.5.csv")  
  
  # add lineage to the shotgun data
  shotgun.lineage_fram2 <- plyr::join(shotgun_fram2, lineage_fram2, by = "taxID")
  
  # load metadata containing sample name, age and depth
  metadata.readin_fram2 <- openxlsx::read.xlsx(xlsxFile = "//projects/biodiv/user/shotgun_data/05.2_shotgun_data_nt05/MSM_NextSeqBHV/MSM_nt05_shotgun_sample_name_change.xlsx")
  unique(metadata.readin_fram2$core) # check how many different cores in the sequencing run
  
  # -------------------------------------------------------------------------------------------------
  
  metadata_fram2 <- metadata.readin_fram2
  
  # remove whitespace at the end of the string:
  metadata_fram2$samples <- str_trim(metadata_fram2$samples)
  
  # add metadata to the shotgun data
  shotgun.join_fram2 <- plyr::join(shotgun.lineage_fram2, metadata_fram2, by = "samples")
  
  shotgun.join_fram2_samples <- shotgun.join_fram2 %>% filter(read_fraction == "merged" & type == "sample")
  
  shotgun.join_fram2_samples <- shotgun.join_fram2_samples %>% 
    select(samples, Percentage, CladeCount, TaxCount, Rank, taxID, Name, superkingdom, phylum, class, order, family, genus, species, read_fraction, site, core, depth, age, lib_id, lib_batch, extract, type, sample_name)
  
  shotgun.join_framstrait_samples <- rbind(shotgun.join_fram_samples, shotgun.join_fram2_samples)
  
  # -------------------------------------------------------------------------------------------------
  
  shotgun.join_framstrait_eukaryota <- shotgun.join_framstrait_samples %>% filter(superkingdom == "Eukaryota")
  
  shotgun.join_framstrait_eukaryota_total <- shotgun.join_framstrait_eukaryota %>% 
    group_by(age) %>% 
    summarise(eukaryota_total = sum(CladeCount))
  
  shotgun.join_framstrait_streptophyta <- shotgun.join_framstrait_eukaryota %>% filter(phylum == "Streptophyta")
  
  shotgun.join_framstrait_streptophyta_total <- shotgun.join_framstrait_streptophyta %>% 
    group_by(age) %>% 
    summarise(streptophyta_total = sum(CladeCount))
  
  shotgun.join_framstrait_other <- shotgun.join_framstrait_eukaryota %>% filter(!phylum == "Streptophyta")
  
  shotgun.join_framstrait_other_total <- shotgun.join_framstrait_other %>% 
    group_by(age) %>% 
    summarise(other_total = sum(CladeCount))
  
  shotgun.join_framstrait_embryophyta  <- shotgun.join_framstrait_streptophyta %>% filter(!class %in% target)
  
  shotgun.join_framstrait_embryophyta_terra <- shotgun.join_framstrait_embryophyta %>% filter(!family %in% coastal)
  shotgun.join_framstrait_embryophyta_aqua  <- shotgun.join_framstrait_embryophyta %>% filter(family %in% coastal)
  
  shotgun.join_framstrait_embryophyta_terra_total <- shotgun.join_framstrait_embryophyta_terra %>% 
    group_by(age) %>% 
    summarise(embryophyta_terrestric_total = sum(CladeCount))
  
  shotgun.join_framstrait_embryophyta_aqua_total <- shotgun.join_framstrait_embryophyta_aqua %>% 
    group_by(age) %>% 
    summarise(embryophyta_aquatic_total = sum(CladeCount))
  
  # -------------------------------------------------------------------------------------------------
  
  framstrait_eukaryota_embryophyta_df1 <- left_join(shotgun.join_framstrait_eukaryota_total, shotgun.join_framstrait_streptophyta_total)
  framstrait_eukaryota_embryophyta_df2 <- left_join(framstrait_eukaryota_embryophyta_df1, shotgun.join_framstrait_other_total)
  framstrait_eukaryota_embryophyta_df3 <- left_join(shotgun.join_framstrait_embryophyta_terra_total, shotgun.join_framstrait_embryophyta_aqua_total)
  
  framstrait_eukaryota_embryophyta_df  <- left_join(framstrait_eukaryota_embryophyta_df2, framstrait_eukaryota_embryophyta_df3)
  
  framstrait_eukaryota_embryophyta_df[is.na(framstrait_eukaryota_embryophyta_df)] <- 0
  
  framstrait_eukaryota_embryophyta_df <- framstrait_eukaryota_embryophyta_df %>% 
    group_by(age) %>% 
    mutate(embryophyta_total            = (embryophyta_terrestric_total + embryophyta_aquatic_total),
           ratio_eukaryota_streptophyta = (streptophyta_total/eukaryota_total)*100,
           ratio_eukaryota_other        = (other_total/eukaryota_total)*100,
           ratio_eukaryota_terrestric   = (embryophyta_terrestric_total/eukaryota_total)*100,
           ratio_eukaryota_aquatic      = (embryophyta_aquatic_total/eukaryota_total)*100,
           ratio_eukaryota_embryophyta  = (embryophyta_total/eukaryota_total)*100) 
  
  framstrait_eukaryota_embryophyta_df$core_name <- "Fram Strait"
  
  # -------------------------------------------------------------------------------------------------
  
  framstrait_eukaryota_embryophyta_timeslices <- NULL
  
  
  for(i in 1:(length(timeslices)-1)){
    
    x1 <- subset(framstrait_eukaryota_embryophyta_df, age >= timeslices[i] & age < timeslices[i+1])
    
    if(dim(x1)[1] > 0){
      
      y1 <- as.data.frame(t(colMeans(x1[,-c(1,13)])))
      
      timeslice_df <-unique(cbind(timeslice=(timeslices[i]/1000),y1,core_name=x1$core_name))
      
      framstrait_eukaryota_embryophyta_timeslices <- rbind(framstrait_eukaryota_embryophyta_timeslices, timeslice_df)
      
    }
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # write.csv(framstrait_eukaryota_embryophyta_timeslices, file="C:/Data/Marine_cores/2023-06-01_eukaryota_embryophyta_framstrait.csv", row.names=FALSE)
  
  # -------------------------------------------------------------------------------------------------
  
  # for PCA - terrestrial:
  
  {
    
    shotgun.join_framstrait_embryophyta_terra_families <- shotgun.join_framstrait_embryophyta_terra %>% filter(Rank == "F")
    
    framstrait_embryophyta_terra_select <- select(shotgun.join_framstrait_embryophyta_terra_families, c("age", "Name", "CladeCount"))
    
    framstrait_embryophyta_terra_taxatab <- pivot_wider(framstrait_embryophyta_terra_select, names_from = Name, values_from = CladeCount)
    framstrait_embryophyta_terra_taxatab[is.na(framstrait_embryophyta_terra_taxatab)] <- 0 # replace NA with zero
    
    framstrait_embryophyta_terra_counts <- as.data.frame(framstrait_embryophyta_terra_taxatab[order(framstrait_embryophyta_terra_taxatab$age), ])
    
    framstrait_embryophyta_terra_perc <- data.frame(age = framstrait_embryophyta_terra_counts[ ,1], 
                                                (framstrait_embryophyta_terra_counts[ ,-1] / rowSums(framstrait_embryophyta_terra_counts[ ,-1]) * 100))
    
    # -------------------------------------------------------------------------------------------------
    
    framstrait_embryophyta_terra_perc_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(framstrait_embryophyta_terra_perc, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-1])))
        
        timeslice_df <-unique(cbind(timeslice = paste0("FS",timeslices[i]/1000),y1))
        
        framstrait_embryophyta_terra_perc_timeslices <- rbind(framstrait_embryophyta_terra_perc_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
  # for PCA - aquatic:
  
  {
    
    shotgun.join_framstrait_embryophyta_aqua_families <- shotgun.join_framstrait_embryophyta_aqua %>% filter(Rank == "F")
    
    framstrait_embryophyta_aqua_select <- select(shotgun.join_framstrait_embryophyta_aqua_families, c("age", "Name", "CladeCount"))
    
    framstrait_embryophyta_aqua_taxatab <- pivot_wider(framstrait_embryophyta_aqua_select, names_from = Name, values_from = CladeCount)
    framstrait_embryophyta_aqua_taxatab[is.na(framstrait_embryophyta_aqua_taxatab)] <- 0 # replace NA with zero
    
    framstrait_embryophyta_aqua_counts <- as.data.frame(framstrait_embryophyta_aqua_taxatab[order(framstrait_embryophyta_aqua_taxatab$age), ])
    
    framstrait_embryophyta_aqua_perc <- data.frame(age = framstrait_embryophyta_aqua_counts[ ,1], 
                                               (framstrait_embryophyta_aqua_counts[ ,-1] / rowSums(framstrait_embryophyta_aqua_counts[ ,-1]) * 100))
    
    # -------------------------------------------------------------------------------------------------
    
    framstrait_embryophyta_aqua_perc_timeslices <- NULL
    
    
    for(i in 1:(length(timeslices)-1)){
      
      x1 <- subset(framstrait_embryophyta_aqua_perc, age >= timeslices[i] & age < timeslices[i+1])
      
      if(dim(x1)[1] > 0){
        
        y1 <- as.data.frame(t(colMeans(x1[,-1])))
        
        timeslice_df <-unique(cbind(timeslice = paste0("FS",timeslices[i]/1000),y1))
        
        framstrait_embryophyta_aqua_perc_timeslices <- rbind(framstrait_embryophyta_aqua_perc_timeslices, timeslice_df)
        
      }
    }
    
    # -------------------------------------------------------------------------------------------------
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### combine marine datasets:

marine_cores_ratios_combined <- rbind(tobago_eukaryota_embryophyta_timeslices,
                                      australia_eukaryota_embryophyta_timeslices,
                                      antarctica_eukaryota_embryophyta_timeslices,
                                      kamchatka_eukaryota_embryophyta_timeslices,
                                      beringia_eukaryota_embryophyta_timeslices,
                                      framstrait_eukaryota_embryophyta_timeslices)

marine_cores_ratios_combined <- as.data.frame(marine_cores_ratios_combined)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### Stratigraphic plots:

{

level_order <- c("Bransfield Strait", "off-Australia", "Fram Strait", "Beringian Sea", "off-Kamchatka", "off-Tobago")  

colramp <- c("#c8586c","#a3ad62","#bd925a","#C6D9F1","#17365D","#86776F") 

# -------------------------------------------------------------------------------------------------

marine_cores_ratios_combined <- marine_cores_ratios_combined %>% 
  mutate(age = case_when(core_name == "Bransfield Strait" ~ timeslice+0.3,
                         core_name == "Fram Strait" ~ timeslice+0.6,
                         core_name == "off-Kamchatka" ~ timeslice+0.9,
                         core_name == "off-Australia" ~ timeslice+1.2,
                         core_name == "Beringian Sea" ~ timeslice+1.5,
                         core_name == "off-Tobago" ~ timeslice+1.8)) 

marine_cores_ratios_combined <- marine_cores_ratios_combined[order(marine_cores_ratios_combined$age), ]


## terrestrial plants:

{

terrestrial_ratio <- ggplot(marine_cores_ratios_combined, aes(x=age, y=ratio_eukaryota_terrestric,fill=core_name)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
  labs(x="",y="terrestrial") + 
  easy_remove_x_axis() + 
  theme(plot.title=element_text(size=10),
        plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title=element_text(size=10),
        axis.text.y=element_text(size=10)) + 
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
  theme(panel.background = element_blank())+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 


terrestrial_ratio1 <- terrestrial_ratio + scale_fill_manual(name = "", values = colramp, breaks = level_order) 

}

# -------------------------------------------------------------------------------------------------

## aquatic plants:

{

  ### highest value (6.67%) is manually set to 0.5%:
  marine_cores_ratios_combined$ratio_eukaryota_aquatic[marine_cores_ratios_combined$ratio_eukaryota_aquatic > 0.5] <- 0.5
  
  
aquatic_ratio <- ggplot(marine_cores_ratios_combined, aes(x=age, y=ratio_eukaryota_aquatic,fill=core_name)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
  labs(x="",y="aquatic") + 
  #easy_remove_x_axis() + 
  theme(plot.title=element_text(size=10),
        plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10)) + 
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
  theme(panel.background = element_blank())+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(axis.line.x = element_line(color="grey50", linewidth = 0.2), 
        axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 


aquatic_ratio1 <- aquatic_ratio + scale_fill_manual(name = "", values = colramp, breaks = level_order) 

}

# -------------------------------------------------------------------------------------------------

##### plot all #######
plot_grid(terrestrial_ratio1,
          aquatic_ratio1,
          ncol=1,nrow=2)


ggsave("C:/Data/2023-06-01_stratigraphic_plot_marine_cores_ratios_eukaryota_embryophyta_terrestrial_aquatic.png", width = 18, height = 9, units = "cm", dpi = 300)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

## terrestrial and aquatic combined as Embryophyta:

{
  
  # embryophyta_ratio <- ggplot(marine_cores_ratios_combined, aes(x=age, y=ratio_eukaryota_embryophyta,fill=core_name)) +
  #   geom_bar(stat="identity", width=0.3) + 
  #   scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
  #   labs(x="",y="% Embryophyta on Eukaryota") + 
  #   easy_remove_x_axis() + 
  #   theme(plot.title=element_text(size=10),
  #         plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
  #         axis.title=element_text(size=7),
  #         axis.text.y=element_text(size=10)) + 
  #   theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
  #   theme(panel.background = element_blank())+ 
  #   theme(panel.grid.major=element_blank(),
  #         panel.grid.minor=element_blank()) +
  #   theme(axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  # 
  # 
  # embryophyta_ratio1 <- embryophyta_ratio + scale_fill_manual(name = "", values = colramp, breaks = level_order) 
  
 
  embryophyta_ratio <- ggplot(marine_cores_ratios_combined, aes(x=age, y=ratio_eukaryota_embryophyta,fill=core_name)) +
    geom_bar(stat="identity", width=0.3) + 
    scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
    labs(x="",y="% Embryophyta on Eukaryota") + 
    #easy_remove_x_axis() + 
    theme(plot.title=element_text(size=10),
          plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
          axis.title=element_text(size=7),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10)) + 
    theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
    theme(panel.background = element_blank())+ 
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank()) +
    theme(axis.line.x = element_line(color="grey50", linewidth = 0.2), 
          axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 
  
  
  embryophyta_ratio1 <- embryophyta_ratio + scale_fill_manual(name = "", values = colramp, breaks = level_order) 
  
  
  ggsave("C:/Data/Marine_cores/2023-06-02_stratigraphic_plot_marine_cores_ratios_eukaryota_embryophyta.png", width = 18, height = 5, units = "cm", dpi = 300)
  
  
  
  
}

# -------------------------------------------------------------------------------------------------

## Aquatic taxa as percentages of Embryophyta:

{

marine_cores <- read.csv("/Volumes/projects/biodiv/user/tboehmer/Marine_cores/2023-05-30_shotgun_marine_cores_families_percentage_timeslices.csv", header=TRUE)

### Preparation marine cores:

  tsl <- marine_cores[,1]
  core_names <- data.frame(timeslice=tsl, startage=NA, core_name=NA)
  
  for(i in 1:length(tsl)){
    
    if(str_detect(tsl[i],"To")){core_names[i,"core_name"] <- "off-Tobago"
    }else if(str_detect(tsl[i],"Au")){core_names[i,"core_name"] <- "off-Australia"
    }else if(str_detect(tsl[i],"Ant")){core_names[i,"core_name"] <- "Bransfield Strait"
    }else if(str_detect(tsl[i],"Ka")){core_names[i,"core_name"] <- "off-Kamchatka"
    }else if(str_detect(tsl[i],"BS")){core_names[i,"core_name"] <- "Beringian Sea"
    }else if(str_detect(tsl[i],"FS")){core_names[i,"core_name"] <- "Fram Strait"
    }
    
    core_names[i,"startage"] <- as.numeric(str_extract(tsl[i], "[0-9]+"))
    
  }
  
  # -------------------------------------------------------------------------------------------------
  
aquatic_marine_cores <- marine_cores[ ,c(which(names(marine_cores) %in% coastal))]
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


Aquatics <- ggplot(Aquatic_df, aes(x=age, y=Aquatic,fill=core_name)) +
  geom_bar(stat="identity", width=0.3) + 
  scale_x_continuous(limits=c(0,50), breaks=round(seq(min(0), max(50), by =2),1), labels=round(seq(min(0), max(50), by =2),1), expand=c(0.025,0)) +
  #scale_y_continuous(limits=c(0,35), breaks=round(seq(min(0), max(40), by =20),1), labels=round(seq(min(0), max(40), by =20),1), expand=c(0,0)) +
  labs(x="",y="% Aquatic taxa on Embryophyta") + 
  #easy_remove_x_axis() + 
  theme(plot.title=element_text(size=10),
        plot.margin=margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title=element_text(size=7),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10)) + 
  theme(plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm")) +
  theme(panel.background = element_blank())+ 
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) +
  theme(axis.line.x = element_line(color="grey50", linewidth = 0.2), 
        axis.line.y = element_line(color="grey50", linewidth = 0.2), legend.position = "none") 


Aquatics1 <- Aquatics + scale_fill_manual(name = "", values = colramp, breaks = level_order) 

}

# -------------------------------------------------------------------------------------------------

##### plot all #######
plot_grid(embryophyta_ratio1,
          Aquatics1,
          ncol=1,nrow=2)


ggsave("C:/Data/Marine_cores/2023-06-02_stratigraphic_plot_marine_cores_ratios_eukaryota_embryophyta_aquatics.png", width = 18, height = 9, units = "cm", dpi = 300)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

### PCA:

{

library(analogue)
library(ggplot2)
library(ggrepel)

## terrestrial:

{

marine_cores_embryophyta_terra_combined <- dplyr::bind_rows(tobago_embryophyta_terra_perc_timeslices,
                                                            australia_embryophyta_terra_perc_timeslices,
                                                            antarctica_embryophyta_terra_perc_timeslices,
                                                            kamchatka_embryophyta_terra_perc_timeslices,
                                                            beringia_embryophyta_terra_perc_timeslices,
                                                            framstrait_embryophyta_terra_perc_timeslices)


marine_cores_embryophyta_terra_combined[is.na(marine_cores_embryophyta_terra_combined)] <- 0 # replace NA with zero

# -------------------------------------------------------------------------------------------------

marine_cores_embryophyta_terra_percmeans <- marine_cores_embryophyta_terra_combined[,-1] # remove column of age

cores <- marine_cores_embryophyta_terra_combined[,1]
core_names <- data.frame(timeslice=cores, core_name=NA)

for(i in 1:length(cores)){

  if(str_detect(cores[i],"To")){core_names[i,"core_name"] <- "off-Tobago"
  }else if(str_detect(cores[i],"Au")){core_names[i,"core_name"] <- "off-Australia"
  }else if(str_detect(cores[i],"Ant")){core_names[i,"core_name"] <- "Bransfield Strait"
  }else if(str_detect(cores[i],"Ka")){core_names[i,"core_name"] <- "off-Kamchatka"
  }else if(str_detect(cores[i],"BS")){core_names[i,"core_name"] <- "Beringian Sea"
  }else if(str_detect(cores[i],"FS")){core_names[i,"core_name"] <- "Fram Strait"
  }
  
}

# -------------------------------------------------------------------------------------------------

# double square-root:

{

ds <- sqrt(sqrt(marine_cores_embryophyta_terra_percmeans)) # double square root to transform

marine_cores_terra_ds_pca.out <- rda(ds, scale = FALSE)

# extract the results from pca
marine_cores_terra_ds_sites <- as.data.frame(scores(marine_cores_terra_ds_pca.out, display = "sites", scaling = 1))
marine_cores_terra_ds_species <- as.data.frame(scores(marine_cores_terra_ds_pca.out, display = "species", scaling = 1))

# bind the data for plotting
marine_cores_terra_ds_bind.age <- as.data.frame(cbind(marine_cores_terra_ds_sites, marine_cores_embryophyta_terra_combined$timeslice))
names(marine_cores_terra_ds_bind.age)[3] <- "timeslice"
marine_cores_terra_ds_bind <- plyr::join(marine_cores_terra_ds_bind.age, core_names, by = "timeslice")
marine_cores_terra_ds_bind$PC1 <- as.numeric(as.character(marine_cores_terra_ds_bind$PC1))
marine_cores_terra_ds_bind$PC2 <- as.numeric(as.character(marine_cores_terra_ds_bind$PC2))


# plot ages and taxa in pca

ggplot(data = marine_cores_terra_ds_bind, aes(x = PC1, y = PC2)) +
  geom_text(data = marine_cores_terra_ds_species, aes(x = PC1, y = PC2, label = rownames(marine_cores_terra_ds_species)), size = 3) + # taxa names
  geom_text_repel(aes(label = timeslice, color = core_name), size = 4, max.overlaps = Inf) + # sample names in age
  scale_x_continuous(name = paste("PC1 (", round(marine_cores_terra_ds_pca.out$CA$eig[1] / marine_cores_terra_ds_pca.out$tot.chi*100, 1), "%)", sep = "")) + # x axis label
  scale_y_continuous(name = paste("PC2 (", round(marine_cores_terra_ds_pca.out$CA$eig[2] / marine_cores_terra_ds_pca.out$tot.chi*100, 1), "%)", sep = "")) + # y axis label
  geom_hline(yintercept = 0, alpha = 0.3) + geom_vline(xintercept = 0, alpha = 0.3) +
  theme_bw() + theme(panel.grid = element_blank()) + # remove background
  theme(text = element_text(size = 20),  #axis text size
        legend.title = element_text(size = 20), #legend title size
        legend.text = element_text(size = 20)) +
  ggtitle("Embryophyta terrestrial families in marine cores nt 0.5 scaling 1")


ggsave("C:/Data/Marine_cores/2023-06-01_pca_marine_cores_nt_0.5_embryophyta_terrestrial_family_scaling1_double_sqrt.png", width = 40, height = 25, units = "cm")

}

# -------------------------------------------------------------------------------------------------

# triple square-root:

{

ds2 <- sqrt(sqrt(sqrt(marine_cores_embryophyta_terra_percmeans))) # triple square root to transform

marine_cores_terra_ds2_pca.out <- rda(ds2, scale = FALSE)

# extract the results from pca
marine_cores_terra_ds2_sites <- as.data.frame(scores(marine_cores_terra_ds2_pca.out, display = "sites", scaling = 1))
marine_cores_terra_ds2_species <- as.data.frame(scores(marine_cores_terra_ds2_pca.out, display = "species", scaling = 1))

# bind the data for plotting
marine_cores_terra_ds2_bind.age <- as.data.frame(cbind(marine_cores_terra_ds2_sites, marine_cores_embryophyta_terra_combined$timeslice))
names(marine_cores_terra_ds2_bind.age)[3] <- "timeslice"
marine_cores_terra_ds2_bind <- plyr::join(marine_cores_terra_ds2_bind.age, core_names, by = "timeslice")
marine_cores_terra_ds2_bind$PC1 <- as.numeric(as.character(marine_cores_terra_ds2_bind$PC1))
marine_cores_terra_ds2_bind$PC2 <- as.numeric(as.character(marine_cores_terra_ds2_bind$PC2))


# plot ages and taxa in pca

ggplot(data = marine_cores_terra_ds2_bind, aes(x = PC1, y = PC2)) +
  geom_text(data = marine_cores_terra_ds2_species, aes(x = PC1, y = PC2, label = rownames(marine_cores_terra_ds2_species)), size = 3) + # taxa names
  geom_text_repel(aes(label = timeslice, color = core_name), size = 4, max.overlaps = Inf) + # sample names in age
  scale_x_continuous(name = paste("PC1 (", round(marine_cores_terra_ds2_pca.out$CA$eig[1] / marine_cores_terra_ds2_pca.out$tot.chi*100, 1), "%)", sep = "")) + # x axis label
  scale_y_continuous(name = paste("PC2 (", round(marine_cores_terra_ds2_pca.out$CA$eig[2] / marine_cores_terra_ds2_pca.out$tot.chi*100, 1), "%)", sep = "")) + # y axis label
  geom_hline(yintercept = 0, alpha = 0.3) + geom_vline(xintercept = 0, alpha = 0.3) +
  theme_bw() + theme(panel.grid = element_blank()) + # remove background
  theme(text = element_text(size = 20),  #axis text size
        legend.title = element_text(size = 20), #legend title size
        legend.text = element_text(size = 20)) +
  ggtitle("Embryophyta terrestrial families in marine cores nt 0.5 scaling 1")


ggsave("C:/Data/Marine_cores/2023-06-01_pca_marine_cores_nt_0.5_embryophyta_terrestrial_family_scaling1_triple_sqrt.png", width = 40, height = 25, units = "cm")

}

# -------------------------------------------------------------------------------------------------

}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

## aquatic:

{

marine_cores_embryophyta_aqua_combined <- dplyr::bind_rows(tobago_embryophyta_aqua_perc_timeslices,
                                                           australia_embryophyta_aqua_perc_timeslices,
                                                           antarctica_embryophyta_aqua_perc_timeslices,
                                                           kamchatka_embryophyta_aqua_perc_timeslices,
                                                           beringia_embryophyta_aqua_perc_timeslices,
                                                           framstrait_embryophyta_aqua_perc_timeslices)


marine_cores_embryophyta_aqua_combined[is.na(marine_cores_embryophyta_aqua_combined)] <- 0 # replace NA with zero

# -------------------------------------------------------------------------------------------------

marine_cores_embryophyta_aqua_percmeans <- marine_cores_embryophyta_aqua_combined[,-1] # remove column of age

cores <- marine_cores_embryophyta_aqua_combined[,1]
core_names <- data.frame(timeslice=cores, core_name=NA)

for(i in 1:length(cores)){
  
  if(str_detect(cores[i],"To")){core_names[i,"core_name"] <- "off-Tobago"
  }else if(str_detect(cores[i],"Au")){core_names[i,"core_name"] <- "off-Australia"
  }else if(str_detect(cores[i],"Ant")){core_names[i,"core_name"] <- "Bransfield Strait"
  }else if(str_detect(cores[i],"Ka")){core_names[i,"core_name"] <- "off-Kamchatka"
  }else if(str_detect(cores[i],"BS")){core_names[i,"core_name"] <- "Beringian Sea"
  }else if(str_detect(cores[i],"FS")){core_names[i,"core_name"] <- "Fram Strait"
  }
  
}

# -------------------------------------------------------------------------------------------------

# double square-root:

{
  
  ds <- sqrt(sqrt(marine_cores_embryophyta_aqua_percmeans)) # double square root to transform
  
  marine_cores_aqua_ds_pca.out <- rda(ds, scale = FALSE)
  
  # extract the results from pca
  marine_cores_aqua_ds_sites <- as.data.frame(scores(marine_cores_aqua_ds_pca.out, display = "sites", scaling = 1))
  marine_cores_aqua_ds_species <- as.data.frame(scores(marine_cores_aqua_ds_pca.out, display = "species", scaling = 1))
  
  # bind the data for plotting
  marine_cores_aqua_ds_bind.age <- as.data.frame(cbind(marine_cores_aqua_ds_sites, marine_cores_embryophyta_aqua_combined$timeslice))
  names(marine_cores_aqua_ds_bind.age)[3] <- "timeslice"
  marine_cores_aqua_ds_bind <- plyr::join(marine_cores_aqua_ds_bind.age, core_names, by = "timeslice")
  marine_cores_aqua_ds_bind$PC1 <- as.numeric(as.character(marine_cores_aqua_ds_bind$PC1))
  marine_cores_aqua_ds_bind$PC2 <- as.numeric(as.character(marine_cores_aqua_ds_bind$PC2))
  
  
  # plot ages and taxa in pca
  
  ggplot(data = marine_cores_aqua_ds_bind, aes(x = PC1, y = PC2)) +
    geom_text(data = marine_cores_aqua_ds_species, aes(x = PC1, y = PC2, label = rownames(marine_cores_aqua_ds_species)), size = 3) + # taxa names
    geom_text_repel(aes(label = timeslice, color = core_name), size = 4, max.overlaps = Inf) + # sample names in age
    scale_x_continuous(name = paste("PC1 (", round(marine_cores_aqua_ds_pca.out$CA$eig[1] / marine_cores_aqua_ds_pca.out$tot.chi*100, 1), "%)", sep = "")) + # x axis label
    scale_y_continuous(name = paste("PC2 (", round(marine_cores_aqua_ds_pca.out$CA$eig[2] / marine_cores_aqua_ds_pca.out$tot.chi*100, 1), "%)", sep = "")) + # y axis label
    geom_hline(yintercept = 0, alpha = 0.3) + geom_vline(xintercept = 0, alpha = 0.3) +
    theme_bw() + theme(panel.grid = element_blank()) + # remove background
    theme(text = element_text(size = 20),  #axis text size
          legend.title = element_text(size = 20), #legend title size
          legend.text = element_text(size = 20)) +
    ggtitle("Embryophyta aquatic families in marine cores nt 0.5 scaling 1")
  
  
  ggsave("C:/Data/Marine_cores/2023-06-01_pca_marine_cores_nt_0.5_embryophyta_aquatic_family_scaling1_double_sqrt.png", width = 40, height = 25, units = "cm")
  
}

# -------------------------------------------------------------------------------------------------

# triple square-root:

{
  
  ds2 <- sqrt(sqrt(sqrt(marine_cores_embryophyta_aqua_percmeans))) # triple square root to transform
  
  marine_cores_aqua_ds2_pca.out <- rda(ds2, scale = FALSE)
  
  # extract the results from pca
  marine_cores_aqua_ds2_sites <- as.data.frame(scores(marine_cores_aqua_ds2_pca.out, display = "sites", scaling = 1))
  marine_cores_aqua_ds2_species <- as.data.frame(scores(marine_cores_aqua_ds2_pca.out, display = "species", scaling = 1))
  
  # bind the data for plotting
  marine_cores_aqua_ds2_bind.age <- as.data.frame(cbind(marine_cores_aqua_ds2_sites, marine_cores_embryophyta_aqua_combined$timeslice))
  names(marine_cores_aqua_ds2_bind.age)[3] <- "timeslice"
  marine_cores_aqua_ds2_bind <- plyr::join(marine_cores_aqua_ds2_bind.age, core_names, by = "timeslice")
  marine_cores_aqua_ds2_bind$PC1 <- as.numeric(as.character(marine_cores_aqua_ds2_bind$PC1))
  marine_cores_aqua_ds2_bind$PC2 <- as.numeric(as.character(marine_cores_aqua_ds2_bind$PC2))
  
  
  # plot ages and taxa in pca
  
  ggplot(data = marine_cores_aqua_ds2_bind, aes(x = PC1, y = PC2)) +
    geom_text(data = marine_cores_aqua_ds2_species, aes(x = PC1, y = PC2, label = rownames(marine_cores_aqua_ds2_species)), size = 3) + # taxa names
    geom_text_repel(aes(label = timeslice, color = core_name), size = 4, max.overlaps = Inf) + # sample names in age
    scale_x_continuous(name = paste("PC1 (", round(marine_cores_aqua_ds2_pca.out$CA$eig[1] / marine_cores_aqua_ds2_pca.out$tot.chi*100, 1), "%)", sep = "")) + # x axis label
    scale_y_continuous(name = paste("PC2 (", round(marine_cores_aqua_ds2_pca.out$CA$eig[2] / marine_cores_aqua_ds2_pca.out$tot.chi*100, 1), "%)", sep = "")) + # y axis label
    geom_hline(yintercept = 0, alpha = 0.3) + geom_vline(xintercept = 0, alpha = 0.3) +
    theme_bw() + theme(panel.grid = element_blank()) + # remove background
    theme(text = element_text(size = 20),  #axis text size
          legend.title = element_text(size = 20), #legend title size
          legend.text = element_text(size = 20)) +
    ggtitle("Embryophyta aquatic families in marine cores nt 0.5 scaling 1")
  
  
  ggsave("C:/Data/Marine_cores/2023-06-01_pca_marine_cores_nt_0.5_embryophyta_aquatic_family_scaling1_triple_sqrt.png", width = 40, height = 25, units = "cm")
  
}

# -------------------------------------------------------------------------------------------------

}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

}

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------



