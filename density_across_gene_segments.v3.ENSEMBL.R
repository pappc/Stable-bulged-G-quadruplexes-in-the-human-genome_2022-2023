library(tidyverse)
library(data.table)
library(doParallel)
library(gridExtra)

setwd("/home/csaba/work/Projects/G4BS_new_results/New_G4BS_data/gene_segment_density/inputs/")
# setwd("/home/csaba/work/Projects/G4BS_new_results/external_datasets/Ensembl_genes/Gene_segment_density_inputs/")
# setwd("D:/Csaba_stuff/Shared_folder/OneDrive - SUNY Upstate Medical University/Projects/G4BS_new_results/New_G4BS_data/gene_segment_density/inputs/")
#Load data
files <- list.files()
data <- list()
for (i in 1:length(files)) {
  print(i)
  print(files[i])
  data[[i]] <- fread(files[i])
  names(data)[i] <- files[i]
  print(dim(data[[i]]))
}

# names(data) <- c("first_int", "first_exo", "second_int", "second_exo", "whole_gene", "UTR3", "UTR5", "last_exo")
names(data) <- c("first_exo", "first_int", "second_exo", "second_int", "whole_gene", "UTR3", "UTR5", "last_exo")

#Create list for holding final results. DO NOT RERUN DURING ANALYSIS OR RESULTS WILL BE LOST
final <- list()

#Load G4 data
# setwd("D:/Csaba_stuff/Shared_folder/OneDrive - SUNY Upstate Medical University/Projects/G4BS_new_results/")
setwd("~/work/Projects/G4BS_new_results/")
# g4bs <- fread("datasets/cG4_L7.Hg38.uniq.no_repEle.merged.txt", sep = "\t")
# g4bs <- fread("datasets/g4bs.Hg38.uniq.no_repEle.merged.txt", sep = "\t")
# g4bs <- fread("datasets/g4bs.Hg38.uniq.G3B1.no_repEle.merged.txt", sep = "\t")
# g4bs <- fread("datasets/g4bs.Hg38.uniq.G3B2.no_repEle.merged.txt", sep = "\t")
# g4bs <- fread("datasets/g4bs.Hg38.uniq.G2B2.no_repEle.merged.txt", sep = "\t")
# g4bs <- fread("datasets/g4seq.NEW.Hg38.no_repEle.bed", sep = "\t")
g4bs <- fread("datasets/RLFS.Hg38.merged.bed", sep = "\t")

g4bs <- g4bs %>% 
  select(chr = 1, start = 2, end = 3, strand = 6) %>%
  mutate(name = paste(chr, start, end, sep = ":")) %>% 
  data.table()

#Preprocess data
#Selecting correct columns and renaming
for (i in length(data)) {
  data[[i]] <- data[[i]] %>% 
    select(chr, strand, start, end, transcript_id, gene_symbol) %>% 
    data.table()
}
#Check for incorrect records
for (i in data) {
  asd <- i$end > i$start
  print(table(asd))
}
#Remove these records
for (i in 1:length(data)) {
  data[[i]] <- data[[i]] %>% 
    filter(end > start) %>% 
    data.table()
}

#Select our window of interest for each dataset
#The "get_coords" function can be used to select intervals around either the start or end of the given sequences
#We'll create separate datasets for the start and end sites of each dataset
get_coords <- function(df, reference = "start", interval = 250) {
  if (reference == "center") {
    res <- df %>% 
      mutate(ref_point = ceiling((start + end) / 2),
             start = ref_point - interval,
             end = ref_point + interval)
    return(res)
  }
  if (reference == "start") {
    items <- c("start", "end")
  } else {
    items <- c("end", "start")
  }
  plus <- df %>% 
    filter(strand == "+")
  minus <- df %>% 
    filter(strand == "-")
  plus <- plus %>% 
    mutate(ref_point = get(items[1]),
           start = ref_point - interval,
           end = ref_point + interval)
  minus <- minus %>% 
    mutate(ref_point = get(items[2]),
           start = ref_point - interval,
           end = ref_point + interval)
  res <- bind_rows(plus, minus) %>%
    arrange(chr, start, end)
  return(res)
}
data_start <- list()
data_end <- list()
for (i in 1:length(data)) {
  data_start[[i]] <- get_coords(data[[i]], reference = "start", interval = 2000)
  data_end[[i]] <- get_coords(data[[i]], reference = "end", interval = 2000)
}
names(data_start) <- names(data)
names(data_end) <- names(data)

#test to see if the interval sizes are all okay
for (i in data_end){
  print(unique(i$start - i$ref_point))
  print(unique(i$end - i$ref_point))
}

#Find the overlaps between each of our datasets and the G4-bulge dataset
find_overlaps <- function(df) {
  tmp <- g4bs %>% 
    setkey(chr, start, end) %>% 
    foverlaps(df, .) %>% 
    na.omit()
  return(tmp)
}

for (i in 1:length(data)) {
  data_start[[i]] <- find_overlaps(data_start[[i]])
  data_end[[i]] <- find_overlaps(data_end[[i]])
}

#### OPTIONAL: only keep overlaps that are unique to any given region ####
#First, remove redundant gene segments
# data_start$UTR5 <- NULL
# data_end$UTR3 <- NULL
# # 
# # 
# tmp_list <- list(data_start, data_end)
# for (k in 1:length(tmp_list)) {
#   dummy <- tmp_list[[k]]
#   for (i in 1:length(dummy)) {
#     print(i)
#     tmp <- dummy
#     tmp1 <- tmp[[i]]
#     tmp[[i]] <- NULL
#     tmp <- bind_rows(tmp)
#     tmp_list[[k]][[i]] <- tmp1 %>%
#       filter(!name %in% tmp$name)
#   }
# }
# 
# data_start <- tmp_list[[1]]
# data_end <- tmp_list[[2]]
# 
# for (i in data_start){
#   print(dim(i))
# }
# data_start$second_exo$name %in% data_start$first_exo$name -> kek
# table(kek)
# for (i in 1:length(data_start)) {
#   tmp <- data_start$second_exo$name %in% data_start[[i]]$name
#   print(names(data_start)[i])
#   print(table(tmp))
# }
# 
# #remove data.frames that contain no records
# data_start <- data_start[sapply(data_start, nrow) > 0]
# data_end <- data_end[sapply(data_end, nrow) > 0]
# rm(tmp_list)

#Next, we determine the positions of the overlapping G4BS in relation to the reference point in each dataset
#Critical step that ensures that "-" strand data are oriented the same way as 
# the "+" strand observations.
#### Continue ####
for (i in 1:length(data_start)) {
  data_start[[i]] <- data_start[[i]] %>% 
    mutate(ref_start = ifelse(i.strand == "+", start - ref_point, ref_point - end),
           ref_end = ifelse(i.strand == "+", end - ref_point, ref_point - start))
} 

for (i in 1:length(data_end)) {
  data_end[[i]] <- data_end[[i]] %>% 
    mutate(ref_start = ifelse(i.strand == "+", start - ref_point, ref_point - end),
           ref_end = ifelse(i.strand == "+", end - ref_point, ref_point - start))
}

main_list <- list(data_start, data_end)
names(main_list) <- c("start", "end")


##### OPTIONS MUST BE SPECIFIED PRIOR TO RUNNING!!! #########
#Results directory
setwd("New_G4BS_data/gene_segment_density/")
#Solo or grouped analysis? 
solo_plots <- FALSE
group_plots <- FALSE
#Normalized or non_normalized frequency distribution?? "Freq" = non_normalized, "norm_Freq" = normalized
freq_type <- "norm_Freq"
#ID for data in the combined table?
ID_name <- "RLFS"

#Separate the overlaps based on whether they are on the same or different strands
#set our working list of dataframes
for (z in 1:length(main_list)) {
  wdf <- main_list[[z]]
  wdf_nt <- list()
  wdf_t <- list()
  for (i in 1:length(wdf)) {
    wdf_nt[[i]] <- wdf[[i]] %>% 
      filter(strand == i.strand)
    wdf_t[[i]] <- wdf[[i]] %>% 
      filter(strand != i.strand)
  }
  names(wdf_nt) <- names(wdf)
  names(wdf_t) <- names(wdf)

  #Check if the number of records in the new datasets are equal to the ones in the original
  for (i in 1:length(wdf)) {
    if (nrow(wdf[[i]]) != (nrow(wdf_nt[[i]]) + nrow(wdf_t[[i]]))) {
      print(paste(names(data)[i], "ERROR", sep = " "))
      break
    }
  }

  #Calculate the per nucleotide frequency distribution of G4BS across our intervals
  cl <- makeCluster(3)
  clusterEvalQ(cl, {library(dplyr)})
  registerDoParallel(cl)
  freq_nt <- calculate_frequency(df_list = wdf_nt)
  freq_t <- calculate_frequency(df_list = wdf_t)
  stopCluster(cl)
  names(freq_nt) <- names(wdf_nt)
  names(freq_t) <- names(wdf_t)

  #calculate the normalized frequency for all datasets
  for (i in 1:length(freq_nt)) {
    freq_nt[[i]] <- freq_nt[[i]] %>% 
      mutate(norm_Freq = Freq / sum(Freq))
  }

  for (i in 1:length(freq_t)) {
    freq_t[[i]] <- freq_t[[i]] %>% 
      mutate(norm_Freq = Freq / sum(Freq))
  }

#Plot our results
#Identify the maximum frequency for our plots, so we can use that as the Y axis limits (for better visual comparison between plots

  if (solo_plots) {
    y_limit <- find_max(freq_nt, list2 = freq_t, column = freq_type) * 1.05
    for (i in 1:length(freq_nt)) {
      p1 <- plot_res(df = freq_nt[[i]], x = "position", y = freq_type)
      p2 <- plot_res(df = freq_t[[i]], x = "position", y = freq_type)
      name1 <- paste(ID_name, names(freq_nt)[i], names(main_list)[z], "nonTemplate", sep = "_")
      name2 <- paste(ID_name, names(freq_nt)[i], names(main_list)[z], "Template", sep = "_")
      
      ggsave(paste0(name1, freq_type, ".png"), p1, dpi = 300, width = 9, height = 7)
      ggsave(paste0(name2, freq_type, ".png"), p2, dpi = 300, width = 9, height = 7)
      
    }
  }
  #Combined analysis of G4 and G4BS data
  #Lots of rerunning required :S
  for (i in 1:length(freq_nt)) {
  freq_nt[[i]] <- freq_nt[[i]] %>% 
    mutate(G4_type = ID_name,
           strand_type = "NT")
  }
  for (i in 1:length(freq_t)) {
    freq_t[[i]] <- freq_t[[i]] %>% 
      mutate(G4_type = ID_name,
             strand_type = "T")
  }
  main <- list()
  for (i in 1:length(freq_nt)) {
    tmp <- list(freq_nt[[i]], freq_t[[i]]) %>% 
      bind_rows()
    main[[i]] <- tmp
    if (names(freq_nt)[i] == names(freq_t)[i]) {
      names(main)[i] <- names(freq_nt)[i] 
    } else {
      print("Something's wrong, I can feel it!")
      break
    }
  }
  tmp_name <- paste(ID_name, names(main_list)[z], sep = "_")
  if (!tmp_name %in% names(final)) {
    counter <- length(final) + 1
    final[[counter]] <- main
    names(final)[counter] <- tmp_name 
  }
}
#Afterwards, there are cycles of unpacking the lists by creating grouping variables then binding the rows of different tables
#First, we will create a "start" and "end" list, each containing the same number of tables for each gene segment, while also creating
#the "ID" variable that will allow for the combining of these two lists
final2 <- list()
for (i in 1:2) {
  kek <- c("start", "end")
  tmp_list <- list()
  for (k in 1:length(final[[i]])) { 
    tmp <- bind_rows(final[[i]][[k]], final[[i + 2]][[k]], final[[i + 4]][[k]], 
                     final[[i + 6]][[k]], final[[i + 8]][[k]], final[[i + 10]][[k]],
                     final[[i + 12]][[k]]) %>% 
      mutate(ID = kek[i])
    tmp_list[[k]] <- tmp
  }
  names(tmp_list) <- names(final[[1]])
  final2[[i]] <- tmp_list
}
names(final2) <- c("start", "end")
#We combine the corresponding dataframes in the "start" and "end" lists
tmp_list <- list()
for (i in 1:length(final2$start)) { 
  tmp <- bind_rows(final2$start[[i]], final2$end[[i]])
  tmp_list[[i]] <- tmp
}
names(tmp_list) <- names(final2$start)
#We create the "segment" variable to be able to combine the dataframes for different gene segments
for (i in 1:length(tmp_list)) {
  tmp_list[[i]] <- tmp_list[[i]] %>% 
    mutate(segment = names(tmp_list)[i])
}
#We finally combine the dataframes into a single table
#We get a single dataframe that has 4 grouping variables ("G4_type", "strand_type", "ID", and "segment")

final_df <- bind_rows(tmp_list)
# write.table(final_df, "All_G4motifs_dist_MAIN_2kb.1129.txt", sep = "\t", col.names = T, row.names = F)
# write.table(final_df, "ensembl_test.txt", sep = "\t", col.names = T, row.names = F)

#Let's plot the different combinations of data
#Incoming giant loop
freq_type <- "norm_Freq"
y_limit <- find_max(final_df, column = freq_type) * 1.05

setwd("/home/csaba/work/Projects/G4BS_new_results/New_G4BS_data/gene_segment_distribution/Distribution_plots/")

for (i in 1:length(unique(final_df$segment))) {
  cur_seg <- unique(final_df$segment)[i]
  for (k in 1:length(unique(final_df$ID))) {
    cur_ID <- unique(final_df$ID)[k]
    non_template <- final_df %>% 
      filter(segment == cur_seg & ID == cur_ID & strand_type == "NT")
    template <- final_df %>% 
      filter(segment == cur_seg & ID == cur_ID & strand_type == "T")
    p1 <- plot_res(non_template, "position", freq_type, "G4_type", limits = 250)
    # name1 <- paste(cur_seg, cur_ID, "nonTemplate", freq_type, sep = "_")
    p2 <- plot_res(template, "position", freq_type, "G4_type", limits = 250)
    # name2 <- paste(cur_seg, cur_ID, "Template", freq_type, sep = "_")
    name <- paste(cur_seg, cur_ID, freq_type, ".png", sep = "_")
    g <- arrangeGrob(p1, p2, ncol = 2)
    ggsave(name, g, width = 31, height = 13, dpi = 300)
  }
}

#### FUNCTIONS ####
calculate_frequency <- function(df_list) {
  res <- list()
  res <- foreach(i = 1:length(df_list)) %dopar% {
    list1 <- list()
    tmp <- df_list[[i]] 
    for (a in 1:nrow(tmp)){
      list1[[a]] <- seq(tmp$ref_start[a], tmp$ref_end[a])
    }
    
    tmp1 <- unlist(list1) %>% 
      table() %>% 
      data.frame() %>% 
      select(position = 1, Freq) %>% 
      mutate(position = as.numeric(levels(position)))
    
    tmp1
  }
  return(res)
}

find_max <- function(list1, list2 = NA, column) {
  tmp1 <- bind_rows(list1)
  if (missing(list2)) {
    highest <- tmp1 %>% 
      select(all_of(column)) %>% 
      max()
  } else {
    tmp2 <- bind_rows(list2)
    highest <- bind_rows(tmp1, tmp2) %>% 
      select(all_of(column)) %>% 
      max()
  }
  return(highest)
}

plot_solo_res <- function(df, x, y) {
  p1 <- ggplot(df, aes(x = get(x), y = get(y))) + geom_point(size = 1) + 
    theme_gray() +
    theme(text = element_text(size = 30)) +
    labs(x = ("Position relative to reference point"), y = 'Frequency') + 
    xlim(-250, 250) +
    ylim(0, y_limit) + 
    guides(colour = guide_legend(override.aes = list(size=7))) +
    geom_vline(xintercept = c(0), size = 0.5)
  
  return(p1)
}

plot_res <- function(df, x, y, group, limits = 250) {
  if (missing(group)) {
    p1 <- ggplot(df, aes(x = get(x), y = get(y))) + 
      geom_point(size = 1) +
      theme_gray() +
      theme(text = element_text(size = 30)) +
      labs(x = ("Position relative to reference point"), y = 'Frequency') + 
      xlim(-limits, limits) +
      ylim(0, y_limit) + 
      guides(colour = guide_legend(override.aes = list(size=7))) #+
      # geom_vline(xintercept = c(0), size = 0.5)
  } else {
    p1 <- ggplot(df, aes(x = get(x), y = get(y))) + 
      geom_point(data = df, size = 1, aes(color = get(group))) + 
      scale_color_manual(values = c("red","green4", "blue", "orange2", "black", "purple")) +
      labs(x = ("Position relative to reference point"), y = 'Frequency',
           color = "Type of G4") +
      theme_gray() +
      theme(text = element_text(size = 30)) +
      xlim(-limits, limits) +
      ylim(0, y_limit) + 
      guides(colour = guide_legend(override.aes = list(size=7))) #+
      # geom_vline(xintercept = c(0), size = 0.5)
  }
  return(p1)
}
 