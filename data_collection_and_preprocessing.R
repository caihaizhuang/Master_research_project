library("pracma")
library(parallel)
library(MASS)
library("ape")
library(Biostrings)
library(seqinr)
library(stringr)
library(dplyr)
###############
setwd("/Volumes/hard_drive_Lucas/haizhuang/research_project_data/target_protein_seq_filter/")
dir_list <- list.files()
head(dir_list,10)
species_list <- str_split(dir_list, "\\.")
head(species_list, 10)
species_list1 <- matrix(NA, 28452, 1)
species_list1 <- as.data.frame(species_list1)
for (i in 1:length (species_list)) {
    species_list1[i,1] <- species_list[[i]][[1]]
}
colnames(species_list1)[1] <- "strain_names"
species_list1$species_names <- "NA"


for (i in 1:length(species_list1$species_names)) { #get the name of different species from strain names
    ans <- str_split(species_list1[i,1], "_")
    if (ans[[1]][[1]]=="") {
        species_list1$species_names[i] <- paste(ans[[1]][[2]], "_", ans[[1]][[3]], sep = "")
    } else {
        species_list1$species_names[i] <- paste(ans[[1]][[1]], "_", ans[[1]][[2]], sep = "")
    }
    
}

species_list2 <- species_list1$species_names[!duplicated(species_list1$species_names)]  #remove the duplicated names
species_list2 <- as.data.frame(species_list2)
species_list2$fasta.count <- NA
colnames(species_list2)[1] <- "species_names"

for (i in 1 :length(species_list2$fasta.count)) { #count how many fasta is there for each species
    species_list2$fasta.count[i] <- length(grep(species_list2$species_names[i], dir_list))
}

species_list1 <- na.omit(species_list1)
dir_list <- as.character(dir_list)
species_list1$fasta_name <- dir_list  # adding the fasta file name to each strain 
class(species_list1$fasta_name)

#calculate the length of protein in all fasta(only the first the AA seq)
species_list1$protein_length <- NA
for (i in 1:length(dir_list)) {
    ans <- readBStringSet(as.character(dir_list[i]))
    ans <- ans[1]
    species_list1$protein_length[i] <-str_length(ans[1])
}
hist(species_list1$protein_length, 100)

###################Test: data extraction of phylum level(6-phylum)#################
#1.read phylum information#
setwd("/Users/caihaizhuang/Desktop/")
Phylum_info <- read.csv("Phylum_genus.csv")
Phylum_info$genus[Phylum_info$genus==""] <- NA
Phylum_info <- na.omit(Phylum_info)

#2.make a fasta count of each Phylum
species_list1$genus_index <- NA
for (i in 1:length(species_list1$species_names)) { #make a genus index for selection
    species_list1$genus_index[i]<- paste("_",str_split(species_list1$species_names[i], "_")[[1]][1], "_", sep = "")
}
colnames(Phylum_info)[5] <- "fasta_count"
Phylum_info$fasta_count <- NA
for (i in 1:length(Phylum_info$genus)) {
    Phylum_info$fasta_count[i] <- length(grep(paste(Phylum_info$genus[i], "_", sep = ""), species_list1$species_names, ignore.case = TRUE))
}
###################Test: data extraction of phylum level(4-phylum) base on sunjae's phylum info#################
###build genus_index
species_list1$genus_index <- NA
for (i in 1:length(species_list1$species_names)) { #make a genus index for selection
    species_list1$genus_index[i]<- paste("_",str_split(species_list1$species_names[i], "_")[[1]][1], "_", sep = "")
}
###
setwd("/Users/caihaizhuang/Desktop/project_data/")
phylum_info <- read.csv("Bacteria_of_4phylum.csv")
phylum_info$"genus_name" <- NA
for (i in 1:length(phylum_info$species)) {
    phylum_info$genus_name[i] <- paste("_",str_split(phylum_info$species[i], " ")[[1]][1], "_", sep = "")
}
phylum_info <- phylum_info[!duplicated(phylum_info$genus_name),]

phylum_info$fasta_count <- NA
for (i in 1:length(phylum_info$genus)) {
    phylum_info$fasta_count[i] <- length(grep(phylum_info$genus_name[i], species_list1$genus_index, ignore.case = TRUE))
}
phylum_info<- phylum_info[!phylum_info$fasta_count==0,]

phylum_info_1 <- phylum_info[rep(seq(nrow(phylum_info)), phylum_info$fasta_count),]
phylum_info_1$fasta_dir <- NA


a <- c() #build a character of fasta file and append value to it, and merge it with "phylum_info_1"
as.character(a)
for (i in 1:length(phylum_info$phylum)) {
    ans <- species_list1[grep(phylum_info$genus_name[i], species_list1$genus_index, ignore.case = TRUE), c(3)]
    a <- append(a, ans, after = length(a))
}
phylum_info_1$fasta_dir <- a

# test the length of AA seq
setwd("/Volumes/hard_drive_Lucas/haizhuang/research_project_data/target_protein_seq_filter/")
phylum_info_1$protein_length <- NA
for (i in 1:length(phylum_info_1$fasta_dir)) {
    ans <- readBStringSet(as.character(phylum_info_1$fasta_dir[i]))
    ans <- ans[1]
    phylum_info_1$protein_length[i] <-str_length(ans[1])
}
setwd("/Users/caihaizhuang/Desktop/")

# eliminate all the protein with a length shorter than 500
phylum_info_1 <- phylum_info_1[phylum_info_1$protein_length>=500, ]

# Cut all AA seq into same length ==> 500 AA seq here!
setwd("/Volumes/hard_drive_Lucas/haizhuang/research_project_data/target_protein_seq_filter/")
phylum_info_1$protein_seq <- NA
phylum_info_1$length_after_cut <- NA
for (i in 1:length(phylum_info_1$protein_length)) {
    ans <- readBStringSet(as.character(phylum_info_1$fasta_dir[i]))
    ans <- ans[1]
    AA_seq <- str_sub(ans, 1, 500)
    LAC <- str_length(AA_seq)
    phylum_info_1$protein_seq[i] <- AA_seq
    phylum_info_1$length_after_cut[i] <- LAC
}
setwd("/Users/caihaizhuang/Desktop/")

#randomize the order of fasta file each phylum
#####Firmicutes#####
Firmicutes_info <- phylum_info_1[phylum_info_1$phylum=="Firmicutes",]
set.seed(1)
size = length(Firmicutes_info$phylum)
Firmicutes_info <- sample_n(Firmicutes_info, size = size, replace = FALSE, weight = NULL, .env = NULL) #randomize the order of data
Firmicutes_info <- Firmicutes_info[1:200,]

#####Bacteroidetes#####
Bacteroidetes_info <- phylum_info_1[phylum_info_1$phylum=="Bacteroidetes",]
set.seed(1)
size = length(Bacteroidetes_info$phylum)
Bacteroidetes_info <- sample_n(Bacteroidetes_info, size = size, replace = FALSE, weight = NULL, .env = NULL)
Bacteroidetes_info <- Bacteroidetes_info[1:200,]

#####Proteobacteria#####
Proteobacteria_info <- phylum_info_1[phylum_info_1$phylum=="Proteobacteria",]
set.seed(1)
size = length(Proteobacteria_info$phylum)
Proteobacteria_info <- sample_n(Proteobacteria_info, size = size, replace = FALSE, weight = NULL, .env = NULL)
Proteobacteria_info <- Proteobacteria_info[1:200,]

#####Actinobacteria#####
Actinobacteria_info <- phylum_info_1[phylum_info_1$phylum=="Actinobacteria",]
set.seed(1)
size = length(Actinobacteria_info$phylum)
Actinobacteria_info <- sample_n(Actinobacteria_info, size = size, replace = FALSE, weight = NULL, .env = NULL)
Actinobacteria_info <- Actinobacteria_info[1:200,]
#write.fasta(as.list(Actinobacteria_info$protein_seq), Actinobacteria_info$fasta_dir, "Actinobacteria_samples.fasta", nbchar = 500)

#####collection####
collection_phylum <- rbind(Proteobacteria_info, Actinobacteria_info,Bacteroidetes_info, Firmicutes_info) #merge four dataframe 
set.seed(1)
size = length(collection_phylum$phylum)
collection_phylum <- sample_n(collection_phylum, size = size, replace = FALSE, weight = NULL, .env = NULL) #randomize the data collection
colnames(meta_phylum_info)[1] <- "phylum"
ML_data <- left_join(collection_phylum, meta_phylum_info, by = "phylum") #join the info table and extract only the sequences and labels
ML_data_P <- ML_data[,c(7,10)]
colnames(ML_data_P)[2] <- "MLID"
#output the machine learning data
setwd("/Users/caihaizhuang/Desktop/")
write.csv(ML_data_P ,"ML_data_test_4phylum.csv")
write.fasta(as.list(ML_data_P$protein_seq), ML_data_P$MLID, "collection_samples.fasta", nbchar = 500)


#build a metadata table
meta_phylum_info <- phylum_info$phylum[!duplicated(phylum_info$phylum)]
meta_phylum_info <- as.data.frame(meta_phylum_info)
colnames(meta_phylum_info)[1] <- "phylum_name"
meta_phylum_info$fasta_count <- NA
for (i in 1:4) {
    meta_phylum_info$fasta_count[i] <- sum(phylum_info$fasta_count[grep(meta_phylum_info$phylum_name[i], phylum_info$phylum, ignore.case = TRUE)])
}
meta_phylum_info$index <- c(0:3) #give each phylum an index

setwd("/Users/caihaizhuang/Desktop/")
##########################################test end###########################################

########################## data extract and method test with 4 species###########################
#build index for 4 example species
setwd("/Volumes/hard_drive_Lucas/haizhuang/research_project_data/target_protein_seq_filter/")
index_1 <- grep("Staphylococcus_aureus", dir_list)[1:2000]
index_2 <- grep("Streptococcus_pneumoniae", dir_list)[1:2000]
index_3 <- grep("Salmonella_enterica", dir_list)[1:2000]
index_4 <- grep("Escherichia_coli", dir_list)[1:2000]
species_list_3 <- species_list1[c(index_1, index_2, index_3, index_4),]
species_list_3$protein_length <- NA

species_list_test <- dir_list[c(index_1, index_2, index_3, index_4)]

# test the length of AA seq
for (i in 1:length(species_list_test)) {
    ans <- readBStringSet(as.character(species_list_test[i]))
    ans <- ans[1]
    species_list_3$protein_length[i] <-str_length(ans[1])
}

# Cut all AA seq into same length ==> 500 AA seq here!

species_list_3$protein_seq <- NA
species_list_3$length_after_cut <- NA
for (i in 1:length(species_list_test)) {
    ans <- readBStringSet(species_list_test[i])
    ans <- ans[1]
    AA_seq <- str_sub(ans, 1, 500)
    LAC <- str_length(AA_seq)
    species_list_3$protein_seq[i] <- AA_seq
    species_list_3$length_after_cut[i] <- LAC
}

#make label for 4 example species
species_list_4example <- species_list2[species_list2[,2]>1000,]
species_list_4example$MLID <- NA
k <- 0
for (i in 1:4) {
    species_list_4example$MLID[i] <- k
    k <- k+1
}

library(dplyr)
ML_data <- left_join(species_list_4example, species_list_3, by = "species_names") #join the info table and extract only the sequences and labels
ML_data <- ML_data[,c(3,8)]

#randomize the AA seq
set.seed(1)
ML_data_S <- sample_n(ML_data, size=8000, replace = FALSE, weight = NULL, .env = NULL)

#output the machine learning data
setwd("/Users/caihaizhuang/Desktop/")
write.csv(ML_data ,"ML_data_test_4example.csv")
###################################end of the species level test###############################

##################################validation test phylum#################
#subset AA seq from the same group(Firmicutes_info) to 4 group 
Firmicutes_info_1 <- Firmicutes_info
# Cut all AA seq into same length ==> 500 AA seq here!

for (i in 1:50) {
    ans <- Firmicutes_info_1$protein_seq[i]
    Firmicutes_info_1$protein_seq[i] <- str_sub(ans, 1, 440)
    Firmicutes_info_1$MLID[i] <- 0
    
}
for (i in 51:100) {
    ans <- Firmicutes_info_1$protein_seq[i]
    Firmicutes_info_1$protein_seq[i] <- str_sub(ans, 21, 460)
    Firmicutes_info_1$MLID[i] <- 1
}
for (i in 101:150) {
    ans <- Firmicutes_info_1$protein_seq[i]
    Firmicutes_info_1$protein_seq[i] <- str_sub(ans, 41, 480)
    Firmicutes_info_1$MLID[i] <- 2
}
for (i in 151:200) {
    ans <- Firmicutes_info_1$protein_seq[i]
    Firmicutes_info_1$protein_seq[i] <- str_sub(ans, 61, 500)
    Firmicutes_info_1$MLID[i] <- 3
}
Firmicutes_info_1$length_after_cut <- "440"
ML_data_V <- Firmicutes_info_1[,c(7,9)]



#randomize the AA seq
set.seed(1)
ML_data_V <- sample_n(ML_data_V, size=200, replace = FALSE, weight = NULL, .env = NULL)

#output the machine learning data
setwd("/Users/caihaizhuang/Desktop/")
write.csv(ML_data_V ,"ML_data_test_Validation.csv")



