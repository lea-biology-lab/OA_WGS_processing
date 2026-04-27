library(googlesheets4)
library(readxl)
library(dplyr)
library(lubridate)
#code from Audrey Arner

######################
# wgs that I fixed thru concordance
######################
fixed_tids=read_sheet("https://docs.google.com/spreadsheets/d/1JF7aLZbsl_JzM9DW8YbcmY55BN2oShxJTPpIVqUrFyw/edit?gid=1855266807#gid=1855266807", sheet = 1)
remove=read_sheet("https://docs.google.com/spreadsheets/d/1JF7aLZbsl_JzM9DW8YbcmY55BN2oShxJTPpIVqUrFyw/edit?gid=65024686#gid=65024686", sheet = 2)

fixed_tids=subset(fixed_tids, data_type=="wgs")
remove=subset(remove, data_type=="wgs")

#add "tid_" to column so will match with genotype_ids
fixed_tids$original_tid <- paste0("tid_", fixed_tids$original_tid)
fixed_tids$corrected_tid <- paste0("tid_", fixed_tids$corrected_tid)
remove$original_tid <- paste0("tid_", remove$original_tid)

#final set of WGS
genotype_ids=read.delim("~/Library/CloudStorage/Box-Box/Audrey/Lab/eqtl_analyses_new/oa_1_phased_subset.fam", header = F)

#look at how many we need to fix
to_fix=subset(fixed_tids, original_tid %in% genotype_ids$V2) 
nrow(to_fix) #11

#here we will want to look at which "fixes" are already available as values
already_in_data=subset(to_fix, corrected_tid %in% genotype_ids$V2)
nrow(already_in_data) #8
#make sure that all genotype_ids are unique
any(duplicated(genotype_ids$V2)) #FALSE

#make a new column of fixed genotype_ids
genotype_ids$corrected_TID=NA

genotype_ids$corrected_TID <- ifelse(
  !is.na(match(genotype_ids$V2, to_fix$original_tid)),
  to_fix$corrected_tid[match(genotype_ids$V2, to_fix$original_tid)],
  genotype_ids$V2
)

any(duplicated(genotype_ids$corrected_TID)) #TRUE
dup_tids <- genotype_ids[genotype_ids$corrected_TID %in% genotype_ids$corrected_TID[duplicated(genotype_ids$corrected_TID)], ]

#we are going to keep the tids with the highest coverage
coverage=read_sheet("https://docs.google.com/spreadsheets/d/1Xb9mUALNtjmLVQRILqbPuijbV-XhmdvCiQwGIwOWrgs/edit?gid=1708865777#gid=1708865777", sheet = 5)
coverage$tid <- paste0("tid_", coverage$tid)

dup_tids=merge(dup_tids, coverage, by.x="V2", by.y="tid")

#identify the TIDs that are the lower coverage out of the duplicates
keep <- with(dup_tids, ave(coverage, corrected_TID, FUN = function(x) x == max(x)))

# subset to keep only those rows so know which ones to delete
dup_tids_delete <- dup_tids[keep == 0, ]

#here is where we actually get rid of the rows we don't want (lower coverage)
genotype_ids <- genotype_ids[ !genotype_ids$V2 %in% dup_tids_delete$V2, ]

#rename column
names(genotype_ids)[2]="sequenced_tid"

genotype_ids=genotype_ids[c("sequenced_tid", "corrected_TID")]

write.table(genotype_ids, "~/Library/CloudStorage/Box-Box/Audrey/Lab/WGS/wgs_samples_touse_25Feb26.txt", row.names=F, quote=F, sep = "\t")

######################
#look at ethnicity groups and ages
######################
medical=read.csv('~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/Malaysia/data/oahelp/data/medical.csv')
personal=read.csv('~/Library/CloudStorage/Box-Box/Lea Lab/Audrey_Arner/Lab/Malaysia/data/oahelp/data/personal_information.csv')

#get ethnicity for meta data
compile_ethnic <- function(personal_info = NULL, med = NULL){
  if(is.null(personal_info)) stop("ERROR: Must include valid dataframe as population register")
  
  # add medical visit village for contextualizing info
  personal_info$interview_location_med <- med$interview_location_med[match(personal_info$rid, med$rid)]
  personal_info$other_ethnicity[which(personal_info$birth_place_state_other == "BugIS INDONESIA")] <- "Bugis"
  
  # First, deal with "Other" category by standardizing to fewer number of categories
  personal_info$other_ethnicity <- tolower(personal_info$other_ethnicity)
  
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("bidayuh"))]<-'Bidayuh'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("bugis"))]<-'Bugis'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("dayak, orang melayu"))]<-'Dayak_Malay'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("dusun"))]<-'Dusun'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("iban"))]<-'Iban'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("father indo", "indo", "indonesia", "indonesian"))]<-'Indonesian'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("jah hut", "jahut"))]<-'Jahut'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("jahut dan cewang"))]<-'Jahut_Cewang'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("kentak", "kintaq"))]<-'Kintaq'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("lano"))]<-'Lanoh'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("mah meri"))]<-'MahMeri'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("melayu", "orang melayu"))]<-'Malay'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("murut"))]<-'Murut'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("kensiu"))]<-'Kensiu'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("kensiu, melayu (father)"))]<-'Kensiu_Malay'
  personal_info$other_ethnicity[which(personal_info$other_ethnicity %in% c("semelai"))]<-'Semelai'
  
  # special case of Jakun
  personal_info$ethnicity___2[which(personal_info$ethnicity___2 == 0 & personal_info$interview_location_med == 25)]<-1  # this is Kg Petoh, so almost certainly jakun
  
  # this case marked only temiar just needs to be moved
  personal_info$ethnicity___9[which(personal_info$other_ethnicity %in% c("temiar"))] <- 1
  personal_info$ethnicity___99[which(personal_info$other_ethnicity %in% c("temiar"))] <- 0
  
  ## Extract numbers associated with ethnicities 
  # Get columns that match the pattern "ethnicity___[number]"
  pattern_cols <- grep("ethnicity___\\d+", names(personal_info), value = TRUE)
  
  # Extract numbers from column names
  numbers <- as.numeric(sub(".*___(\\d+)", "\\1", pattern_cols))
  
  # Create list to store results
  result_list <- vector("list", nrow(personal_info))
  
  # For each row
  for (i in 1:nrow(personal_info)) {
    
    # Get indices where value is 1
    ones_indices <- which(personal_info[i, pattern_cols] == 1)
    
    # Store corresponding numbers in result list
    result_list[[i]] <- numbers[ones_indices]
  }
  
  ## Replace numbers with the ethnic group names
  result_list <- lapply(result_list, FUN=function(x){
    case_match(x,
               0 ~ "Batek",
               1 ~ "Jehai",
               2 ~ "Jakun",
               3 ~ "MahMeri",
               4 ~ "Mendriq",
               5 ~ "OrangKuala",
               6 ~ "OrangSeletar",
               7 ~ "Semai",
               8 ~ "SemaqBeri",
               9 ~ "Temiar",
               10 ~ "Temuan",
               99 ~ "Other")
  })
  
  ## Replace "Other" designations in list
  for(i in 1:length(result_list)){
    if(sum(result_list[[i]] == "Other") > 0){
      result_list[[i]] <- case_match(result_list[[i]], 
                                     "Other" ~ personal_info$other_ethnicity[[i]])
    }
    
    result_list[[i]][which(is.na(result_list[[i]]))] <- "Unknown"
  }
  
  
  
  ## Assign 'Multiple' to ethnicity if more than one ethnicity is indicated
  personal_info$number_ethnics <- unlist(lapply(result_list, length)) 
  
  ## Make unknown for anyone without ethnicity information
  result_list[which(personal_info$number_ethnics == 0)] <- "Unknown"
  
  
  ## Create a general "ethnicity" column for use
  personal_info$ethnicity <- unlist(lapply(result_list, paste, collapse="_"))
  
  # make broader groups
  personal_info$ethnicity_groups<-'Other/Mixed'
  personal_info$ethnicity_groups[which(personal_info$ethnicity %in% c('Batek','Jehai','Mendriq', "Kensiu", "Lanoh", "Kintaq"))]<-'Semang'
  personal_info$ethnicity_groups[which(personal_info$ethnicity %in% c('Temuan','Jakun','OrangSeletar','OrangKuala','Semelai'))]<-'Proto_Malay'
  personal_info$ethnicity_groups[which(personal_info$ethnicity %in% c('Temiar','Semai','SemaqBeri','MahMeri',"Jahut"))]<-'Senoi'
  
  return(personal_info)
}
meta=compile_ethnic(personal, medical)

medical = medical[!is.na(medical$tid), ]

oa_data <- left_join(medical, meta, by = "rid", suffix = c(".medical", ".personal")) #join data

oa_data <- oa_data %>%
  mutate(age = as.numeric(ymd(med_date) - ymd(date_of_birth)) / 365.25)

oa_data$tid <- paste0("tid_", oa_data$tid)

oa_data_wgs=oa_data[oa_data$tid %in% genotype_ids$corrected_TID,]

table(oa_data_wgs$ethnicity)

#look into individuals with ethnicity "unknown"
oa_data_19=oa_data[oa_data$interview_location.medical==19,c("rid", "tid", "med_name", "ethnicity", "father_name")]
oa_data_wgs_see=oa_data_wgs[oa_data_wgs$ethnicity=="Unknown",] 
#both individuals are from Tementong and are daughters of Rosli. All other children of Rosli are Batek.

oa_data_wgs$ethnicity=ifelse(oa_data_wgs$ethnicity=="Unknown", "Batek", oa_data_wgs$ethnicity)

#remove individuals who aren't OA
oa_data_wgs=oa_data_wgs[!oa_data_wgs$ethnicity %in% c("Bugis", "Dayak_Malay", "Iban", "Malay", "Murut", "Unknown_Indonesian"),]

#remove individuals who are younger than 17.5
oa_data_wgs=oa_data_wgs[oa_data_wgs$age > 17.5,]

genotype_ids2=genotype_ids[genotype_ids$corrected_TID %in% oa_data_wgs$tid,]

write.table(genotype_ids2, "~/Library/CloudStorage/Box-Box/Audrey/Lab/WGS/wgs_samples_touse_22Apr26.txt", row.names=F, quote=F, sep = "\t")

#list of all high coverage data to use
genotype_ids2=read.delim("~/Library/CloudStorage/Box-Box/Audrey/Lab/WGS/wgs_samples_touse_22Apr26.txt")

genotype_ids2=merge(genotype_ids2, coverage, by.x="sequenced_tid", by.y="tid")

genotype_ids_high=genotype_ids2[genotype_ids2$coverage > 15,]

write.table(genotype_ids_high, "~/Library/CloudStorage/Box-Box/Audrey/Lab/WGS/wgs_highcov_samples_touse_22Apr26.txt", row.names=F, quote=F, sep = "\t")
