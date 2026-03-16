library(googlesheets4)
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
