################################################
##    Created by ZAINAB SIDDIQUI              ##
##    On 11/04/2026                           ##
##    Test of SRA Accession embedding in R    ##
##      & Draft Preprocessing Pipeline        ##
################################################


library(data.table)

# -------------------- Loading SRA ---------------------

# test on one isolated accession
test_id <- read.delim("../breastcancer_kmers/data/abc_ids.txt",header = FALSE)[1,]

# leads to data directory; storage for SRA prefetch
## LIMITATION: relative paths are weak/not recommended in case of reorganization
data_dir <- file.path("..",       #../breastcancer_kmers/data/prefetched
                      "breastcancer_kmers",
                      "data",
                      "prefetched") 

message("Downloading SRA file...")
system(paste("prefetch", test_id,"-O",data_dir))

accession_id  <- list.files(data_dir) # uses ID that's being used to label download
accession_dir <- file.path(data_dir,accession_id)

fetched_file <- file.path(accession_dir,list.files(accession_dir))

# extract with fastq-dump
message("Extracting FASTQ data...")
system(paste("fasterq-dump --split-files", fetched_file,"-O",accession_dir))
# delete R2
message("Deleting reverse read file...")
system(paste("rm",file.path(accession_dir,
                            list.files(accession_dir,pattern="_2\\.fastq$")))) 

# clean dumped FASTQ
uncleaned_fastq_file <- file.path(accession_dir, 
                                 list.files(accession_dir,
                                            pattern = "_1\\.fastq$")
                                 )

cleaned_fastq_file <- file.path(accession_dir,paste0(accession_id,"_seqonly.txt"))
#dir.create(cleaned_fastq_dir, recursive = TRUE) # make sure path exists

# take the 2nd row of every 4 to isolate sequences
message("Extracting sequences from FASTQ dump...")
system(paste("awk 'NR % 4 == 2'", uncleaned_fastq_file,">",cleaned_fastq_file))
# 9.2GB .fastq --> 3.25GB .txt

# delete original file
message("Deleting original fastq file...")
system(paste("rm",uncleaned_fastq_file)) 

# random sampling of sequence-lines
n = 1e6 # num lines to sample
sampled_fastq_file <- file.path(accession_dir,
                                paste0(accession_id,"_sample",
                                       gsub("\\+0", "", 
                                            format(n, scientific = TRUE)),
                                       ".txt"))

lines <- fread(cleaned_fastq_file, header = FALSE)[[1]]
sampled_lines <- sample(lines,n)

writeLines(sampled_lines, sampled_fastq_file)
# 3.25GB --> 1.51 GB











