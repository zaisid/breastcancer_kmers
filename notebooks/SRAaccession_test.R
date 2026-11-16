################################################
##    Created by ZAINAB SIDDIQUI              ##
##    On 11/04/2026                           ##
##    Test of SRA Accession embedding in R    ##
##      & Draft Preprocessing Pipeline        ##
################################################


library(data.table)
library(Biostrings)

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


# -------------------- Cleaning SRA ---------------------

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
rm(lines) 

writeLines(sampled_lines, sampled_fastq_file)
# 3.25GB --> 1.51 GB


# -------------------- K-merizing ---------------------

k <- 2
chunk_size <- 50000  # lighten memory load

# Initialize aggregate counts
kmer_sums <- numeric(4^k)
names(kmer_sums) <- oligonucleotideFrequency(
    DNAString("ACGT"),
    width = k
) |> names()

# Open sampled sequence file
con <- file(sampled_fastq_file, open = "r")

chunk_num <- 0

repeat {
    # Read next chunk
    seq_chunk <- readLines(con, n = chunk_size)
    
    # Stop when file is exhausted
    if (length(seq_chunk) == 0) break
    
    chunk_num <- chunk_num + 1
    
    # Convert to DNA sequences
    dna_chunk <- DNAStringSet(seq_chunk)
    
    # Aggregate k-mer counts across chunk
    chunk_counts <- oligonucleotideFrequency(
        dna_chunk,
        width = k,
        simplify.as = "collapsed"
    )
    
    # Accumulate counts
    kmer_sums[names(chunk_counts)] <-
        kmer_sums[names(chunk_counts)] + chunk_counts
    
    # Progress indicator
    if (chunk_num %% 10 == 0) {
        message("Processed ", chunk_num * chunk_size,
                " sequences")
    }
    
    rm(seq_chunk, dna_chunk, chunk_counts)
}

close(con)

# Normalize into frequencies 
kmer_freq <- kmer_sums / sum(kmer_sums)

kmer_freq



